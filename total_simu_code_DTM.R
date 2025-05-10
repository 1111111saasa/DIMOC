rm(list=ls())
gc()
library(doParallel)
library(foreach)
library(lars)
library(MASS)


##########################################
######## Settings: Parallel ##############
#######################################
ns =8       ###simulation times
no_cores <- detectCores()
cl <- makeCluster(no_cores-8)
registerDoParallel(cl)

##########################################
######## Settings: Parameters ############
##########################################
n =400  ##sample size
B=30    ##refinement times
p =4000  ##dimension
ratio=0.5 ##splitting ratio
pi10 = 0.03  ##portion of alpha!=0,beta=0
pi01 = 0.01
pi00 = 1-pi10-pi01*2
pi = c(pi00,pi10,pi01,1-sum(pi00,pi10,pi01))
nonzero = (sum(pi[1:3])*p+1):p
q = 0.1 ##FDR level
tau = 3 ##signal strength
corr=0.5  ##dependence

#########  The input for X should be an array of n*1, 
#########   M for n*p matrix, Y for n*1 array, n for the sample size.
#########   p for dimension, pi for the portion of pi00, pi10, pi01, pi11.
#########   nonzero for the index of pi11, q for the target fdr level, B
#########   for refinement times.



Process <- function(n,p,pi,nonzero,q,tau,B,correlation=corr){
  #####Model: Mk=ak*X+ek;Y=r*X+bk*Mk+epsk,k=1,...,p########
  
  
  
  #------------Functions for generating data--------------------------#
  
  ####################################################################
  ############################generating data#########################
  ######Model: Mk=ak*X+ek;Y=r*X+sum_{k=1}^p bk*Mk+eps,k=1,...,p#######
  
  ModelMultiReg <- function(n,p,pi,matrix,tau,index){
    
    ####generate X####
    # X = rnorm(n) #X~N(0,1)
    X = rbinom(n,1,0.2) #X~Ber(0.2)
    
    ####generate alpha and beta######
    m = pi*p
    coeff = data.frame(alpha=rep(0,p),beta=rep(0,p))
    
    beta_Ind=c()
    non_beta=(sum(pi[1:2])*p+1):p
    for(i in 1:(pi[3]*p) ){
      beta_Ind[2*i-1]=(1:p)[(sum(pi[1:2])*p)+i]
    }
    for(i in 1:(pi[4]*p) ){
      beta_Ind[2*i]=(1:p)[(p*sum(pi[1:3]))+i]
    }
    Ind=1:p
    Ind[(sum(pi[1:2])*p+1):p]=beta_Ind
    temp_frame=data.frame(Ind=Ind, temp_ind=1:p)
    temp_frame=split(temp_frame,ceiling(temp_frame$temp_ind/2))
    Ind=lapply(temp_frame, function(x){x$Ind})
    Ind=Ind[index]
    Index1=c()
    for(i in 1:(p/2)){Index1=c(Index1,Ind[[i]])}
    
    aInd = (m[1]+1):(m[1]+m[2])
    if(m[4]!=0){aInd = c(aInd,(sum(m[1:3])+1):p)}
    coeff$alpha[aInd] = 0.8*tau/3#runif(length(aInd),0.5,2.5) #alpha_j ~ 0.2Unif(0.5,2.5)
    coeff$alpha=coeff$alpha[Index1]
    bInd = (sum(m[1:2])+1):p
    coeff$beta[bInd] = 1.2*tau/3#*runif(length(bInd),0.5,2.5) #beta_j ~ 0.3Unif(0.5,2.5)
    coeff$beta=coeff$beta[Index1]
    ####generate M: a matrix--n*p#####
    Em =mvrnorm(n,rep(0,times=p),matrix) #Em ~ N(0,1)
    # Em = matrix(rt(n*p, df = 3), n, p)*sqrt(1/3)
    # Em = matrix(rchisq(n*p, df=1)-1,n,p)*sqrt(2)
    # Em = matrix(rexp(n*p)-1,n,p) # Em ~ exp(1)
    M = matrix(X,n,1)%*%coeff$alpha+Em
    
    ####generate Y: a matrix--n*p#####
    betaX = 1
    #beta_j=0 for M_j
    EpY0 = rnorm(n) #Var(Y|x)=1
    Y = betaX*X + M%*%coeff$beta + EpY0
    
    return(list(X=X,Y=Y,M=M,Index1=Index1))
  }
  


  
  #-----------------Functions for DTM----------------------#
  

  ############################################################
  ############   Calculate the empirical FDR, power ##########
  ############################################################
  
  fdr_ap_fun <- function(S01,S02,W_rej,Index1){
    S01=S01[Index1]
    S02=S02[Index1]
    S_free <- 2*((S01<S02)+0)
    S_null <- (S01>S02)+0+S_free        ##here S_null=0 for alpha=beta;S_null=1 for alpha=0,beta!=0;S_null=2 for alpha!=0,beta=0 
    false_rej <- which((S_null==0)&(W_rej==1))
    true_rej <- which((S_null==1)&(W_rej==1))
    power <- length(true_rej)/max(sum(S_null==1),1)
    FDR <- length(false_rej)/max(sum(W_rej),1)
    fdr_ap_result <- list(FDR,power)
    return(fdr_ap_result)
  }
  
  ############################################################
  ####### Given the L2 threshold, it returns L1 threshold#####
  ############################################################
  
  W1_thre_new_fun <- function(W1,W2,L2,alpha,options,p){
    W1_abs_sort = sort(abs(W1))
    W1_abs_sort <- W1_abs_sort[which(W1_abs_sort!=0)]
    W_frame <- data.frame(ind=c(1:p),W1=W1,W2=W2)
    
    if(options=='+'){
      Ta = sapply(W1_abs_sort,function(x){
        (1+dim(W_frame[(W_frame$W1<=(-x))&(W_frame$W2<=L2),])[1] +  dim(W_frame[(W_frame$W1>=(x))&(W_frame$W2<=L2),])[1] -
           dim(W_frame[((W_frame$W2>=(-L2))&(W_frame$W2<0)&(W_frame$W1>=x)),])[1] - dim(W_frame[(W_frame$W2<=0)&(W_frame$W1>=x),])[1] )/
          max((dim(W_frame[((W_frame$W2<=(L2))&(W_frame$W1>=x)),])[1]),1)
      })
    }else{
      Ta = sapply(W1_abs_sort,function(x){
        (dim(W_frame[(W_frame$W1<=(-x))&(W_frame$W2<=L2),])[1] +  dim(W_frame[(W_frame$W1>=(x))&(W_frame$W2<=L2),])[1] -
           dim(W_frame[((W_frame$W2>=(-L2))&(W_frame$W2<0)&(W_frame$W1>=x)),])[1] - dim(W_frame[(W_frame$W2<=0)&(W_frame$W1>=x),])[1] )/
          max((dim(W_frame[((W_frame$W2<=(L2))&(W_frame$W1>=x)),])[1]),1)
      })
    }
    if(min(Ta)>alpha){
      W1_thre <- 1000000
      estimate_DTFDR=100
    }else{
      W1_thre = min(W1_abs_sort[which(Ta<=alpha)])
      ind=which(W1_abs_sort==min(W1_abs_sort[which(Ta<=alpha)]))
      estimate_DTFDR=Ta[[ind]]
    }
    
    
    det=W_frame[(W_frame$W1>=W1_thre) & (W_frame$W2<=L2) ,]$ind
    rej = rep(0,p)
    rej[det] = 1
    W_rej_list <- list(W1_thre,L2,rej) #return L1,L2,reject set
    return(list(W_rej_list=W_rej_list,estimate_DTFDR=estimate_DTFDR))
  }
  
  ############################################################
  # It inputs each value of W2 as the L2 threshold into the ##
  # above function and returns the best threshold pairs.   ###
  ############################################################
  
  opt_L1_L2_fun <- function(W1,W2,alpha,options,p){
    L2_can <- sort(abs(W2))
    L2_can <- L2_can[which(L2_can!=0)]
    num_rej <- sapply(L2_can,function(x){
      sum(W1_thre_new_fun(W1,W2,x,alpha,options,p)$W_rej_list[[3]]) 
    })
    
    L2 <- max(L2_can[which(num_rej==max(num_rej))])
    L1_L2_rej_list <- W1_thre_new_fun(W1,W2,L2,alpha,options,p)  
    return(L1_L2_rej_list)
  }
  
  ###########################################################
  ########     Function for linear regression      ##########
  ###########################################################
  
  regression=function(y,x){
    lm(y~x)
  }
  
  
  ###########################################################
  ##################generate the AR matrix###################
  ###########################################################
  
  AR_matrix<- function(n, rho) {
    exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                      (1:n - 1))
    rho^exponent
  }
  
  
  
  
  #----------Functions for detecting confounders------------------#
  
  
  DTM=function(X,M,Y,n,p,q,pi,B,Index1,ratio=0.5){
    result_rej2<-rep(0,p)
    All_result <- vector("list",length = B)
    for (b in 1:B)
    {
      #######DTM#############
      data=list("Y"=Y,'covariate'=cbind(X,M))
      S01=rep(0,time=p);S02=rep(0,time=p)
      support_S02=c((pi[1]*p+1):(sum(pi[1:2])*p),(sum(pi[1:3])*p+1):p)      ##True nonzero locations of beta
      S02[support_S02]=1
      support_S01=(sum(pi[1:2])*p+1):p       ##True nonzero locations of alpha
      S01[support_S01]=1
      
      ###Construct filter for beta
      D_1= sample(1:n,floor(n/2))
      D_2= (1:n)[-D_1]  
      part1 = glmnet::cv.glmnet(cbind(X[D_2],M[D_2,]),as.matrix(Y[D_2,]))
      ind0 = which(part1$nzero>1 & part1$nzero<(nrow(M[D_1,])*0.8))
      ind = which.min(part1$cvm[ind0]) #remove the case: no variable is significant
      beta = part1$glmnet.fit$beta[,ind0[ind]]
      
      lasso_coefficients1=beta
      support_1=which(lasso_coefficients1!=0)
      Y1=data$Y[D_2]
      x1=data$covariate[D_2,1]
      M1=data$covariate[D_2,setdiff(support_1,1)]
      reg_coefficients1=lm(Y1~x1+M1)$coefficients[-c(1,2)]
      W_11=rep(0,time=p+1)
      W_11[setdiff(support_1,1)]=reg_coefficients1
      W_11=W_11[-1]
      
      
      Y2=data$Y[D_1]
      x2=data$covariate[D_1,1]
      M2=data$covariate[D_1,setdiff(support_1,1)]
      reg_coefficients2=lm(Y2~x2+M2)$coefficients[-c(1,2)]
      W_12=rep(0,time=p+1)
      W_12[setdiff(support_1,1)]=reg_coefficients2
      W_12=W_12[-1]
      
      W_1=W_11*W_12
      
      ###Construct filter for alpha
      regression_list21=apply(data$covariate[D_1,-1],2,regression,x=data$covariate[D_1,1])                 #fit the data with  model(1)
      regression_list22=apply(data$covariate[D_2,-1],2,regression,x=data$covariate[D_2,1])
      coefficients21=lapply(regression_list21,function(x) x$coefficients[2])
      coefficients22=lapply(regression_list22,function(x) x$coefficients[2])
      W_2=unlist(coefficients21)*unlist(coefficients22)
      
      #choose proper thresholds L1,L2
      L1_L2_rej <- opt_L1_L2_fun(W_1,W_2,q,"",p)$W_rej_list
      L1 <- L1_L2_rej[[1]]
      L2 <- L1_L2_rej[[2]]
      W_rej <- L1_L2_rej[[3]]
      
      emp_FDR_power <- fdr_ap_fun(S01,S02,W_rej,Index1)
      emp_FDR <- emp_FDR_power[[1]]
      emp_power <- emp_FDR_power[[2]]
      
      result_DT<- list("W1_DT"=W_1,"W2_DT"=W_2,"L1_DT"=L1,"L2_DT"=L2,"reject_indexDT"=W_rej,
                       "FDR_DT"=emp_FDR,"Power_DT"=emp_power,"Support1_DT"=support_1,Index1=Index1
      )
      result_rej_b2<-result_DT[[5]]
      result_rej2<-cbind(result_rej2,result_rej_b2)
      All_result[[b]] <- result_DT
    }
    result_rej2<-result_rej2[,-1]
    sum_rej2<-rowSums(result_rej2)
    bag_rej2<-(sum_rej2>=ceiling(B/2))+0
    
    num<-sapply(1:B, function(x){
      sum((bag_rej2+result_rej2[,x]==2)+0) + sum((bag_rej2+result_rej2[,x]==0)+0)
    })
    final_decide<-which(num==max(num)) 
    final_decide<-final_decide[1]
    
    result_list <- All_result[[final_decide]] #list(emp_FDR,emp_power,result_rej)
    return(result_list) 
  }
  

  
  ############################################################
  ###############   start simulation  ########################
  ############################################################
  
  index=sample(p/2,p/2)
  AR_matrix=AR_matrix(p,correlation)
  Dat = ModelMultiReg(n,p,pi,AR_matrix,tau,index)     ###generate the data
  Index1=Dat$Index1
 
  Res=DTM(Dat$X,Dat$M,Dat$Y,n,p,q,pi,B,Index1,ratio=0.5)
  
  

  return(list(Dat = Dat,Res = Res,Index1=Index1))
  closeAllConnections()
}

##########################################
######## Simulation: Process #############
##########################################

set.seed(1000)

nonzero=(sum(pi[1:3])*p+1):p

print("running...")
start_time=proc.time()

RES0 = foreach(i=1:ns,.packages = c("MASS","lars")) %dopar% Process(n,p,pi,nonzero,q,tau,B)    
cat("模拟 time:", (proc.time()-start_time)[3][[1]],"秒","\n")

saveRDS(RES0, file=paste0("./Simulation/pro=_",pi01,"_corr=_",corr,"_n=_",n,"_p=_",p,"_t=_",tau,"_ratio=_",ratio,"_refine=_",B,"_mediation.RData"))

stopCluster(cl)

print("Down")
