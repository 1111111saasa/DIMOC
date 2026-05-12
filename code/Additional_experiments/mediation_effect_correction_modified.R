rm(list=ls())
gc()
library(doParallel)
library(foreach)
library(MASS)
library(hdi)
library(cp4p)

##########################################
######## Settings: Parallel ##############
#######################################
ns =100     ###simulation times
no_cores <- detectCores()
cl <- makeCluster(no_cores-2)
registerDoParallel(cl)

##########################################
######## Settings: Parameters ############
##########################################
n =400  ##sample size 
p =2000  ##dimension
ratio=0.5 ##splitting ratio
pi10 = 0.03  ##portion of alpha!=0,beta=0
pi01 = 0.01
pi00 = 1-pi10-pi01*2
pi = c(pi00,pi10,pi01,1-sum(pi00,pi10,pi01))
B=30
nonzero = (sum(pi[1:3])*p+1):p
q = 0.1 ##FDR level
iota = 3 ##tau=c(0.1,0.3,0.7)signal strength
corr=0.5  ##dependence



Process=function(n,p,pi,nonzero,q,iota,B,mode=0,correlation=corr,noise_sd=1){
  #-----------------Functions for DACT----------------------#
  
  
  source("C://Users//Administrator//Documents//DACT_master.R")
  
  EstimateCoef <- function(X,M,Y,p=ncol(M)){
    
    EsAlpha= sapply(1:p,function(t){summary(lm(M[,t]~X))$coefficients[2,]})
    EsAlpha = data.frame(t(EsAlpha)) #4 columns: 1-estimate coefficient; 2- std. error; 3-t value; 4-p value
    
    if(ncol(Y)==p){
      EsBeta = sapply(1:p,function(t){summary(lm(Y[,t]~X+M[,t]))$coefficients[3,]})    
    }else{
      EsBeta = sapply(1:p,function(t){summary(lm(Y~X+M[,t]))$coefficients[3,]})
    }
    
    EsBeta = data.frame(t(EsBeta))
    
    return(list(EsAlpha=EsAlpha,EsBeta=EsBeta))  
  }
  
  BH <- function(pvalue,p,q){
    orderp = order(pvalue)
    sortp = sort(pvalue)/(1:p)
    k = max(which(sortp<=(q/p)))
    det = NULL
    if(!is.infinite(k)){det = orderp[1:k]}
    
    return(det)
  }
  
  
  
  DACT2020 <- function(pA, pB, p, q){
    
    pvalue <- DACT(pA,pB,correction= "NULL")
    det = BH(pvalue,p,q)
    
    return(det)
  }
  
  
  
  #-----------------Functions for DIMOC----------------------#
  
  
  ####################################################################
  ############################generating data#########################
  ######Model: Mk=ak*X+ek;Y=r*X+sum_{k=1}^p bk*Mk+eps,k=1,...,p#######
  
  ModelMultiReg <- function(n,p,pi,matrix,iota,index){
    
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
    coeff$alpha[aInd] = 0.8*3/3#runif(length(aInd),0.5,2.5) #alpha_j ~ 0.2Unif(0.5,2.5)
    coeff$alpha=coeff$alpha[Index1]
    bInd = (sum(m[1:2])+1):p
    bInd1= (sum(m[1:2])+1):sum(m[1:3])
    bInd2= (sum(m[1:3])+1):p
    coeff$beta[bInd1] = 1.2*3*iota/3#*runif(length(bInd),0.5,2.5) #beta_j ~ 0.3Unif(0.5,2.5)
    coeff$beta[bInd2] = 1.2
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
  
  
  ###########################################################
  ##################generate the AR matrix###################
  ###########################################################
  
  AR_matrix<- function(n, rho) {
    exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                      (1:n - 1))
    rho^exponent
  }
  
  ############################################################
  ############   Calculate the empirical FDR, power ##########
  ############################################################
  
  fdr_ap_fun <- function(S01,S02,W_rej,Index1){
    S01=S01[Index1]
    S02=S02[Index1]
    S_free <- 2*((S01<S02)+0)
    S_null <- (S01>S02)+0+S_free        ##here S_null=0 for alpha=beta;S_null=1 for alpha!=0,beta=0;S_null=2 for alpha=0,beta!=0 
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
  
  
  DIMOC=function(X,M,Y,n,p,q,pi,Index1,ratio=0.5){
    result_rej2<-rep(0,p)
    All_result <- vector("list",length = B)
    for (b in 1:B)
    {
      #######DT#############
      data=list("Y"=Y,'covariate'=cbind(X,M))
      S01=rep(0,time=p);S02=rep(0,time=p)
      support_S02=c((pi[1]*p+1):(sum(pi[1:2])*p),(sum(pi[1:3])*p+1):p)      ##beta
      S02[support_S02]=1
      support_S01=(sum(pi[1:2])*p+1):p       ##alpha
      S01[support_S01]=1
      
      D_1= sample(1:n,floor(n/2))
      D_2= (1:n)[-D_1]  ##beta
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
      
      ###alpha
      regression_list21=apply(data$covariate[D_1,-1],2,regression,x=data$covariate[D_1,1])                 #fit the data with  model(1)
      regression_list22=apply(data$covariate[D_2,-1],2,regression,x=data$covariate[D_2,1])
      coefficients21=lapply(regression_list21,function(x) x$coefficients[2])
      coefficients22=lapply(regression_list22,function(x) x$coefficients[2])
      W_2=unlist(coefficients21)*unlist(coefficients22)
      
      #choose L1,L2
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
      result=result_DT
      All_result[[b]] <- result
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
  
  
  
  EstimateCoef2 <- function(X,M,Y,p=ncol(M),nonbeta){
    
    EsAlpha= sapply(1:p,function(t){
      temp = summary(lm(M[,t]~X))
      temp$coefficients[,2] = temp$coefficients[,2]/temp$sigma #no need to estimate sigma
      return(temp$coefficients[2,])
    })
    
    EsAlpha = data.frame(t(EsAlpha)) #4 columns: 1-estimate coefficient; 2- std. error; 3-t value; 4-p value
    
    ###Initialize EsBeta###
    EsBeta = EsAlpha
    EsBeta$Estimate = 0;
    EsBeta$Std..Error = 1;
    EsBeta$t.value = 0;
    EsBeta$Pr...t.. = 1; ###
    
    lmYXM = summary(lm(Y~X+M[,nonbeta]))
    lmYXM$coefficients[,2] = lmYXM$coefficients[,2]/lmYXM$sigma
    EsBeta[nonbeta, ] = lmYXM$coefficients[-c(1,2),]
    
    return(list(EsAlpha=EsAlpha,EsBeta=EsBeta))  
  }
  
  
  
  
  
  
  #----------------------HIMA2 class----------------------#
  HIMA2<-function(X,Y,M)
  {
    n <- dim(M)[1]  # number of samples
    p <- dim(M)[2]  # number of mediators
    d <- dim(X)[2]  # number of exposures
    q <- 0  # number of covariates
    Z=data.frame(matrix(ncol = 0, nrow = n))
    MZX<-cbind(M,Z,X)
    
    #########################################################################
    ########################### (Step 1) SIS step ########################### 
    #########################################################################
    message("Step 1: Sure Independent Screening ...", "  (", Sys.time(), ")")
    
    d_0 <- 2*round(n/log(n)) 
    beta_SIS <- matrix(0,1,p) 
    
    # Estimate the regression coefficients beta (mediators --> outcome)
    for (i in 1:p){
      ID_S <- c(i, (p+1):(p+q+1))
      MZX_SIS <- MZX[,ID_S]
      fit <- lsfit(MZX_SIS,Y,intercept = TRUE)
      beta_SIS[i] <- fit$coefficients[2]
    }
    
    # Estimate the regression coefficients alpha (exposure --> mediators)
    alpha_SIS <- matrix(0,1,p)
    XZ <- cbind(X,Z)
    for (i in 1:p){
      fit_a  <- lsfit(XZ,M[,i],intercept = TRUE)
      est_a <- matrix(coef(fit_a))[2]
      alpha_SIS[i] <- est_a
    }
    
    # Select the d_0 number of mediators with top largest effect 
    ab_SIS <- alpha_SIS*beta_SIS
    ID_SIS  <- which(-abs(ab_SIS) <= sort(-abs(ab_SIS))[d_0])
    d <- length(ID_SIS)
    
    #########################################################################
    ################### (Step 2) De-biased Lasso Estimates ##################
    #########################################################################
    message("Step 2: De-biased Lasso Estimates ...", "   (", Sys.time(), ")")
    
    P_beta_SIS <- matrix(0,1,d)
    beta_DLASSO_SIS_est <- matrix(0,1,d)
    beta_DLASSO_SIS_SE <- matrix(0,1,d)
    MZX_SIS <- MZX[,c(ID_SIS, (p+1):(p+q+1))]
    
    DLASSO_fit <- lasso.proj(x=MZX_SIS, y=Y, family = "gaussian",Z = NULL)
    beta_DLASSO_SIS_est <- DLASSO_fit$bhat[1:d]
    beta_DLASSO_SIS_SE <- DLASSO_fit$se
    P_beta_SIS <- t(DLASSO_fit$pval[1:d])
    
    ################### Estimate alpha ################
    alpha_SIS_est <- matrix(0,1,d)
    alpha_SIS_SE <- matrix(0,1,d)
    P_alpha_SIS <- matrix(0,1,d)
    
    XZ <- cbind(X,Z)
    for (i in 1:d){
      fit_a  <- lsfit(XZ,M[,ID_SIS[i]],intercept = TRUE)
      est_a <- matrix(coef(fit_a))[2]
      se_a <- ls.diag(fit_a)$std.err[2]
      sd_1 <- abs(est_a)/se_a
      P_alpha_SIS[i] <- 2*(1-pnorm(sd_1,0,1))  ## the SIS for alpha
      alpha_SIS_est[i] <- est_a
      alpha_SIS_SE[i] <- se_a
    }
    
    #########################################################################
    ################ (step 3) The multiple-testing  procedure ###############
    #########################################################################
    message("Step 3: Joint significance test ...", "     (", Sys.time(), ")")
    
    PA <- cbind(t(P_alpha_SIS),(t(P_beta_SIS)))
    P_value <- apply(PA,1,max)  #The joint p-values for SIS variable
    
    N0 <- dim(PA)[1]*dim(PA)[2]
    input_pvalues <- PA + matrix(runif(N0,0,10^{-10}),dim(PA)[1],2)
    input_pvalues=input_pvalues <- as.data.frame(apply(input_pvalues,2, unlist))
    
    # Estimate the proportions of the three component nulls
    nullprop <- null_estimation(input_pvalues)
    
    # Compute the estimated pointwise FDR for every observed p-max
    fdrcut  <- fdr_est(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10, nullprop$alpha1,nullprop$alpha2,input_pvalues,method=0)
    
    ID_fdr <- which(fdrcut <= 0.05)
    
    # Following codes extract the estimates for mediators with fdrcut<=0.05
    beta_hat_est <- beta_DLASSO_SIS_est[ID_fdr]
    beta_hat_SE  <- beta_DLASSO_SIS_SE[ID_fdr]
    
    alpha_hat_est <-  alpha_SIS_est[ID_fdr]
    alpha_hat_SE  <-  alpha_SIS_SE[ID_fdr]
    
    P.value_raw <- P_value[ID_fdr]
    
    # Indirect effect
    IDE <- beta_hat_est*alpha_hat_est # mediation(indirect) effect
    
    # Here we name the mediators as M1-Mp and extract the names of significant ones.
    M<-(sprintf("M%d", 1:p))[ID_SIS[ID_fdr]]
    
    # create a data frame with output values
    output<-data.frame(cbind(M, alpha=alpha_hat_est,alpha_SE=alpha_hat_SE,beta=beta_hat_est,beta_SE=beta_hat_SE,"alpha*beta"=IDE, 
                             p_val=P.value_raw,ID_med=ID_SIS[ID_fdr] ))
    
    message("Done!", "     (", Sys.time(), ")")
    
    return(output)
  }
  
  null_estimation <-function(input_pvalues,lambda=0.5) {
    
    ## input_pvalues is a matrix with 2 columns of p-values, the first column is p-value for exposure-mediator association, the second column is p-value for mediator-outcome association adjusted for exposure
    ## lambda is the threshold for pi_{00} estimation, default 0.5
    #check input
    #if (!class(input_pvalues) %in% c("matrix","data.frame"))
    #stop("input_pvalues should be a matrix or data frame")
    if (ncol(input_pvalues) !=2)
      stop("inpute_pvalues should have 2 column")
    #input_pvalues <- matrix(input_pvalues,nrow=nrow(input_pvalues)*2)
    if (sum(complete.cases(input_pvalues))<nrow(input_pvalues))
      warning("input_pvalues contains NAs to be removed from analysis")
    input_pvalues <- input_pvalues[complete.cases(input_pvalues),]
    
    
    pcut <- seq(0.1,0.8,0.1) 
    frac1 <- rep(0,8)
    frac2 <- rep(0,8)
    frac12<- rep(0,8)
    for (i in 1:8) {
      frac1[i] <- mean(input_pvalues[,1]>=pcut[i])/(1-pcut[i])
      frac2[i] <- mean(input_pvalues[,2]>=pcut[i])/(1-pcut[i]) 
      frac12[i]<- mean(input_pvalues[,2]>=pcut[i] & input_pvalues[,1]>=pcut[i])/(1-pcut[i])^2
    }  
    
    ## use the median estimates for pi00 ##
    
    alpha00 <- min(frac12[pcut==lambda],1)
    
    ## alpha1 is the proportion of nulls for first p-value 
    ## alpha2 is the proportion of nulls for second p-value 
    #library(cp4p)
    if (ks.test(input_pvalues[,1],"punif",0,1)$p>0.05) alpha1 <- 1 else     alpha1 <- estim.pi0(p=input_pvalues[,1], pi0.method="slim")$pi0.SLIM
    if (ks.test(input_pvalues[,2],"punif",0,1)$p>0.05) alpha2 <- 1 else     alpha2 <- estim.pi0(p=input_pvalues[,2], pi0.method="slim")$pi0.SLIM
    
    if (alpha1 == 1 | alpha2==1) {
      alpha11 <- 0
      alpha10 <- 1-alpha1
      alpha01 <- 1-alpha2
      if (alpha1==1 & alpha2==1) alpha00 <- 1 else {
        if (sum(alpha00+alpha10+alpha01)>1) {
          err <- sum(alpha00+alpha10+alpha01)-1
          alpha00 <- alpha00- err 
        }  else {
          err <- 1- sum(alpha00+alpha10+alpha01)
          if (alpha10 !=0) alpha10 <- alpha10 + err
          if (alpha01 !=0) alpha01 <- alpha01 + err
        }      
      }
    } else {
      
      alpha11 <- max(1-alpha1+1-alpha2+alpha00-1,0)
      if (alpha11 ==0) {
        alpha10 <- 1-alpha1 
        alpha01 <- 1-alpha2 
        err <- 1-(alpha00+alpha10+alpha01)
        alpha00 <- alpha00 + err
        
      }  else {
        if (alpha11> (1-alpha1)| alpha11> (1-alpha2)) {
          alpha11 <- 0
          alpha10 <- 1-alpha1
          alpha01 <- 1-alpha2
          
          err <- 1-(alpha00+alpha10+alpha01)
          alpha00 <- alpha00 + err
          
        } else {
          alpha10 <- max(1-alpha1-alpha11,0)
          alpha01 <- max(1-alpha2-alpha11,0)
          alpha00 <- 1-alpha10-alpha01-alpha11
        }
      }
    }
    alpha.null <- list(alpha10=alpha10,alpha01=alpha01,alpha00=alpha00,alpha1=alpha1,alpha2=alpha2)
    return(alpha.null)
  }
  
  fdr_est <-function(alpha00,alpha01,alpha10,alpha1,alpha2,input_pvalues, method=0){
    
    ## alpha10,alpha01,alpha00 are null proportions
    ## alpha1 is the marginal null proportion for first p-value
    ## alpha2 is the marginal null proportion for second p-value
    ## input pvalues are two columns of p-values
    ## alpha is the level of FWER to be control at   
    ## method=0 corresponding to the approximation used in section 2.2-2.3 in the paper, the default value is 0  
    ## method=1 corresponding to the exact used in section 2.4 in the paper
    #check input
    
    #input_pvalues <- input_pvalues[complete.cases(input_pvalues),]
    # if (!is.null(nrow(input_pvalues)) & nrow(input_pvalues)<1)
    # stop("input_pvalues doesn't have valid p-values")
    
    pmax <- apply(input_pvalues,1,max)
    nmed <- length(pmax)
    efdr1 <- rep(0,nmed)
    
    for (i in 1:nmed) {
      fdr11 <-  (pmax[i]*alpha01)/mean(pmax<=pmax[i])
      fdr12 <-  (pmax[i]*alpha10)/mean(pmax<=pmax[i])          
      fdr2  <-  (pmax[i]*pmax[i]*alpha00)/mean(pmax<=pmax[i])   
      efdr1[i] <- fdr11+fdr12+fdr2
    }  
    
    if (method==1) {
      #library(fdrtool)
      nmed  <- nrow(input_pvalues)  
      cdf12 <- input_pvalues
      
      xx1 <- c(0,input_pvalues[order(input_pvalues[,1]),1])
      yy1 <- c(0,seq(1,nmed,by=1)/nmed)
      gfit1<- gcmlcm(xx1,yy1,type="lcm")
      xknots1 <- gfit1$x.knots[-1]
      Fknots1 <- cumsum(diff(gfit1$x.knots)*gfit1$slope.knots)
      
      xx2 <- c(0,input_pvalues[order(input_pvalues[,2]),2])
      yy2 <- c(0,seq(1,nmed,by=1)/nmed)
      gfit2<- gcmlcm(xx2,yy2,type="lcm")
      xknots2 <- gfit2$x.knots[-1]
      Fknots2 <- cumsum(diff(gfit2$x.knots)*gfit2$slope.knots)
      
      if (alpha1!=1) Fknots1 <- (Fknots1 - alpha1*xknots1)/(1-alpha1) else Fknots1 <- rep(0,length(xknots1))
      if (alpha2!=1) Fknots2 <- (Fknots2 - alpha2*xknots2)/(1-alpha2) else Fknots2 <- rep(0,length(xknots2))
      
      
      orderq1 <- pmax
      orderq2 <- pmax
      
      gcdf1 <- pmax
      gcdf2 <- pmax
      for (i in 1:length(xknots1)) {
        if (i==1) {
          gcdf1[orderq1<=xknots1[i]] <- (Fknots1[i]/xknots1[i])*orderq1[orderq1<=xknots1[i]] 
        } else {   
          if (sum(orderq1>xknots1[i-1] & orderq1<=xknots1[i])>0){
            temp <- orderq1[orderq1>xknots1[i-1] & orderq1<=xknots1[i]] 
            gcdf1[orderq1>xknots1[i-1] & orderq1<=xknots1[i]] <- Fknots1[i-1] + (Fknots1[i]-Fknots1[i-1])/(xknots1[i]-xknots1[i-1])*(temp-xknots1[i-1])
          }
        }
      }
      
      for (i in 1:length(xknots2)) {
        if (i==1) {
          gcdf2[orderq2<=xknots2[i]] <- (Fknots2[i]/xknots2[i])*orderq2[orderq2<=xknots2[i]] 
        } else {   
          if (sum(orderq2>xknots2[i-1] & orderq2<=xknots2[i])>0){
            temp <- orderq2[orderq2>xknots2[i-1] & orderq2<=xknots2[i]] 
            gcdf2[orderq2>xknots2[i-1] & orderq2<=xknots2[i]] <- Fknots2[i-1] + (Fknots2[i]-Fknots2[i-1])/(xknots2[i]-xknots2[i-1])*(temp-xknots2[i-1])
          } 
        }
      }
      
      
      gcdf1 <- ifelse(gcdf1>1,1,gcdf1)
      gcdf2 <- ifelse(gcdf2>1,1,gcdf2)
      
      cdf12[,1] <- gcdf1
      cdf12[,2] <- gcdf2
      
      for (i in 1:nmed) {
        fdr11 <-  (pmax[i]*cdf12[i,2]*alpha01)/mean(pmax<=pmax[i])
        fdr12 <-  (pmax[i]*cdf12[i,1]*alpha10)/mean(pmax<=pmax[i])          
        fdr2  <-  (pmax[i]*pmax[i]*alpha00)/mean(pmax<=pmax[i])   
        efdr1[i] <- fdr11+fdr12+fdr2
      }  
    }
    
    efdr1.order <- efdr1[order(pmax,decreasing=T)]
    
    for (i in 2:nmed)  {
      efdr1.order[i] <- min(efdr1.order[i],efdr1.order[i-1])
    }
    
    efdr1 <- efdr1.order[rank(-pmax)]
    
    return(efdr=efdr1)
  }
  
  
  
  
  
  
  #-----------------compared methods----------------------#
  
  
  DACT_DIMOC<-function(X,M,Y,nonbeta,nonzero_med,nonzero_conf,Index1){
    p=dim(M)[2]
    n=dim(M)[1]
    EsAB = EstimateCoef2(X,M,Y, nonbeta=nonbeta)
    ind_mediator=DACT2020(EsAB$EsAlpha[,4],EsAB$EsBeta[,4], p, q)
    res_DIMOC=DIMOC(X,M,Y,n,p,q,pi,Index1,ratio=0.5)
    ind_confounder=res_DIMOC[["reject_indexDT"]]
    ind_confounder=which(ind_confounder!=0)
    ind_mediator=setdiff(ind_mediator,ind_confounder)
    
    
    coefficient_beta=rep(0,p)
    model_beta=lm(Y~X+M[,sort(union(ind_mediator,ind_confounder)) ] )
    coefficient_beta[sort(union(ind_mediator,ind_confounder) )]=summary(model_beta)$coefficients[-c(1,2),1]
    
    
    nonzero_conf=match(nonzero_conf,Index1)
    nonzero_med =match(nonzero_med,Index1)
    ind_mediator_H0=setdiff((1:p),nonzero_med )
    fdr_conf=length(setdiff(ind_confounder,nonzero_conf))/max(1,length(ind_confounder))
    fdr_med=length(setdiff(ind_mediator,nonzero_med))/max(1,length(ind_mediator))
    power_conf=length(intersect(ind_confounder,nonzero_conf))/length(nonzero_conf)
    power_med=length(intersect(ind_mediator,nonzero_med))/length(nonzero_med)
    TNR_med=length(intersect(ind_mediator_H0,ind_mediator))/max(1,length(ind_mediator_H0))
    return(list(TNR_med=TNR_med,coefficient_beta=coefficient_beta,ind_confounder=ind_confounder,ind_mediator=ind_mediator, fdr_conf=fdr_conf,fdr_med=fdr_med,power_conf=power_conf,power_med=power_med) ) 
  }
  
  truemed_DIMOC<-function(X,M,Y,ind_confounder,Index1){
    p=dim(M)[2]
    ind_mediator=(sum(pi[1:3])*p+1):p
    ind_mediator=match(ind_mediator,Index1)
    ind_confounder_true=(sum(pi[1:2])*p+1):(sum(pi[1:3])*p)
    ind_confounder_true=match(ind_confounder_true,Index1)
    ind_confounder_H0=setdiff((1:p),ind_confounder_true )
    ind_confounder_H0=match(ind_confounder_H0,Index1)
    coefficient_beta=rep(0,p)
    model_beta=lm(Y~X+M[,sort(union(ind_mediator,ind_confounder))  ]  )
    coefficient_beta[sort(union(ind_mediator,ind_confounder))]=summary(model_beta)$coefficients[-c(1,2),1]
    ind_confounder=setdiff(ind_confounder,ind_mediator)
    TNR_DIMOC=length(setdiff(ind_confounder,ind_confounder_H0))/max(1,length(ind_confounder_H0))
    FDR_DIMOC=length(setdiff(ind_confounder,ind_confounder_true))/max(1,length(ind_confounder))
    Power_DIMOC=length(intersect(ind_confounder,ind_confounder_true))/length(ind_confounder_true)
    
    
    
    return(list(coefficient_beta=coefficient_beta,TNR_DIMOC=TNR_DIMOC,FDR_DIMOC=FDR_DIMOC,Power_DIMOC=Power_DIMOC,ind_confounder=ind_confounder) )
  }
  
  oracle<-function(X,M,Y,ind_mediator,ind_confounder,Index1){
    n=dim(M)[1]
    p=dim(M)[2]
    ind_mediator=match(ind_mediator,Index1)
    ind_confounder=match(ind_confounder,Index1)
    model_beta=lm(Y~X+M[,c(ind_mediator,ind_confounder)])
    coefficient_beta=rep(0,p)
    coefficient_beta[c(ind_mediator,ind_confounder)]=summary(model_beta)$coefficients[-c(1,2),1]
    return(coefficient_beta)  
  }
  
  
  DACT_coe<-function(X,M,Y,Index1){
    p=dim(M)[2]
    n=dim(M)[1]
    ind_mediator_true=(sum(pi[1:3])*p+1):p
    ind_mediator_H0=setdiff((1:p),ind_mediator_true)
    ind_mediator_true=match(ind_mediator_true,Index1)
    ind_mediator_H0=match(ind_mediator_H0,Index1)
    EsAB = EstimateCoef2(X,M,Y, nonbeta=nonbeta)
    ind_mediator=DACT2020(EsAB$EsAlpha[,4],EsAB$EsBeta[,4], p, q)
    model_beta=lm(Y~X+M[,sort(ind_mediator) ] )
    coefficient_beta=rep(0,p)
    coefficient_beta[sort(ind_mediator )]=summary(model_beta)$coefficients[-c(1,2),1]
    TNR_DACT=length(setdiff(ind_mediator,ind_mediator_H0))/max(1,length(ind_mediator_H0))
    FDR_DACT=length(setdiff(ind_mediator,ind_mediator_true))/max(1,length(ind_mediator_true))
    Power_DACT=length(intersect(ind_mediator,ind_mediator_true))/length(ind_mediator_true)
    return(list(coefficient_beta=coefficient_beta,ind_mediator=ind_mediator,Power_DACT=Power_DACT,FDR_DACT=FDR_DACT,TNR_DACT=TNR_DACT) ) 
  }
  
  
  DACT_penalty_coe <- function(X, M, Y, Index1, nonbeta,
                               q = 0.1,
                               alpha = 1,          # lasso
                               nfolds = 10,
                               use_lambda = c("lambda.1se", "lambda.min"),
                               standardize = TRUE,
                               seed = NULL) {
    use_lambda <- match.arg(use_lambda)
    
    n <- nrow(M)
    p <- ncol(M)
    
    EsAB <- EstimateCoef2(X, M, Y, nonbeta = nonbeta)
    ind_mediator <- DACT2020(EsAB$EsAlpha[,4], EsAB$EsBeta[,4], p, q)
    
    if (length(ind_mediator) == 0) {
      xmat <- cbind(X, M)
      penalty.factor <- c(0, rep(1, p))  # X 不惩罚
    } else {
      xmat <- cbind(X, M)
      penalty.factor <- rep(1, 1 + p)
      penalty.factor[1] <- 0                          # X unpenalized
      penalty.factor[1 + ind_mediator] <- 0           # selected mediators unpenalized
    }
    
    if (!is.null(seed)) set.seed(seed)
    
    fit_cv <- glmnet::cv.glmnet(
      x = xmat,
      y = as.numeric(Y),
      family = "gaussian",
      alpha = alpha,
      penalty.factor = penalty.factor,
      standardize = standardize,
      nfolds = nfolds
    )
    
    s <- if (use_lambda == "lambda.1se") fit_cv$lambda.1se else fit_cv$lambda.min
    beta_all <- as.vector(glmnet::coef.glmnet(fit_cv$glmnet.fit, s = s)) # (Intercept, X, M1..Mp)
  
    coefficient_beta <- rep(0, p)
    coefficient_beta <- beta_all[-c(1,2)] 
    ind_nonzero_final <- which(abs(coefficient_beta) > 0)
  
    return(list(
      coefficient_beta   = coefficient_beta,
      ind_mediator       = ind_mediator,
      ind_nonzero_final  = ind_nonzero_final,
      penalty_factor_M   = penalty.factor[-1], 
      lambda_used        = s,
      cvm_min            = min(fit_cv$cvm),
      fit_cv             = fit_cv
    ))
  }
  
  
  
  eye_true<-function(X,M,Y,nonbeta,Index1=Index1){
    p=dim(M)[2]
    n=dim(M)[1]
    ind_mediator=(sum(pi[1:3])*p+1):p
    ind_mediator=match(ind_mediator,Index1)
    
    regression=function(y,x){
      summary(lm(y~x))
    }
    support_1=nonbeta
    
    regression_list21=apply(M,2,regression,x=X)                 #fit the data with  model(1)
    coefficients21=lapply(regression_list21,function(x) ifelse(x$coefficients[2,4]<=q,x$coefficients[2,1],0))
    nonsupport_2=which(unlist(coefficients21)==0)
    ind_confounder <- intersect(support_1,nonsupport_2)
    coefficient_beta=rep(0,p)
    model_beta=lm(Y~X+M[,sort(union(ind_mediator,ind_confounder))  ]  )
    coefficient_beta[c(ind_mediator,ind_confounder)]=summary(model_beta)$coefficients[-c(1,2),1]
    return(list(coefficient_beta=coefficient_beta,ind_confounder=ind_confounder))
    
  }
  
  
  #-----------------simulation---------------------#
  index=sample(p/2,p/2)
  AR_matrix=AR_matrix(p,correlation)
  Dat = ModelMultiReg(n,p,pi,AR_matrix,iota,index)     ###generate the data
  Index1=Dat$Index1
  
  
  part1 = glmnet::cv.glmnet(cbind(Dat$X,Dat$M),Dat$Y)
  ind0 = which(part1$nzero>1 & part1$nzero<(nrow(Dat$M)*0.8))
  ind = which.min(part1$cvm[ind0]) #remove the case: no variable is significant
  beta = part1$glmnet.fit$beta[,ind0[ind]]
  nonbeta = which(beta[-1]!=0)
  ind_mediator=(sum(pi[1:3])*p+1):p
  ind_confounder=(sum(pi[1:2])*p+1):(sum(pi[1:3])*p)
  
  
  coefficient_alpha=unlist(lapply(1:p,function(x){ summary(lm(Dat$M[,x]~Dat$X ))$coefficients[2,1]  }   ))
  coe_DACT_DIMOC=DACT_DIMOC(Dat$X,Dat$M,Dat$Y,nonbeta,nonzero_med=ind_mediator,nonzero_conf =ind_confounder ,Index1 = Index1)
  coe_truemed_DIMOC=truemed_DIMOC(Dat$X,Dat$M,Dat$Y,coe_DACT_DIMOC$ind_confounder,Index1=Index1)
  coe_oracle=oracle(Dat$X,Dat$M,Dat$Y,ind_mediator,ind_confounder,Index1=Index1)
  coe_HIMA2_res=HIMA2(Dat$X,Dat$Y,Dat$M)
  coe_DACT=DACT_coe(Dat$X,Dat$M,Dat$Y,Index1)
  coe_eye_true=eye_true(Dat$X,Dat$M,Dat$Y,nonbeta,Index1=Index1)
  coe_DACT_penalty <- DACT_penalty_coe(Dat$X, Dat$M, Dat$Y, Index1, nonbeta, q=q)
  return(list(Dat=Dat,coefficient_alpha=coefficient_alpha,coe_DACT_DIMOC=coe_DACT_DIMOC,coe_truemed_DIMOC=coe_truemed_DIMOC,coe_oracle=coe_oracle,coe_HIMA2=coe_HIMA2_res,coe_DACT=coe_DACT,coe_DACT_penalty=coe_DACT_penalty,coe_eye_true=coe_eye_true,Index1=Index1))
}




##########################################
######## Simulation: Process #############
##########################################

set.seed(1000)

nonzero=(sum(pi[1:3])*p+1):p
RES0 = foreach(i=1:ns,.packages = c("MASS","lars","cp4p","hdi")) %dopar% Process(n,p,pi,nonzero,q,iota,B)   
saveRDS(RES0,"n=400_p=2000.RData")




