rm(list=ls())
gc()
library(doParallel)
library(foreach)

# Settings
ns =100
no_cores <- detectCores()
cl <- makeCluster(no_cores-4)
clusterEvalQ(cl, {
  userlib <- path.expand("~/R/library")
  .libPaths(c(userlib, .libPaths()))
  NULL
})
registerDoParallel(cl)

n =400
p =2000
ratio=0.5
pi10 = 0.03
pi01 = 0.01
pi00 = 1-pi10-pi01*2
pi = c(pi00,pi10,pi01,1-sum(pi00,pi10,pi01))
B=30
nonzero = (sum(pi[1:3])*p+1):p
q = 0.1
iota = 2
corr=0.5

# Main simulation function
Process=function(n,p,pi,nonzero,q,iota,B,mode=0,correlation=corr,noise_sd=1){

  source("DACT_master.R")

  # Helper functions

  EstimateCoef <- function(X,M,Y,p=ncol(M)){

    EsAlpha= sapply(1:p,function(t){summary(lm(M[,t]~X))$coefficients[2,]})
    EsAlpha = data.frame(t(EsAlpha))

    if(ncol(Y)==p){
      EsBeta = sapply(1:p,function(t){summary(lm(Y[,t]~X+M[,t]))$coefficients[3,]})
    }else{
      EsBeta = sapply(1:p,function(t){summary(lm(Y~X+M[,t]))$coefficients[3,]})
    }

    EsBeta = data.frame(t(EsBeta))

    return(list(EsAlpha=EsAlpha,EsBeta=EsBeta))
  }

  BH <- function(pvalue, q){
    pvalue <- as.numeric(pvalue)
    p <- length(pvalue)
    o <- order(pvalue)
    sp <- pvalue[o]
    k <- max(which(sp <= (q * (1:p))/p), 0)
    if(k == 0) return(integer(0))
    o[1:k]
  }

  Estimate_zAB_candidates <- function(X, M, Y, cand){
    p <- ncol(M)
    zA <- rep(0, p)
    zB <- rep(0, p)

    for(j in cand){
      zA[j] <- summary(lm(M[,j] ~ X))$coefficients["X","t value"]
      zB[j] <- summary(lm(Y ~ X + M[,j]))$coefficients[2+1, "t value"]
    }
    list(zA=zA, zB=zB)
  }

  get_medScan_p <- function(res, p){

    if(!is.null(res$pvalues)) return(as.numeric(res$pvalues))
    if(!is.null(res$pvalue))  return(as.numeric(res$pvalue))
    if(!is.null(res$p))       return(as.numeric(res$p))
    if(is.data.frame(res) && any(grepl("p", names(res), ignore.case=TRUE))){
      nm <- names(res)[which(grepl("p", names(res), ignore.case=TRUE))[1]]
      return(as.numeric(res[[nm]]))
    }

    if(is.numeric(res) && length(res)==p) return(as.numeric(res))
    stop("Cannot find p-values in medScan result. Please print(str(res)).")
  }

  ols_beta_from_set <- function(X, M, Y, idx){
    p <- ncol(M)
    beta <- rep(0, p)
    idx <- sort(unique(idx))
    idx <- idx[idx>=1 & idx<=p]
    if(length(idx)==0) return(beta)

    Msel <- M[, idx, drop=FALSE]
    colnames(Msel) <- paste0("M", idx)

    df <- data.frame(Y=Y, X=X)
    df <- cbind(df, as.data.frame(Msel))

    fit <- lm(Y ~ ., data=df)
    coefs <- coef(fit)

    mnames <- paste0("M", idx)
    keep <- intersect(names(coefs), mnames)
    beta[as.integer(sub("^M","",keep))] <- coefs[keep]
    beta
  }

  DACT2020 <- function(pA, pB, q){

    pvalue <- DACT(pA,pB,correction= "NULL")
    det = BH(pvalue,q)

    return(det)
  }

  Estimate_zAB <- function(X, M, Y, cand){
    p <- ncol(M)
    zA <- rep(0, p)
    zB <- rep(0, p)

    for(j in cand){
      fit_a <- summary(lm(M[,j] ~ X))
      zA[j] <- fit_a$coefficients["X","t value"]

      fit_b <- summary(lm(Y ~ X + M[,j]))
      zB[j] <- fit_b$coefficients["M[, j]","t value"]
    }

    list(zA=zA, zB=zB)
  }

  # Data generation
  ModelMultiReg <- function(n,p,pi,matrix,iota,index){

    X = rbinom(n,1,0.2)

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
    coeff$alpha[aInd] = 0.8*3/3
    coeff$alpha=coeff$alpha[Index1]
    bInd = (sum(m[1:2])+1):p
    bInd1= (sum(m[1:2])+1):sum(m[1:3])
    bInd2= (sum(m[1:3])+1):p
    coeff$beta[bInd1] = 1.2*3*iota/3
    coeff$beta[bInd2] = 1.2
    coeff$beta=coeff$beta[Index1]

    Em =mvrnorm(n,rep(0,times=p),matrix)

    M = matrix(X,n,1)%*%coeff$alpha+Em

    betaX = 1

    EpY0 = rnorm(n)
    Y = betaX*X + M%*%coeff$beta + EpY0

    return(list(X=X,Y=Y,M=M,Index1=Index1))
  }

  AR_matrix<- function(n, rho) {
    exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) -
                      (1:n - 1))
    rho^exponent
  }

  fdr_ap_fun <- function(S01,S02,W_rej,Index1){
    S01=S01[Index1]
    S02=S02[Index1]
    S_free <- 2*((S01<S02)+0)
    S_null <- (S01>S02)+0+S_free
    false_rej <- which((S_null==0)&(W_rej==1))
    true_rej <- which((S_null==1)&(W_rej==1))
    power <- length(true_rej)/max(sum(S_null==1),1)
    FDR <- length(false_rej)/max(sum(W_rej),1)
    fdr_ap_result <- list(FDR,power)
    return(fdr_ap_result)
  }

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
    W_rej_list <- list(W1_thre,L2,rej)
    return(list(W_rej_list=W_rej_list,estimate_DTFDR=estimate_DTFDR))
  }

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

  regression=function(y,x){
    lm(y~x)
  }

  # DIMOC procedure
  DIMOC=function(X,M,Y,n,p,q,pi,Index1,ratio=0.5){
    result_rej2<-rep(0,p)
    All_result <- vector("list",length = B)
    for (b in 1:B)
    {

      data=list("Y"=Y,'covariate'=cbind(X,M))
      S01=rep(0,time=p);S02=rep(0,time=p)
      support_S02=c((pi[1]*p+1):(sum(pi[1:2])*p),(sum(pi[1:3])*p+1):p)
      S02[support_S02]=1
      support_S01=(sum(pi[1:2])*p+1):p
      S01[support_S01]=1

      D_1= sample(1:n,floor(n/2))
      D_2= (1:n)[-D_1]
      part1 = glmnet::cv.glmnet(cbind(X[D_2],M[D_2,]),as.matrix(Y[D_2,]))
      ind0 = which(part1$nzero>1 & part1$nzero<(nrow(M[D_1,])*0.8))
      ind = which.min(part1$cvm[ind0])
      beta = part1$glmnet.fit$beta[,ind0[ind]]

      lasso_coefficients1=beta
      support_1=which(lasso_coefficients1!=0)

      sel <- setdiff(support_1, 1)

      Y1 <- data$Y[D_2]
      x1 <- data$covariate[D_2, 1]

      if(length(sel) == 0){
        reg_coefficients1 <- numeric(0)
        W_11 <- rep(0, p)
      } else {
        M1 <- data$covariate[D_2, sel, drop=FALSE]
        fit1 <- lm(Y1 ~ x1 + M1)
        reg_coefficients1 <- coef(fit1)[-(1:2)]
        W_11 <- rep(0, p)
        W_11[sel - 1] <- reg_coefficients1
      }

      Y2 <- data$Y[D_1]
      x2 <- data$covariate[D_1, 1]

      if(length(sel) == 0){
        reg_coefficients2 <- numeric(0)
        W_12 <- rep(0, p)
      } else {
        M2 <- data$covariate[D_1, sel, drop=FALSE]
        fit2 <- lm(Y2 ~ x2 + M2)
        reg_coefficients2 <- coef(fit2)[-(1:2)]
        W_12 <- rep(0, p)
        W_12[sel - 1] <- reg_coefficients2
      }

      W_1=W_11*W_12

      regression_list21=apply(data$covariate[D_1,-1],2,regression,x=data$covariate[D_1,1])
      regression_list22=apply(data$covariate[D_2,-1],2,regression,x=data$covariate[D_2,1])
      coefficients21=lapply(regression_list21,function(x) x$coefficients[2])
      coefficients22=lapply(regression_list22,function(x) x$coefficients[2])
      W_2=unlist(coefficients21)*unlist(coefficients22)

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

    result_list <- All_result[[final_decide]]
    return(result_list)
  }

  EstimateCoef2 <- function(X,M,Y,p=ncol(M),nonbeta){

    EsAlpha= sapply(1:p,function(t){
      temp = summary(lm(M[,t]~X))
      temp$coefficients[,2] = temp$coefficients[,2]/temp$sigma
      return(temp$coefficients[2,])
    })

    EsAlpha = data.frame(t(EsAlpha))

    EsBeta = EsAlpha
    EsBeta$Estimate = 0;
    EsBeta$Std..Error = 1;
    EsBeta$t.value = 0;
    EsBeta$Pr...t.. = 1;

    lmYXM = summary(lm(Y~X+M[,nonbeta]))
    lmYXM$coefficients[,2] = lmYXM$coefficients[,2]/lmYXM$sigma
    EsBeta[nonbeta, ] = lmYXM$coefficients[-c(1,2),]

    return(list(EsAlpha=EsAlpha,EsBeta=EsBeta))
  }

  # Compared methods
  DACT_DIMOC<-function(X,M,Y,nonbeta,nonzero_med,nonzero_conf,Index1){
    p=dim(M)[2]
    n=dim(M)[1]
    EsAB = EstimateCoef2(X,M,Y, nonbeta=nonbeta)
    ind_mediator=DACT2020(EsAB$EsAlpha[,4],EsAB$EsBeta[,4], q)
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

  DACT_coe<-function(X,M,Y,Index1){
    p=dim(M)[2]
    n=dim(M)[1]
    ind_mediator_true=(sum(pi[1:3])*p+1):p
    ind_mediator_H0=setdiff((1:p),ind_mediator_true)
    ind_mediator_true=match(ind_mediator_true,Index1)
    ind_mediator_H0=match(ind_mediator_H0,Index1)
    EsAB = EstimateCoef2(X,M,Y, nonbeta=nonbeta)
    ind_mediator=DACT2020(EsAB$EsAlpha[,4],EsAB$EsBeta[,4], q)
    model_beta=lm(Y~X+M[,sort(ind_mediator) ] )
    coefficient_beta=rep(0,p)
    coefficient_beta[sort(ind_mediator )]=summary(model_beta)$coefficients[-c(1,2),1]
    TNR_DACT=length(setdiff(ind_mediator,ind_mediator_H0))/max(1,length(ind_mediator_H0))
    FDR_DACT=length(setdiff(ind_mediator,ind_mediator_true))/max(1,length(ind_mediator_true))
    Power_DACT=length(intersect(ind_mediator,ind_mediator_true))/length(ind_mediator_true)
    return(list(coefficient_beta=coefficient_beta,ind_mediator=ind_mediator,Power_DACT=Power_DACT,FDR_DACT=FDR_DACT,TNR_DACT=TNR_DACT) )
  }

  Sobel_coe <- function(X, M, Y, nonbeta, q=0.1){
    n <- nrow(M);
    p <- ncol(M)

    zA <- sapply(nonbeta, function(j) summary(lm(M[,j] ~ X))$coefficients["X","t value"] )
    zB <- sapply(nonbeta, function(j) summary(lm(Y ~ X + M[,j]))$coefficients[3,"t value"])

    sobel_res <- medScan::medScan(z.alpha=zA, z.beta=zB, method="Sobel")
    p_sub <- as.numeric(get_medScan_p(sobel_res, length(nonbeta)))

    sel_sub <- BH(p_sub, q=q)
    ind_mediator <- nonbeta[sel_sub]
    coefficient_beta <- ols_beta_from_set(X, M, Y, ind_mediator)
    list( coefficient_beta = coefficient_beta, ind_mediator = ind_mediator )
  }

  JST_coe <- function(X, M, Y, nonbeta, q=0.1){
    n <- nrow(M); p <- ncol(M)

    zA <- sapply(nonbeta, function(j){
      dfA <- data.frame(Mj = M[, j], X = X)
      summary(lm(Mj ~ X, data=dfA))$coefficients["X", "t value"] })
    zB <- sapply(nonbeta, function(j){ dfB <- data.frame(Y = Y, X = X, Mj = M[, j])

    summary(lm(Y ~ X + Mj, data=dfB))$coefficients["Mj", "t value"] })

    jst_res <- medScan::medScan(z.alpha = zA, z.beta = zB, method = "JT_comp")
    p_sub <- as.numeric(get_medScan_p(jst_res, length(nonbeta)))
    p_sub[is.na(p_sub)] <- 1

    sel_sub <- BH(p_sub, q = q)
    ind_mediator <- nonbeta[sel_sub]
    coefficient_beta <- ols_beta_from_set(X, M, Y, ind_mediator)
    list( coefficient_beta = coefficient_beta, ind_mediator = ind_mediator )
  }

  Sobel_DIMOC <- function(X, M, Y, nonbeta, Index1, q=0.1, pi, ratio=0.5){
    n <- nrow(M); p <- ncol(M)

    zA <- sapply(nonbeta, function(j) summary(lm(M[,j] ~ X))$coefficients["X","t value"] )
    zB <- sapply(nonbeta, function(j) summary(lm(Y ~ X + M[,j]))$coefficients[3,"t value"] )
    sobel_res <- medScan::medScan(z.alpha=zA, z.beta=zB, method="Sobel")
    p_sub <- as.numeric(get_medScan_p(sobel_res, length(nonbeta)))
    sel_sub <- BH(p_sub, q=q)
    ind_mediator <- nonbeta[sel_sub]

    res_DIMOC <- DIMOC(X, M, Y, n, p, q, pi, Index1, ratio=ratio)
    ind_confounder <- which(res_DIMOC[["reject_indexDT"]] != 0)
    ind_mediator <- setdiff(ind_mediator, ind_confounder)

    use_set <- sort(union(ind_mediator, ind_confounder))
    coefficient_beta <- ols_beta_from_set(X, M, Y, use_set)
    list( coefficient_beta = coefficient_beta, ind_mediator = ind_mediator, ind_confounder = ind_confounder )
  }

  JST_DIMOC<- function(X, M, Y, nonbeta, Index1, q=0.1, pi, ratio=0.5){
    if(!requireNamespace("medScan", quietly=TRUE)){ stop("Package 'medScan' not installed/loaded. Please install it first.") }
    n <- nrow(M); p <- ncol(M)

    zA <- sapply(nonbeta, function(j){ dfA <- data.frame(Mj = M[, j], X = X)
    summary(lm(Mj ~ X, data=dfA))$coefficients["X", "t value"] })
    zB <- sapply(nonbeta, function(j){ dfB <- data.frame(Y = Y, X = X, Mj = M[, j])
    summary(lm(Y ~ X + Mj, data=dfB))$coefficients["Mj", "t value"] })

    jst_res <- medScan::medScan(z.alpha = zA, z.beta = zB, method = "JT_comp")
    p_sub <- as.numeric(get_medScan_p(jst_res, length(nonbeta)))
    p_sub[is.na(p_sub)] <- 1

    sel_sub <- BH(p_sub, q = q)
    ind_mediator <- nonbeta[sel_sub]

    res_DIMOC <- DIMOC(X, M, Y, n, p, q, pi, Index1, ratio=ratio)
    ind_confounder <- which(res_DIMOC[["reject_indexDT"]] != 0)
    ind_mediator <- setdiff(ind_mediator, ind_confounder)

    use_set <- sort(union(ind_mediator, ind_confounder))
    coefficient_beta <- ols_beta_from_set(X, M, Y, use_set)
    list( coefficient_beta = coefficient_beta, ind_mediator = ind_mediator, ind_confounder = ind_confounder )
  }

  # One simulation run
  index=sample(p/2,p/2)
  AR_matrix=AR_matrix(p,correlation)
  Dat = ModelMultiReg(n,p,pi,AR_matrix,iota,index)
  Index1=Dat$Index1

  part1 = glmnet::cv.glmnet(cbind(Dat$X,Dat$M),Dat$Y)
  ind0 = which(part1$nzero>1 & part1$nzero<(nrow(Dat$M)*0.8))
  ind = which.min(part1$cvm[ind0])
  beta = part1$glmnet.fit$beta[,ind0[ind]]
  nonbeta = which(beta[-1]!=0)
  ind_mediator=(sum(pi[1:3])*p+1):p
  ind_confounder=(sum(pi[1:2])*p+1):(sum(pi[1:3])*p)

  coefficient_alpha=unlist(lapply(1:p,function(x){ summary(lm(Dat$M[,x]~Dat$X ))$coefficients[2,1] } ))
  coe_DACT_DIMOC=DACT_DIMOC(Dat$X,Dat$M,Dat$Y,nonbeta,nonzero_med=ind_mediator,nonzero_conf =ind_confounder ,Index1 = Index1)
  coe_DACT=DACT_coe(Dat$X,Dat$M,Dat$Y,Index1)
  coe_Sobel <- Sobel_coe(Dat$X, Dat$M, Dat$Y, nonbeta, q=q)
  coe_JST <- JST_coe(Dat$X, Dat$M, Dat$Y, nonbeta, q=q)
  coe_Sobel_DIMOC <- Sobel_DIMOC(Dat$X, Dat$M, Dat$Y, nonbeta=nonbeta, Index1=Index1, q=q, pi=pi, ratio=0.5)
  coe_JST_DIMOC <- JST_DIMOC(Dat$X, Dat$M, Dat$Y, nonbeta=nonbeta, Index1=Index1, q=q, pi=pi, ratio=0.5)

  return(list(coefficient_alpha=coefficient_alpha,coe_DACT_DIMOC=coe_DACT_DIMOC,coe_DACT=coe_DACT,coe_Sobel_DIMOC=coe_Sobel_DIMOC,coe_Sobel=coe_Sobel,coe_JST_DIMOC=coe_JST_DIMOC,coe_JST=coe_JST, Index1=Index1))
  }

# Simulation
set.seed(1000)
nonzero=(sum(pi[1:3])*p+1):p
print("running...")
start_time=proc.time()
RES0 <- foreach(i=1:ns, .packages=c("MASS","glmnet","medScan"), .errorhandling="pass") %dopar% { tryCatch( Process(n,p,pi,nonzero,q,iota,B), error=function(e){ list( error = TRUE, i = i, msg = conditionMessage(e), call = deparse(conditionCall(e)), calls = lapply(sys.calls(), deparse) ) } ) }

cat("Simulation time:", (proc.time()-start_time)[3][[1]], "seconds", "\n")
saveRDS(RES0,file=paste0("pro=_",pi01,"_corr=_",corr,"_n=_",n,"_p=_",p,"_t=_",iota,"_ratio=_",ratio,"_refine=_",B,"_mediation.RData"))

stopCluster(cl)
print("Down")
