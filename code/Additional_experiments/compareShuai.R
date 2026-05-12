rm(list=ls())
gc()
library(doParallel)
library(foreach)

# Parallel settings
ns =100
no_cores <- detectCores()
cl <- makeCluster(no_cores-4)
clusterEvalQ(cl, {
  userlib <- path.expand("~/R/library")
  .libPaths(c(userlib, .libPaths()))
  NULL
})
registerDoParallel(cl)

# Simulation settings
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
iota = 1
corr=0.7

# Main process
Process=function(n,p,pi,nonzero,q,iota,B,mode=0,correlation=corr,noise_sd=1){

  # DACT functions
  source("DACT_master.R")

  DACT2020 <- function(pA, pB, q){

    pvalue <- DACT(pA,pB,correction= "NULL")
    det = BH(pvalue,q)

    return(det)
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

  # Utility functions
  AR_matrix<- function(n, rho) {
    exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) -
                      (1:n - 1))
    rho^exponent
  }

  # DIMOC functions
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
  DACT_DIMOC_drop <- function(X, M, Y,
                               nonbeta,
                               nonzero_med, nonzero_conf,
                               Index1,
                               drop_idx,
                               q = 0.1,
                               pi, ratio=0.5){

    p <- ncol(M); n <- nrow(M)

    drop_idx <- sort(unique(drop_idx))

    conf_true <- match(nonzero_conf, Index1)

    keep <- setdiff(1:p, drop_idx)
    M_obs <- M[, keep, drop=FALSE]

    nonbeta_obs <- match(nonbeta, keep)
    nonbeta_obs <- nonbeta_obs[!is.na(nonbeta_obs)]
    if(length(nonbeta_obs) == 0) nonbeta_obs <- 1

    EsAB <- EstimateCoef2(X, M_obs, Y, nonbeta = nonbeta_obs)
    ind_mediator <- DACT2020(EsAB$EsAlpha[,4], EsAB$EsBeta[,4], q)

    res_DIMOC <- DIMOC(X, M_obs, Y, n=n, p=ncol(M_obs), q=q, pi=pi, Index1=1:ncol(M_obs), ratio=ratio)
    ind_confounder <- which(res_DIMOC[["reject_indexDT"]] != 0)

    ind_mediator <- setdiff(ind_mediator, ind_confounder)

    coefficient_beta_obs <- rep(0, ncol(M_obs))
    use_set <- sort(union(ind_mediator, ind_confounder))
    if(length(use_set) > 0){
      fit <- lm(Y ~ X + M_obs[, use_set, drop=FALSE])
      coefficient_beta_obs[use_set] <- coef(fit)[-(1:2)]
    }

    coefficient_beta <- rep(0, p)
    coefficient_beta[keep] <- coefficient_beta_obs
    ind_mediator_full <- keep[ind_mediator]
    ind_confounder_full <- keep[ind_confounder]

    list(
      drop_idx = drop_idx,
      keep_idx = keep,
      coefficient_beta = coefficient_beta,
      ind_mediator = ind_mediator_full,
      ind_confounder = ind_confounder_full
    )
  }

  Shuai <- function(X, M, Y,
                    Index1,
                    drop_idx = integer(0),
                    r = 1,
                    nfolds = 10,
                    alpha_glmnet = 1,
                    tau = 2,
                    fm = "minres",
                    alpha_p = 0.05) {

    suppressPackageStartupMessages({
      library(glmnet)
      library(psych)
    })

    n <- nrow(M)
    p <- ncol(M)
    X <- as.numeric(X)

    drop_idx <- sort(unique(drop_idx))
    drop_idx <- drop_idx[drop_idx >= 1 & drop_idx <= p]
    keep <- setdiff(1:p, drop_idx)

    M_obs <- as.matrix(M[, keep, drop = FALSE])
    p_obs <- ncol(M_obs)
    if(p_obs < 2) stop("Too few observed mediators after dropping.")

    sds <- apply(M_obs, 2, sd, na.rm = TRUE)
    keep2 <- which(is.finite(sds) & sds > 1e-8)
    if(length(keep2) < 2) stop("Too few non-constant mediators after filtering.")
    M2 <- M_obs[, keep2, drop = FALSE]
    keep_obs_cols <- keep[keep2]
    p2 <- ncol(M2)

    alpha_pv <- rep(1, p2)
    for(j in 1:p2){
      s <- summary(lm(M2[, j] ~ X))
      alpha_pv[j] <- s$coefficients["X", "Pr(>|t|)"]
    }
    ind_alpha <- which(alpha_pv < alpha_p)

    M_fa <- matrix(0, n, p2)
    for(j in 1:p2){
      M_fa[, j] <- residuals(lm(M2[, j] ~ X))
    }

    EU_hat <- rep(0, n)
    if(r > 0){
      fa_fit <- psych::fa(scale(M_fa), nfactors = r, fm = fm,
                          rotate = "none", scores = "regression")
      U_hat <- as.matrix(fa_fit$scores)
      if(!is.null(U_hat) && all(is.finite(U_hat)) && ncol(U_hat) >= 1){
        EU_hat <- U_hat[, 1]
      }
    }

    data1 <- cbind(X, M2, EU_hat)

    pen1 <- c(0, rep(1, p2), 0)
    fit1 <- cv.glmnet(
      x = data1 / sqrt(n),
      y = as.numeric(Y) / sqrt(n),
      nfolds = nfolds,
      alpha = alpha_glmnet,
      penalty.factor = pen1,
      standardize = TRUE,
      intercept = TRUE
    )
    b1 <- as.numeric(coef(fit1, s = "lambda.min"))
    beta_init <- b1[3:(p2 + 2)]

    pen2 <- c(0, 1 / (abs(beta_init) + 1/n)^tau, 0)
    fit2 <- cv.glmnet(
      x = data1 / sqrt(n),
      y = as.numeric(Y) / sqrt(n),
      nfolds = nfolds,
      alpha = alpha_glmnet,
      penalty.factor = pen2,
      standardize = TRUE,
      intercept = TRUE
    )
    b2 <- as.numeric(coef(fit2, s = "lambda.min"))
    beta_hat <- b2[3:(p2 + 2)]
    ind_beta <- which(beta_hat != 0)

    ind_mediator_local <- intersect(ind_alpha, ind_beta)

    coefficient_beta <- rep(0, p)
    coefficient_beta[keep_obs_cols] <- beta_hat

    ind_mediator <- keep_obs_cols[ind_mediator_local]

    list(
      coefficient_beta = coefficient_beta,
      ind_mediator = ind_mediator,
      Index1 = Index1
    )
  }

  # Single simulation
  index=sample(p/2,p/2)
  AR_matrix=AR_matrix(p,correlation)
  Dat = ModelMultiReg(n,p,pi,AR_matrix,iota,index)
  Index1=Dat$Index1

  confounder_raw <- (sum(pi[1:2]) * p + 1):(sum(pi[1:3]) * p)
  confounder_index1 <- match(confounder_raw, Index1)
  confounder_index1 <- confounder_index1[!is.na(confounder_index1)]
  drop_idx <- sample(confounder_index1, 5, replace = FALSE)

  part1 = glmnet::cv.glmnet(cbind(Dat$X,Dat$M),Dat$Y)
  ind0 = which(part1$nzero>1 & part1$nzero<(nrow(Dat$M)*0.8))
  ind = which.min(part1$cvm[ind0])
  beta = part1$glmnet.fit$beta[,ind0[ind]]
  nonbeta = which(beta[-1]!=0)
  ind_mediator=(sum(pi[1:3])*p+1):p
  ind_confounder=(sum(pi[1:2])*p+1):(sum(pi[1:3])*p)

  coefficient_alpha=unlist(lapply(1:p,function(x){ summary(lm(Dat$M[,x]~Dat$X ))$coefficients[2,1]  }   ))
  coe_DACT_DIMOC <- DACT_DIMOC_drop(
    Dat$X, Dat$M, Dat$Y,
    nonbeta = nonbeta,
    nonzero_med = ind_mediator,
    nonzero_conf = ind_confounder,
    Index1 = Index1,
    drop_idx = drop_idx,
    q = q, pi = pi, ratio = ratio
  )

  coe_Shuai <- Shuai(
    X = Dat$X,
    M = Dat$M,
    Y = Dat$Y,
    Index1 = Index1,
    drop_idx = drop_idx,
    r = 5,
    nfolds = 10
  )

  return(list(coefficient_alpha=coefficient_alpha,coe_DACT_DIMOC=coe_DACT_DIMOC, coe_Shuai=coe_Shuai,Index1=Index1))
}

# Run simulations
set.seed(1000)

nonzero=(sum(pi[1:3])*p+1):p

print("running...")
start_time=proc.time()
RES0 <- foreach(i=1:ns,
                .packages=c("MASS","glmnet","medScan","psych"),
                .errorhandling="pass") %dopar% {
                  tryCatch(
                    Process(n,p,pi,nonzero,q,iota,B),
                    error=function(e){
                      list(
                        error = TRUE,
                        i = i,
                        msg = conditionMessage(e),
                        call = deparse(conditionCall(e)),
                        calls = lapply(sys.calls(), deparse)
                      )
                    }
                  )
                }

cat("Simulation time:", (proc.time()-start_time)[3][[1]], "seconds", "\n")

saveRDS(RES0,file=paste0("compareShuai_pro=_",pi01,"_corr=_",corr,"_n=_",n,"_p=_",p,"_t=_",iota,"_ratio=_",ratio,"_refine=_",B,"_mediation.RData"))

stopCluster(cl)

print("Down")
