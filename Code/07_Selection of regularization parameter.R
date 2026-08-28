################################################################################
# Selection of regularization parameter lambda_xi
# For estimating the derivative-based scores
################################################################################

# lambda_xi setting: lambda_xi=0, 0.5*hat{lambda_xi}, hat{lambda_xi}, 2*hat{lambda_xi} 
# Aim1: regularization improves the estimation of derivative-based scores
# Aim2: estimation of derivative-based scores remains stable for different level of lambda_xi

library(fda)
library(funData)
library(fdapace)
library(refund)
library(foreach)
library(doParallel)
source("01_dgp.R")
source("02_method.R")
source("03_errors.R")

# new1 ----
# lambda_xi=0 
deriv_fpca_fun_new1 <- function(yraw, argvals, d1, d2, K.p){
  # yraw: Lxn observations
  # argvals: a vector containing observed locations on the functional domain
  # d1/d2: difference penalty order
  # K.p: the number of knots;
  
  # Pre-setting
  n <- ncol(yraw)
  L <- nrow(yraw)
  List <- pre_setting(argvals = argvals, d1 = d1, d2 = d2, K.p = K.p)
  B <- List$B
  Bder <- List$Bder
  Index.miss <- is.na(yraw)
  totalmiss <- 0
  
  # Sparse
  if (sum(Index.miss) > 0) {
    num.miss <- colSums(is.na(yraw))
    totalmiss <- mean(Index.miss)
    for (i in 1:n) {
      if (num.miss[i] > 0) {
        y <- yraw[, i]
        seq <- (1:L)[!is.na(y)] 
        seq2 <- (1:L)[is.na(y)] #na
        t1 <- argvals[seq]
        t2 <- argvals[seq2] #na
        fit <- c(JOPS::psNormal(x=t1, y=y[seq], xl=min(argvals), 
                                xr=max(argvals), xgrid=t2)$ygrid)
        yraw[seq2,i] <- fit
      }
    }
  }
  
  # Mean centered
  muRes <- muder_fun(yraw = yraw, argvals = argvals)
  mu_hat <- muRes$mu_hat
  Ycentered <- sweep(yraw, 1, mu_hat, "-")
  muder_hat <- muRes$muder_hat
  
  # Theta
  ThetaRes <- Theta_fun(Y = t(Ycentered), argvals = argvals, 
                        d1 = d1, d2 = d2, K.p = K.p, totalmiss = totalmiss)
  Theta_mat <- ThetaRes$Theta
  A0 <- ThetaRes$A0
  Sigmas <- ThetaRes$Sigmas
  
  # Estimate derivative-based eigfun and eigval
  EigRes <- dfpca_eigen_fun(argvals = argvals, Bder = Bder, Theta = Theta_mat)
  Phi <- EigRes$eigfuns
  lambda <- EigRes$eigvals
  # Estimate derivative-based scores
  # Smoothed curves
  A0Sigmas <- A0 %*% Sigmas %*% t(A0)
  Yder_hat <- Bder %*% A0Sigmas %*% t(B) %*% Ycentered
  
  # 'raw' derivative covariance 
  Kder_raw <- Yder_hat %*% t(Yder_hat) /n 
  # 'smoothed' derivative covariance
  Kder_smooth <- Bder %*% ThetaRes$Theta %*% t(Bder)
  gap <- diag(Kder_raw) - diag(Kder_smooth)
  lambda_xi <- 0
  
  # Derivative-based scores
  if (length(lambda) == 1) {
    #note that Phi and lambda are derivatived-based
    Phi <- matrix(Phi, ncol=1)
    tmp <- crossprod(Phi) + lambda_xi/lambda
  } else{
    tmp <- crossprod(Phi) + lambda_xi*diag(1/lambda)
  }
  xi <- solve(tmp) %*% t(Phi) %*% Yder_hat
  
  # Derivative-based refitted curves
  Refit <- Phi %*% xi + muder_hat
  
  return(list(mu = muder_hat, Phi = Phi, lambda = lambda, lambda_xi=lambda_xi,
              FVE = EigRes$FVE, xi = t(xi), Refit = Refit,
              rho = ThetaRes$rho, w = ThetaRes$w, gcv = ThetaRes$gcv))
}



# new2 ----
# lambda_xi=0.5*hat{lambda_xi} 
deriv_fpca_fun_new2 <- function(yraw, argvals, d1, d2, K.p){
  # yraw: Lxn observations
  # argvals: a vector containing observed locations on the functional domain
  # d1/d2: difference penalty order
  # K.p: the number of knots;
  
  # Pre-setting
  n <- ncol(yraw)
  L <- nrow(yraw)
  List <- pre_setting(argvals = argvals, d1 = d1, d2 = d2, K.p = K.p)
  B <- List$B
  Bder <- List$Bder
  Index.miss <- is.na(yraw)
  totalmiss <- 0
  
  # Sparse
  if (sum(Index.miss) > 0) {
    num.miss <- colSums(is.na(yraw))
    totalmiss <- mean(Index.miss)
    for (i in 1:n) {
      if (num.miss[i] > 0) {
        y <- yraw[, i]
        seq <- (1:L)[!is.na(y)] 
        seq2 <- (1:L)[is.na(y)] #na
        t1 <- argvals[seq]
        t2 <- argvals[seq2] #na
        fit <- c(JOPS::psNormal(x=t1, y=y[seq], xl=min(argvals), 
                                xr=max(argvals), xgrid=t2)$ygrid)
        yraw[seq2,i] <- fit
      }
    }
  }
  
  # Mean centered
  muRes <- muder_fun(yraw = yraw, argvals = argvals)
  mu_hat <- muRes$mu_hat
  Ycentered <- sweep(yraw, 1, mu_hat, "-")
  muder_hat <- muRes$muder_hat
  
  # Theta
  ThetaRes <- Theta_fun(Y = t(Ycentered), argvals = argvals, 
                        d1 = d1, d2 = d2, K.p = K.p, totalmiss = totalmiss)
  Theta_mat <- ThetaRes$Theta
  A0 <- ThetaRes$A0
  Sigmas <- ThetaRes$Sigmas
  
  # Estimate derivative-based eigfun and eigval
  EigRes <- dfpca_eigen_fun(argvals = argvals, Bder = Bder, Theta = Theta_mat)
  Phi <- EigRes$eigfuns
  lambda <- EigRes$eigvals
  # Estimate derivative-based scores
  # Smoothed curves
  A0Sigmas <- A0 %*% Sigmas %*% t(A0)
  Yder_hat <- Bder %*% A0Sigmas %*% t(B) %*% Ycentered
  
  # 'raw' derivative covariance 
  Kder_raw <- Yder_hat %*% t(Yder_hat) /n 
  # 'smoothed' derivative covariance
  Kder_smooth <- Bder %*% ThetaRes$Theta %*% t(Bder)
  gap <- diag(Kder_raw) - diag(Kder_smooth)
  hat_lambda_xi <- fdapace:::trapzRcpp(argvals, pmax(gap, 0))
  lambda_xi <- 0.5*hat_lambda_xi
  
  # Derivative-based scores
  if (length(lambda) == 1) {
    #note that Phi and lambda are derivatived-based
    Phi <- matrix(Phi, ncol=1)
    tmp <- crossprod(Phi) + lambda_xi/lambda
  } else{
    tmp <- crossprod(Phi) + lambda_xi*diag(1/lambda)
  }
  xi <- solve(tmp) %*% t(Phi) %*% Yder_hat
  
  # Derivative-based refitted curves
  Refit <- Phi %*% xi + muder_hat
  
  return(list(mu = muder_hat, Phi = Phi, lambda = lambda, lambda_xi=lambda_xi,
              FVE = EigRes$FVE, xi = t(xi), Refit = Refit,
              rho = ThetaRes$rho, w = ThetaRes$w, gcv = ThetaRes$gcv))
}



# new3 ----
# lambda_xi=2*hat{lambda_xi} 
deriv_fpca_fun_new3 <- function(yraw, argvals, d1, d2, K.p){
  # yraw: Lxn observations
  # argvals: a vector containing observed locations on the functional domain
  # d1/d2: difference penalty order
  # K.p: the number of knots;
  
  # Pre-setting
  n <- ncol(yraw)
  L <- nrow(yraw)
  List <- pre_setting(argvals = argvals, d1 = d1, d2 = d2, K.p = K.p)
  B <- List$B
  Bder <- List$Bder
  Index.miss <- is.na(yraw)
  totalmiss <- 0
  
  # Sparse
  if (sum(Index.miss) > 0) {
    num.miss <- colSums(is.na(yraw))
    totalmiss <- mean(Index.miss)
    for (i in 1:n) {
      if (num.miss[i] > 0) {
        y <- yraw[, i]
        seq <- (1:L)[!is.na(y)] 
        seq2 <- (1:L)[is.na(y)] #na
        t1 <- argvals[seq]
        t2 <- argvals[seq2] #na
        fit <- c(JOPS::psNormal(x=t1, y=y[seq], xl=min(argvals), 
                                xr=max(argvals), xgrid=t2)$ygrid)
        yraw[seq2,i] <- fit
      }
    }
  }
  
  # Mean centered
  muRes <- muder_fun(yraw = yraw, argvals = argvals)
  mu_hat <- muRes$mu_hat
  Ycentered <- sweep(yraw, 1, mu_hat, "-")
  muder_hat <- muRes$muder_hat
  
  # Theta
  ThetaRes <- Theta_fun(Y = t(Ycentered), argvals = argvals, 
                        d1 = d1, d2 = d2, K.p = K.p, totalmiss = totalmiss)
  Theta_mat <- ThetaRes$Theta
  A0 <- ThetaRes$A0
  Sigmas <- ThetaRes$Sigmas
  
  # Estimate derivative-based eigfun and eigval
  EigRes <- dfpca_eigen_fun(argvals = argvals, Bder = Bder, Theta = Theta_mat)
  Phi <- EigRes$eigfuns
  lambda <- EigRes$eigvals
  # Estimate derivative-based scores
  # Smoothed curves
  A0Sigmas <- A0 %*% Sigmas %*% t(A0)
  Yder_hat <- Bder %*% A0Sigmas %*% t(B) %*% Ycentered
  
  # 'raw' derivative covariance 
  Kder_raw <- Yder_hat %*% t(Yder_hat) /n 
  # 'smoothed' derivative covariance
  Kder_smooth <- Bder %*% ThetaRes$Theta %*% t(Bder)
  gap <- diag(Kder_raw) - diag(Kder_smooth)
  hat_lambda_xi <- fdapace:::trapzRcpp(argvals, pmax(gap, 0))
  lambda_xi <- 2*hat_lambda_xi
  
  # Derivative-based scores
  if (length(lambda) == 1) {
    #note that Phi and lambda are derivatived-based
    Phi <- matrix(Phi, ncol=1)
    tmp <- crossprod(Phi) + lambda_xi/lambda
  } else{
    tmp <- crossprod(Phi) + lambda_xi*diag(1/lambda)
  }
  xi <- solve(tmp) %*% t(Phi) %*% Yder_hat
  
  # Derivative-based refitted curves
  Refit <- Phi %*% xi + muder_hat
  
  return(list(mu = muder_hat, Phi = Phi, lambda = lambda, lambda_xi=lambda_xi,
              FVE = EigRes$FVE, xi = t(xi), Refit = Refit,
              rho = ThetaRes$rho, w = ThetaRes$w, gcv = ThetaRes$gcv))
}


# Simulation X1 ----
n_obs <- 100
n_points <- 101
argvals <- seq(0, 1, length.out = n_points)
noise_variance <- 0.25
min_obs <- 50
max_obs <- 60
d1 = 3
d2 = 2
K.p = 35
Nrep <- 500
cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)
MSE_list <- MSE_list_new1 <- MSE_list_new2 <- MSE_list_new3 <- list()


## No noise ----
Reg_lambdaxi_case1 <- foreach(rep=1:Nrep, .packages = c("fda", "fdapace", "funData", "mgcv")) %dopar% {
  # Data generation 
  set.seed(rep)
  fdata <- dgp(n_obs, argvals, noise_variance, min_obs, max_obs)
  
  # Ground truth
  Lt <- rep(list(argvals), n_obs)
  LDy1 <- split(fdata$data_der[[1]]@X, 1:n_obs)
  FPCA_Dy1 <- FPCA(Ly = LDy1, Lt = Lt, optns = list("FVEthreshold" = 0.95))
  
  # Different levels of regularization parameter
  y1_dfpcaRes <-  deriv_fpca_fun(yraw = t(fdata$data[[1]]@X), argvals = argvals, 
                                 d1 = d1, d2 = d2, K.p = K.p)
  y1_dfpcaRes_new1 <-  deriv_fpca_fun_new1(yraw = t(fdata$data[[1]]@X), argvals = argvals, 
                                 d1 = d1, d2 = d2, K.p = K.p)
  y1_dfpcaRes_new2 <-  deriv_fpca_fun_new2(yraw = t(fdata$data[[1]]@X), argvals = argvals, 
                                           d1 = d1, d2 = d2, K.p = K.p)
  y1_dfpcaRes_new3 <-  deriv_fpca_fun_new3(yraw = t(fdata$data[[1]]@X), argvals = argvals, 
                                           d1 = d1, d2 = d2, K.p = K.p)
  
  K1 <- min(ncol(FPCA_Dy1$phi), ncol(y1_dfpcaRes$Phi))
  
  # Evaluation
  MSE_list <- MSE_fun(FPCA_Dy1$xiEst[,1:K1], y1_dfpcaRes$xi[,1:K1])
  MSE_list_new1 <- MSE_fun(FPCA_Dy1$xiEst[,1:K1], y1_dfpcaRes_new1$xi[,1:K1])
  MSE_list_new2 <- MSE_fun(FPCA_Dy1$xiEst[,1:K1], y1_dfpcaRes_new2$xi[,1:K1])
  MSE_list_new3 <- MSE_fun(FPCA_Dy1$xiEst[,1:K1], y1_dfpcaRes_new3$xi[,1:K1])
  
  list(MSE_list = MSE_list, MSE_list_new1 = MSE_list_new1, 
       MSE_list_new2 = MSE_list_new2, MSE_list_new3 = MSE_list_new3)
}
saveRDS(Reg_lambdaxi_case1, file = "./SimRes/Reg_lambdaxi_case1.rds")




## Noisy ----
Reg_lambdaxi_case2 <- foreach(rep=1:Nrep, .packages = c("fda", "fdapace", "funData", "mgcv")) %dopar% {
  # Data generation 
  set.seed(rep)
  fdata <- dgp(n_obs, argvals, noise_variance, min_obs, max_obs)
  
  # Ground truth
  Lt <- rep(list(argvals), n_obs)
  LDy1 <- split(fdata$data_der[[1]]@X, 1:n_obs)
  FPCA_Dy1 <- FPCA(Ly = LDy1, Lt = Lt, optns = list("FVEthreshold" = 0.95))
  
  # Different levels of regularization parameter
  y1_dfpcaRes <-  deriv_fpca_fun(yraw = t(fdata$data_noisy[[1]]@X), argvals = argvals, 
                                 d1 = d1, d2 = d2, K.p = K.p)
  y1_dfpcaRes_new1 <-  deriv_fpca_fun_new1(yraw = t(fdata$data_noisy[[1]]@X), argvals = argvals, 
                                           d1 = d1, d2 = d2, K.p = K.p)
  y1_dfpcaRes_new2 <-  deriv_fpca_fun_new2(yraw = t(fdata$data_noisy[[1]]@X), argvals = argvals, 
                                           d1 = d1, d2 = d2, K.p = K.p)
  y1_dfpcaRes_new3 <-  deriv_fpca_fun_new3(yraw = t(fdata$data_noisy[[1]]@X), argvals = argvals, 
                                           d1 = d1, d2 = d2, K.p = K.p)
  
  K1 <- min(ncol(FPCA_Dy1$phi), ncol(y1_dfpcaRes$Phi))
  
  # Evaluation
  MSE_list <- MSE_fun(FPCA_Dy1$xiEst[,1:K1], y1_dfpcaRes$xi[,1:K1])
  MSE_list_new1 <- MSE_fun(FPCA_Dy1$xiEst[,1:K1], y1_dfpcaRes_new1$xi[,1:K1])
  MSE_list_new2 <- MSE_fun(FPCA_Dy1$xiEst[,1:K1], y1_dfpcaRes_new2$xi[,1:K1])
  MSE_list_new3 <- MSE_fun(FPCA_Dy1$xiEst[,1:K1], y1_dfpcaRes_new3$xi[,1:K1])
  
  list(MSE_list = MSE_list, MSE_list_new1 = MSE_list_new1, 
       MSE_list_new2 = MSE_list_new2, MSE_list_new3 = MSE_list_new3)
}
saveRDS(Reg_lambdaxi_case2, file = "./SimRes/Reg_lambdaxi_case2.rds")



## Sparse ----
Reg_lambdaxi_case3 <- foreach(rep=1:Nrep, .packages = c("fda", "fdapace", "funData", "mgcv")) %dopar% {
  # Data generation 
  set.seed(rep)
  fdata <- dgp(n_obs, argvals, noise_variance, min_obs, max_obs)
  
  # Ground truth
  Lt <- rep(list(argvals), n_obs)
  LDy1 <- split(fdata$data_der[[1]]@X, 1:n_obs)
  FPCA_Dy1 <- FPCA(Ly = LDy1, Lt = Lt, optns = list("FVEthreshold" = 0.95))
  
  # Different levels of regularization parameter
  y1_dfpcaRes <-  deriv_fpca_fun(yraw = t(fdata$data_sparse[[1]]@X), argvals = argvals, 
                                 d1 = d1, d2 = d2, K.p = K.p)
  y1_dfpcaRes_new1 <-  deriv_fpca_fun_new1(yraw = t(fdata$data_sparse[[1]]@X), argvals = argvals, 
                                           d1 = d1, d2 = d2, K.p = K.p)
  y1_dfpcaRes_new2 <-  deriv_fpca_fun_new2(yraw = t(fdata$data_sparse[[1]]@X), argvals = argvals, 
                                           d1 = d1, d2 = d2, K.p = K.p)
  y1_dfpcaRes_new3 <-  deriv_fpca_fun_new3(yraw = t(fdata$data_sparse[[1]]@X), argvals = argvals, 
                                           d1 = d1, d2 = d2, K.p = K.p)
  
  K1 <- min(ncol(FPCA_Dy1$phi), ncol(y1_dfpcaRes$Phi))
  
  # Evaluation
  MSE_list <- MSE_fun(FPCA_Dy1$xiEst[,1:K1], y1_dfpcaRes$xi[,1:K1])
  MSE_list_new1 <- MSE_fun(FPCA_Dy1$xiEst[,1:K1], y1_dfpcaRes_new1$xi[,1:K1])
  MSE_list_new2 <- MSE_fun(FPCA_Dy1$xiEst[,1:K1], y1_dfpcaRes_new2$xi[,1:K1])
  MSE_list_new3 <- MSE_fun(FPCA_Dy1$xiEst[,1:K1], y1_dfpcaRes_new3$xi[,1:K1])
  
  list(MSE_list = MSE_list, MSE_list_new1 = MSE_list_new1, 
       MSE_list_new2 = MSE_list_new2, MSE_list_new3 = MSE_list_new3)
}
saveRDS(Reg_lambdaxi_case3, file = "./SimRes/Reg_lambdaxi_case3.rds")


# Simulation X2 ----

## No noise ----
y2_Reg_lambdaxi_case1 <- foreach(rep=1:Nrep, .packages = c("fda", "fdapace", "funData", "mgcv")) %dopar% {
  # Data generation 
  set.seed(rep)
  fdata <- dgp(n_obs, argvals, noise_variance, min_obs, max_obs)
  
  # Ground truth
  Lt <- rep(list(argvals), n_obs)
  LDy1 <- split(fdata$data_der[[2]]@X, 1:n_obs)
  FPCA_Dy1 <- FPCA(Ly = LDy1, Lt = Lt, optns = list("FVEthreshold" = 0.95))
  
  # Different levels of regularization parameter
  y1_dfpcaRes <-  deriv_fpca_fun(yraw = t(fdata$data[[2]]@X), argvals = argvals, 
                                 d1 = d1, d2 = d2, K.p = K.p)
  y1_dfpcaRes_new1 <-  deriv_fpca_fun_new1(yraw = t(fdata$data[[2]]@X), argvals = argvals, 
                                           d1 = d1, d2 = d2, K.p = K.p)
  y1_dfpcaRes_new2 <-  deriv_fpca_fun_new2(yraw = t(fdata$data[[2]]@X), argvals = argvals, 
                                           d1 = d1, d2 = d2, K.p = K.p)
  y1_dfpcaRes_new3 <-  deriv_fpca_fun_new3(yraw = t(fdata$data[[2]]@X), argvals = argvals, 
                                           d1 = d1, d2 = d2, K.p = K.p)
  
  K1 <- min(ncol(FPCA_Dy1$phi), ncol(y1_dfpcaRes$Phi))
  
  # Evaluation
  MSE_list <- MSE_fun(FPCA_Dy1$xiEst[,1:K1], y1_dfpcaRes$xi[,1:K1])
  MSE_list_new1 <- MSE_fun(FPCA_Dy1$xiEst[,1:K1], y1_dfpcaRes_new1$xi[,1:K1])
  MSE_list_new2 <- MSE_fun(FPCA_Dy1$xiEst[,1:K1], y1_dfpcaRes_new2$xi[,1:K1])
  MSE_list_new3 <- MSE_fun(FPCA_Dy1$xiEst[,1:K1], y1_dfpcaRes_new3$xi[,1:K1])
  
  list(MSE_list = MSE_list, MSE_list_new1 = MSE_list_new1, 
       MSE_list_new2 = MSE_list_new2, MSE_list_new3 = MSE_list_new3)
}
saveRDS(y2_Reg_lambdaxi_case1, file = "./SimRes/y2_Reg_lambdaxi_case1.rds")




## Noisy ----
y2_Reg_lambdaxi_case2 <- foreach(rep=1:Nrep, .packages = c("fda", "fdapace", "funData", "mgcv")) %dopar% {
  # Data generation 
  set.seed(rep)
  fdata <- dgp(n_obs, argvals, noise_variance, min_obs, max_obs)
  
  # Ground truth
  Lt <- rep(list(argvals), n_obs)
  LDy1 <- split(fdata$data_der[[2]]@X, 1:n_obs)
  FPCA_Dy1 <- FPCA(Ly = LDy1, Lt = Lt, optns = list("FVEthreshold" = 0.95))
  
  # Different levels of regularization parameter
  y1_dfpcaRes <-  deriv_fpca_fun(yraw = t(fdata$data_noisy[[2]]@X), argvals = argvals, 
                                 d1 = d1, d2 = d2, K.p = K.p)
  y1_dfpcaRes_new1 <-  deriv_fpca_fun_new1(yraw = t(fdata$data_noisy[[2]]@X), argvals = argvals, 
                                           d1 = d1, d2 = d2, K.p = K.p)
  y1_dfpcaRes_new2 <-  deriv_fpca_fun_new2(yraw = t(fdata$data_noisy[[2]]@X), argvals = argvals, 
                                           d1 = d1, d2 = d2, K.p = K.p)
  y1_dfpcaRes_new3 <-  deriv_fpca_fun_new3(yraw = t(fdata$data_noisy[[2]]@X), argvals = argvals, 
                                           d1 = d1, d2 = d2, K.p = K.p)
  
  K1 <- min(ncol(FPCA_Dy1$phi), ncol(y1_dfpcaRes$Phi))
  
  # Evaluation
  MSE_list <- MSE_fun(FPCA_Dy1$xiEst[,1:K1], y1_dfpcaRes$xi[,1:K1])
  MSE_list_new1 <- MSE_fun(FPCA_Dy1$xiEst[,1:K1], y1_dfpcaRes_new1$xi[,1:K1])
  MSE_list_new2 <- MSE_fun(FPCA_Dy1$xiEst[,1:K1], y1_dfpcaRes_new2$xi[,1:K1])
  MSE_list_new3 <- MSE_fun(FPCA_Dy1$xiEst[,1:K1], y1_dfpcaRes_new3$xi[,1:K1])
  
  list(MSE_list = MSE_list, MSE_list_new1 = MSE_list_new1, 
       MSE_list_new2 = MSE_list_new2, MSE_list_new3 = MSE_list_new3)
}
saveRDS(y2_Reg_lambdaxi_case2, file = "./SimRes/y2_Reg_lambdaxi_case2.rds")



## Sparse ----
y2_Reg_lambdaxi_case3 <- foreach(rep=1:Nrep, .packages = c("fda", "fdapace", "funData", "mgcv")) %dopar% {
  # Data generation 
  set.seed(rep)
  fdata <- dgp(n_obs, argvals, noise_variance, min_obs, max_obs)
  
  # Ground truth
  Lt <- rep(list(argvals), n_obs)
  LDy1 <- split(fdata$data_der[[2]]@X, 1:n_obs)
  FPCA_Dy1 <- FPCA(Ly = LDy1, Lt = Lt, optns = list("FVEthreshold" = 0.95))
  
  # Different levels of regularization parameter
  y1_dfpcaRes <-  deriv_fpca_fun(yraw = t(fdata$data_sparse[[2]]@X), argvals = argvals, 
                                 d1 = d1, d2 = d2, K.p = K.p)
  y1_dfpcaRes_new1 <-  deriv_fpca_fun_new1(yraw = t(fdata$data_sparse[[2]]@X), argvals = argvals, 
                                           d1 = d1, d2 = d2, K.p = K.p)
  y1_dfpcaRes_new2 <-  deriv_fpca_fun_new2(yraw = t(fdata$data_sparse[[2]]@X), argvals = argvals, 
                                           d1 = d1, d2 = d2, K.p = K.p)
  y1_dfpcaRes_new3 <-  deriv_fpca_fun_new3(yraw = t(fdata$data_sparse[[2]]@X), argvals = argvals, 
                                           d1 = d1, d2 = d2, K.p = K.p)
  
  K1 <- min(ncol(FPCA_Dy1$phi), ncol(y1_dfpcaRes$Phi))
  
  # Evaluation
  MSE_list <- MSE_fun(FPCA_Dy1$xiEst[,1:K1], y1_dfpcaRes$xi[,1:K1])
  MSE_list_new1 <- MSE_fun(FPCA_Dy1$xiEst[,1:K1], y1_dfpcaRes_new1$xi[,1:K1])
  MSE_list_new2 <- MSE_fun(FPCA_Dy1$xiEst[,1:K1], y1_dfpcaRes_new2$xi[,1:K1])
  MSE_list_new3 <- MSE_fun(FPCA_Dy1$xiEst[,1:K1], y1_dfpcaRes_new3$xi[,1:K1])
  
  list(MSE_list = MSE_list, MSE_list_new1 = MSE_list_new1, 
       MSE_list_new2 = MSE_list_new2, MSE_list_new3 = MSE_list_new3)
}
saveRDS(y2_Reg_lambdaxi_case3, file = "./SimRes/y2_Reg_lambdaxi_case3.rds")
