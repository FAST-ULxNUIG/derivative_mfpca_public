################################################################################
# Selection of penalty orders (l1, l2)
################################################################################

# Simulation setting: (l1,l2)={(1,2), (1,3), (2,3)}
# Aim: illustrate the selection of (l1, l2) will not affect the final results


library(fda)
library(funData)
library(fdapace)
library(refund)
library(mgcv)
library(foreach)
library(doParallel)
source("01_dgp.R")
source("02_method.R")
source("03_errors.R")


# Pre-setting ----
n_obs <- 100
n_points <- 101
argvals <- seq(0, 1, length.out = n_points)
noise_variance <- 0.25
min_obs <- 50
max_obs <- 60
K.p <- 35


# Simulation ----
Nrep <- 500
cl <- parallel::makeCluster(6)
doParallel::registerDoParallel(cl)

# No noise ----
## d1=1 d2=2 ----
d1=1; d2=2
PenOrder_set1_nonoise <- foreach(rep=1:Nrep, .packages = c("fda", "fdapace", "funData", "mgcv")) %dopar% {
  # Data generation 
  set.seed(rep)
  fdata <- dgp(n_obs, argvals, noise_variance, min_obs, max_obs)
  
  # Ground truth
  Lt <- rep(list(argvals), n_obs)
  LDy1 <- split(fdata$data_der[[1]]@X, 1:n_obs)
  FPCA_Dy1 <- FPCA(Ly = LDy1, Lt = Lt, optns = list("FVEthreshold" = 0.95))
  
  # Additive penalty
  y1_dfpcaRes <-  deriv_fpca_fun(yraw = t(fdata$data[[1]]@X), argvals = argvals, 
                                 d1 = d1, d2 = d2, K.p = K.p)
  
  K1 <- min(ncol(FPCA_Dy1$phi), ncol(y1_dfpcaRes$Phi))
  
  # Refit
  Refit_DFPCA <- 
    tcrossprod(y1_dfpcaRes$Phi[,1:K1], y1_dfpcaRes$xi[,1:K1])+y1_dfpcaRes$mu
  
  # Evaluation
  ISE_list <- ISE_fun(argvals, FPCA_Dy1$phi[,1:K1], y1_dfpcaRes$Phi[,1:K1])
  RE_list <- RE_fun(FPCA_Dy1$lambda[1:K1], y1_dfpcaRes$lambda[1:K1])
  MSE_list <- MSE_fun(FPCA_Dy1$xiEst[,1:K1], y1_dfpcaRes$xi[,1:K1])
  RMISE_list <- RMISE_fun(argvals, t(fdata$data_der[[1]]@X), Refit_DFPCA)
  
  list(RE_list = RE_list, ISE_list = ISE_list, MSE_list = MSE_list, RMISE_list = RMISE_list)
}
saveRDS(PenOrder_set1_nonoise, file = "./SimRes/PenOrder_set1_nonoise.rds")


## d1=1 d2=3 ----
d1=1; d2=3
PenOrder_set2_nonoise <- foreach(rep=1:Nrep, .packages = c("fda", "fdapace", "funData", "mgcv")) %dopar% {
  # Data generation 
  set.seed(rep)
  fdata <- dgp(n_obs, argvals, noise_variance, min_obs, max_obs)
  
  # Ground truth
  Lt <- rep(list(argvals), n_obs)
  LDy1 <- split(fdata$data_der[[1]]@X, 1:n_obs)
  FPCA_Dy1 <- FPCA(Ly = LDy1, Lt = Lt, optns = list("FVEthreshold" = 0.95))
  
  # Additive penalty
  y1_dfpcaRes <-  deriv_fpca_fun(yraw = t(fdata$data[[1]]@X), argvals = argvals, 
                                 d1 = d1, d2 = d2, K.p = K.p)
  
  K1 <- min(ncol(FPCA_Dy1$phi), ncol(y1_dfpcaRes$Phi))
  
  # Refit
  Refit_DFPCA <- 
    tcrossprod(y1_dfpcaRes$Phi[,1:K1], y1_dfpcaRes$xi[,1:K1])+y1_dfpcaRes$mu
  
  # Evaluation
  ISE_list <- ISE_fun(argvals, FPCA_Dy1$phi[,1:K1], y1_dfpcaRes$Phi[,1:K1])
  RE_list <- RE_fun(FPCA_Dy1$lambda[1:K1], y1_dfpcaRes$lambda[1:K1])
  MSE_list <- MSE_fun(FPCA_Dy1$xiEst[,1:K1], y1_dfpcaRes$xi[,1:K1])
  RMISE_list <- RMISE_fun(argvals, t(fdata$data_der[[1]]@X), Refit_DFPCA)
  
  list(RE_list = RE_list, ISE_list = ISE_list, MSE_list = MSE_list, RMISE_list = RMISE_list)
}
saveRDS(PenOrder_set2_nonoise, file = "./SimRes/PenOrder_set2_nonoise.rds")



## d1=2 d2=3 ----
d1=2; d2=3
PenOrder_set3_nonoise <- foreach(rep=1:Nrep, .packages = c("fda", "fdapace", "funData", "mgcv")) %dopar% {
  # Data generation 
  set.seed(rep)
  fdata <- dgp(n_obs, argvals, noise_variance, min_obs, max_obs)
  
  # Ground truth
  Lt <- rep(list(argvals), n_obs)
  LDy1 <- split(fdata$data_der[[1]]@X, 1:n_obs)
  FPCA_Dy1 <- FPCA(Ly = LDy1, Lt = Lt, optns = list("FVEthreshold" = 0.95))
  
  # Additive penalty
  y1_dfpcaRes <-  deriv_fpca_fun(yraw = t(fdata$data[[1]]@X), argvals = argvals, 
                                 d1 = d1, d2 = d2, K.p = K.p)
  
  K1 <- min(ncol(FPCA_Dy1$phi), ncol(y1_dfpcaRes$Phi))
  
  # Refit
  Refit_DFPCA <- 
    tcrossprod(y1_dfpcaRes$Phi[,1:K1], y1_dfpcaRes$xi[,1:K1])+y1_dfpcaRes$mu
  
  # Evaluation
  ISE_list <- ISE_fun(argvals, FPCA_Dy1$phi[,1:K1], y1_dfpcaRes$Phi[,1:K1])
  RE_list <- RE_fun(FPCA_Dy1$lambda[1:K1], y1_dfpcaRes$lambda[1:K1])
  MSE_list <- MSE_fun(FPCA_Dy1$xiEst[,1:K1], y1_dfpcaRes$xi[,1:K1])
  RMISE_list <- RMISE_fun(argvals, t(fdata$data_der[[1]]@X), Refit_DFPCA)
  
  list(RE_list = RE_list, ISE_list = ISE_list, MSE_list = MSE_list, RMISE_list = RMISE_list)
}
saveRDS(PenOrder_set3_nonoise, file = "./SimRes/PenOrder_set3_nonoise.rds")





# Noisy ----
## d1=1 d2=2 ----
d1=1; d2=2
PenOrder_set1_noisy <- foreach(rep=1:Nrep, .packages = c("fda", "fdapace", "funData", "mgcv")) %dopar% {
  # Data generation 
  set.seed(rep)
  fdata <- dgp(n_obs, argvals, noise_variance, min_obs, max_obs)
  
  # Ground truth
  Lt <- rep(list(argvals), n_obs)
  LDy1 <- split(fdata$data_der[[1]]@X, 1:n_obs)
  FPCA_Dy1 <- FPCA(Ly = LDy1, Lt = Lt, optns = list("FVEthreshold" = 0.95))
  
  # Additive penalty
  y1_dfpcaRes <-  deriv_fpca_fun(yraw = t(fdata$data_noisy[[1]]@X), argvals = argvals, 
                                 d1 = d1, d2 = d2, K.p = K.p)
  
  K1 <- min(ncol(FPCA_Dy1$phi), ncol(y1_dfpcaRes$Phi))
  
  # Refit
  Refit_DFPCA <- 
    tcrossprod(y1_dfpcaRes$Phi[,1:K1], y1_dfpcaRes$xi[,1:K1])+y1_dfpcaRes$mu
  
  # Evaluation
  ISE_list <- ISE_fun(argvals, FPCA_Dy1$phi[,1:K1], y1_dfpcaRes$Phi[,1:K1])
  RE_list <- RE_fun(FPCA_Dy1$lambda[1:K1], y1_dfpcaRes$lambda[1:K1])
  MSE_list <- MSE_fun(FPCA_Dy1$xiEst[,1:K1], y1_dfpcaRes$xi[,1:K1])
  RMISE_list <- RMISE_fun(argvals, t(fdata$data_der[[1]]@X), Refit_DFPCA)
  
  list(RE_list = RE_list, ISE_list = ISE_list, MSE_list = MSE_list, RMISE_list = RMISE_list)
}
saveRDS(PenOrder_set1_noisy, file = "./SimRes/PenOrder_set1_noisy.rds")


## d1=1 d2=3 ----
d1=1; d2=3
PenOrder_set2_noisy <- foreach(rep=1:Nrep, .packages = c("fda", "fdapace", "funData", "mgcv")) %dopar% {
  # Data generation 
  set.seed(rep)
  fdata <- dgp(n_obs, argvals, noise_variance, min_obs, max_obs)
  
  # Ground truth
  Lt <- rep(list(argvals), n_obs)
  LDy1 <- split(fdata$data_der[[1]]@X, 1:n_obs)
  FPCA_Dy1 <- FPCA(Ly = LDy1, Lt = Lt, optns = list("FVEthreshold" = 0.95))
  
  # Additive penalty
  y1_dfpcaRes <-  deriv_fpca_fun(yraw = t(fdata$data_noisy[[1]]@X), argvals = argvals, 
                                 d1 = d1, d2 = d2, K.p = K.p)
  
  K1 <- min(ncol(FPCA_Dy1$phi), ncol(y1_dfpcaRes$Phi))
  
  # Refit
  Refit_DFPCA <- 
    tcrossprod(y1_dfpcaRes$Phi[,1:K1], y1_dfpcaRes$xi[,1:K1])+y1_dfpcaRes$mu
  
  # Evaluation
  ISE_list <- ISE_fun(argvals, FPCA_Dy1$phi[,1:K1], y1_dfpcaRes$Phi[,1:K1])
  RE_list <- RE_fun(FPCA_Dy1$lambda[1:K1], y1_dfpcaRes$lambda[1:K1])
  MSE_list <- MSE_fun(FPCA_Dy1$xiEst[,1:K1], y1_dfpcaRes$xi[,1:K1])
  RMISE_list <- RMISE_fun(argvals, t(fdata$data_der[[1]]@X), Refit_DFPCA)
  
  list(RE_list = RE_list, ISE_list = ISE_list, MSE_list = MSE_list, RMISE_list = RMISE_list)
}
saveRDS(PenOrder_set2_noisy, file = "./SimRes/PenOrder_set2_noisy.rds")



## d1=2 d2=3 ----
d1=2; d2=3
PenOrder_set3_noisy <- foreach(rep=1:Nrep, .packages = c("fda", "fdapace", "funData", "mgcv")) %dopar% {
  # Data generation 
  set.seed(rep)
  fdata <- dgp(n_obs, argvals, noise_variance, min_obs, max_obs)
  
  # Ground truth
  Lt <- rep(list(argvals), n_obs)
  LDy1 <- split(fdata$data_der[[1]]@X, 1:n_obs)
  FPCA_Dy1 <- FPCA(Ly = LDy1, Lt = Lt, optns = list("FVEthreshold" = 0.95))
  
  # Additive penalty
  y1_dfpcaRes <-  deriv_fpca_fun(yraw = t(fdata$data_noisy[[1]]@X), argvals = argvals, 
                                 d1 = d1, d2 = d2, K.p = K.p)
  
  K1 <- min(ncol(FPCA_Dy1$phi), ncol(y1_dfpcaRes$Phi))
  
  # Refit
  Refit_DFPCA <- 
    tcrossprod(y1_dfpcaRes$Phi[,1:K1], y1_dfpcaRes$xi[,1:K1])+y1_dfpcaRes$mu
  
  # Evaluation
  ISE_list <- ISE_fun(argvals, FPCA_Dy1$phi[,1:K1], y1_dfpcaRes$Phi[,1:K1])
  RE_list <- RE_fun(FPCA_Dy1$lambda[1:K1], y1_dfpcaRes$lambda[1:K1])
  MSE_list <- MSE_fun(FPCA_Dy1$xiEst[,1:K1], y1_dfpcaRes$xi[,1:K1])
  RMISE_list <- RMISE_fun(argvals, t(fdata$data_der[[1]]@X), Refit_DFPCA)
  
  list(RE_list = RE_list, ISE_list = ISE_list, MSE_list = MSE_list, RMISE_list = RMISE_list)
}
saveRDS(PenOrder_set3_noisy, file = "./SimRes/PenOrder_set3_noisy.rds")




# Sparse ----
## d1=1 d2=2 ----
d1=1; d2=2
PenOrder_set1_sparse <- foreach(rep=1:Nrep, .packages = c("fda", "fdapace", "funData", "mgcv")) %dopar% {
  # Data generation 
  set.seed(rep)
  fdata <- dgp(n_obs, argvals, noise_variance, min_obs, max_obs)
  
  # Ground truth
  Lt <- rep(list(argvals), n_obs)
  LDy1 <- split(fdata$data_der[[1]]@X, 1:n_obs)
  FPCA_Dy1 <- FPCA(Ly = LDy1, Lt = Lt, optns = list("FVEthreshold" = 0.95))
  
  # Additive penalty
  y1_dfpcaRes <-  deriv_fpca_fun(yraw = t(fdata$data_sparse[[1]]@X), argvals = argvals, 
                                 d1 = d1, d2 = d2, K.p = K.p)
  
  K1 <- min(ncol(FPCA_Dy1$phi), ncol(y1_dfpcaRes$Phi))
  
  # Refit
  Refit_DFPCA <- 
    tcrossprod(y1_dfpcaRes$Phi[,1:K1], y1_dfpcaRes$xi[,1:K1])+y1_dfpcaRes$mu
  
  # Evaluation
  ISE_list <- ISE_fun(argvals, FPCA_Dy1$phi[,1:K1], y1_dfpcaRes$Phi[,1:K1])
  RE_list <- RE_fun(FPCA_Dy1$lambda[1:K1], y1_dfpcaRes$lambda[1:K1])
  MSE_list <- MSE_fun(FPCA_Dy1$xiEst[,1:K1], y1_dfpcaRes$xi[,1:K1])
  RMISE_list <- RMISE_fun(argvals, t(fdata$data_der[[1]]@X), Refit_DFPCA)
  
  list(RE_list = RE_list, ISE_list = ISE_list, MSE_list = MSE_list, RMISE_list = RMISE_list)
}
saveRDS(PenOrder_set1_sparse, file = "./SimRes/PenOrder_set1_sparse.rds")


## d1=1 d2=3 ----
d1=1; d2=3
PenOrder_set2_sparse <- foreach(rep=1:Nrep, .packages = c("fda", "fdapace", "funData", "mgcv")) %dopar% {
  # Data generation 
  set.seed(rep)
  fdata <- dgp(n_obs, argvals, noise_variance, min_obs, max_obs)
  
  # Ground truth
  Lt <- rep(list(argvals), n_obs)
  LDy1 <- split(fdata$data_der[[1]]@X, 1:n_obs)
  FPCA_Dy1 <- FPCA(Ly = LDy1, Lt = Lt, optns = list("FVEthreshold" = 0.95))
  
  # Additive penalty
  y1_dfpcaRes <-  deriv_fpca_fun(yraw = t(fdata$data_sparse[[1]]@X), argvals = argvals, 
                                 d1 = d1, d2 = d2, K.p = K.p)
  
  K1 <- min(ncol(FPCA_Dy1$phi), ncol(y1_dfpcaRes$Phi))
  
  # Refit
  Refit_DFPCA <- 
    tcrossprod(y1_dfpcaRes$Phi[,1:K1], y1_dfpcaRes$xi[,1:K1])+y1_dfpcaRes$mu
  
  # Evaluation
  ISE_list <- ISE_fun(argvals, FPCA_Dy1$phi[,1:K1], y1_dfpcaRes$Phi[,1:K1])
  RE_list <- RE_fun(FPCA_Dy1$lambda[1:K1], y1_dfpcaRes$lambda[1:K1])
  MSE_list <- MSE_fun(FPCA_Dy1$xiEst[,1:K1], y1_dfpcaRes$xi[,1:K1])
  RMISE_list <- RMISE_fun(argvals, t(fdata$data_der[[1]]@X), Refit_DFPCA)
  
  list(RE_list = RE_list, ISE_list = ISE_list, MSE_list = MSE_list, RMISE_list = RMISE_list)
}
saveRDS(PenOrder_set2_sparse, file = "./SimRes/PenOrder_set2_sparse.rds")



## d1=2 d2=3 ----
d1=2; d2=3
PenOrder_set3_sparse <- foreach(rep=1:Nrep, .packages = c("fda", "fdapace", "funData", "mgcv")) %dopar% {
  # Data generation 
  set.seed(rep)
  fdata <- dgp(n_obs, argvals, noise_variance, min_obs, max_obs)
  
  # Ground truth
  Lt <- rep(list(argvals), n_obs)
  LDy1 <- split(fdata$data_der[[1]]@X, 1:n_obs)
  FPCA_Dy1 <- FPCA(Ly = LDy1, Lt = Lt, optns = list("FVEthreshold" = 0.95))
  
  # Additive penalty
  y1_dfpcaRes <-  deriv_fpca_fun(yraw = t(fdata$data_sparse[[1]]@X), argvals = argvals, 
                                 d1 = d1, d2 = d2, K.p = K.p)
  
  K1 <- min(ncol(FPCA_Dy1$phi), ncol(y1_dfpcaRes$Phi))
  
  # Refit
  Refit_DFPCA <- 
    tcrossprod(y1_dfpcaRes$Phi[,1:K1], y1_dfpcaRes$xi[,1:K1])+y1_dfpcaRes$mu
  
  # Evaluation
  ISE_list <- ISE_fun(argvals, FPCA_Dy1$phi[,1:K1], y1_dfpcaRes$Phi[,1:K1])
  RE_list <- RE_fun(FPCA_Dy1$lambda[1:K1], y1_dfpcaRes$lambda[1:K1])
  MSE_list <- MSE_fun(FPCA_Dy1$xiEst[,1:K1], y1_dfpcaRes$xi[,1:K1])
  RMISE_list <- RMISE_fun(argvals, t(fdata$data_der[[1]]@X), Refit_DFPCA)
  
  list(RE_list = RE_list, ISE_list = ISE_list, MSE_list = MSE_list, RMISE_list = RMISE_list)
}
saveRDS(PenOrder_set3_sparse, file = "./SimRes/PenOrder_set3_sparse.rds")

