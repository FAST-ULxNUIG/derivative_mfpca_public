################################################################################
# Univariate Derivative-Based FPCA Simulation
# y1 and y2 under three settings: No noise; Noisy; Sparse
# Methods: Single penalty vs Additive penalty vs FPCAder
################################################################################

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
source("SinglePenalty.R")


# Pre-setting ----
n_obs <- 100
n_points <- 101
argvals <- seq(0, 1, length.out = n_points)
noise_variance <- 0.25
min_obs <- 50
max_obs <- 60
d1 = 3
d2 = 2
K.p = 35


# Simulation ----
Nrep <- 500
cl <- parallel::makeCluster(6)
doParallel::registerDoParallel(cl)
RE_list <- ISE_list <- RMISE_list <- MSE_list <- list()


# No noise ----
## y1: SP + no noise ----
start_time <- proc.time()
y1_case1_sp <- foreach(rep=1:Nrep, .packages = c("fda", "fdapace", "funData", "mgcv")) %dopar% {
  # Data generation 
  set.seed(rep)
  fdata <- dgp(n_obs, argvals, noise_variance, min_obs, max_obs)
  
  # Ground truth
  Lt <- rep(list(argvals), n_obs)
  LDy1 <- split(fdata$data_der[[1]]@X, 1:n_obs)
  FPCA_Dy1 <- FPCA(Ly = LDy1, Lt = Lt, optns = list("FVEthreshold" = 0.95))
  
  # Single penalty
  y1_dfpcaRes <-  deriv_fpca_fun_sp(yraw = t(fdata$data[[1]]@X), argvals = argvals, 
                                    d1 = d1, K.p = K.p)
  # Truncation
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
end_time <- proc.time()
time_case1_sp <- end_time - start_time
saveRDS(y1_case1_sp, file = "./SimRes/y1_case1_sp.rds")
saveRDS(time_case1_sp, file = "./SimRes/time_case1_sp.rds")

## y1: AP + no noise ----
start_time <- proc.time()
y1_case1_ap <- foreach(rep=1:Nrep, .packages = c("fda", "fdapace", "funData", "mgcv")) %dopar% {
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
end_time <- proc.time()
time_case1_ap <- end_time - start_time
saveRDS(y1_case1_ap, file = "./SimRes/y1_case1_ap.rds")
saveRDS(time_case1_ap, file = "./SimRes/time_case1_ap.rds")


## y1: FPCAder + no noise ----
start_time <- proc.time()
y1_case1_FPCAder <- foreach(rep=1:Nrep, .packages = c("fda", "fdapace", "funData", "mgcv")) %dopar% {
  # Data generation 
  set.seed(rep)
  fdata <- dgp(n_obs, argvals, noise_variance, min_obs, max_obs)
  
  # Ground truth
  Lt <- rep(list(argvals), n_obs)
  LDy1 <- split(fdata$data_der[[1]]@X, 1:n_obs)
  FPCA_Dy1 <- FPCA(Ly = LDy1, Lt = Lt, optns = list("FVEthreshold" = 0.95))
  
  # FPCAder
  Ly1 <- split(fdata$data[[1]]@X, 1:n_obs)
  fpcaObj1 <- FPCA(Ly = Ly1, Lt = Lt, optns = list("FVEthreshold" = 0.95))
  DPC1 <- FPCAder(fpcaObj1, derOptns = list(method='DPC'))
  
  K1 <- min(ncol(FPCA_Dy1$phi), ncol(DPC1$phiDer))
  
  muRes <- muder_fun(yraw = t(fdata$data[[1]]@X), argvals = argvals)
  muDer <- muRes$muder_hat
  Refit_DPC1 <- tcrossprod(DPC1$phiDer[,1:K1], DPC1$xiDer[,1:K1]) + muDer
  
  # Evaluation
  K1 <- min(ncol(FPCA_Dy1$phi), ncol(DPC1$phiDer))
  ISE_list <- ISE_fun(argvals, FPCA_Dy1$phi[,1:K1], DPC1$phiDer[,1:K1])
  RE_list <- RE_fun(FPCA_Dy1$lambda[1:K1], DPC1$lambdaDer[1:K1])
  MSE_list <- MSE_fun(FPCA_Dy1$xiEst[,1:K1], DPC1$xiDer[,1:K1])
  RMISE_list <- RMISE_fun(argvals, t(fdata$data_der[[1]]@X), Refit_DPC1)
  
  list(RE_list = RE_list, ISE_list = ISE_list, MSE_list = MSE_list, RMISE_list = RMISE_list)
}
end_time <- proc.time()
time_case1_fpcader <- end_time - start_time
saveRDS(y1_case1_FPCAder, file = "./SimRes/y1_case1_fpcader.rds")
saveRDS(time_case1_fpcader, file = "./SimRes/time_case1_fpcader.rds")


# Noisy ----
## y1: SP + noisy ----
start_time <- proc.time()
y1_case2_sp <- foreach(rep=1:Nrep, .packages = c("fda", "fdapace", "funData", "mgcv")) %dopar% {
  # Data generation 
  set.seed(rep)
  fdata <- dgp(n_obs, argvals, noise_variance, min_obs, max_obs)
  
  # Ground truth
  Lt <- rep(list(argvals), n_obs)
  LDy1 <- split(fdata$data_der[[1]]@X, 1:n_obs)
  FPCA_Dy1 <- FPCA(Ly = LDy1, Lt = Lt, optns = list("FVEthreshold" = 0.95))
  
  # Single penalty
  y1_dfpcaRes <-  deriv_fpca_fun_sp(yraw = t(fdata$data_noisy[[1]]@X), argvals = argvals, 
                                    d1 = d1, K.p = K.p)
  
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
end_time <- proc.time()
time_case2_sp <- end_time - start_time
saveRDS(y1_case2_sp, file = "./SimRes/y1_case2_sp.rds")
saveRDS(time_case2_sp, file = "./SimRes/time_case2_sp.rds")


## y1: AP + noisy ----
start_time <- proc.time()
y1_case2_ap <- foreach(rep=1:Nrep, .packages = c("fda", "fdapace", "funData", "mgcv")) %dopar% {
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
end_time <- proc.time()
time_case2_ap <- end_time - start_time
saveRDS(y1_case2_ap, file = "./SimRes/y1_case2_ap.rds")
saveRDS(time_case2_ap, file = "./SimRes/time_case2_ap.rds")



## y1: FPCAder + noisy ----
start_time <- proc.time()
y1_case2_FPCAder <- foreach(rep=1:Nrep, .packages = c("fda", "fdapace", "funData", "mgcv")) %dopar% {
  # Data generation 
  set.seed(rep)
  fdata <- dgp(n_obs, argvals, noise_variance, min_obs, max_obs)
  
  # Ground truth
  Lt <- rep(list(argvals), n_obs)
  LDy1 <- split(fdata$data_der[[1]]@X, 1:n_obs)
  FPCA_Dy1 <- FPCA(Ly = LDy1, Lt = Lt, optns = list("FVEthreshold" = 0.95))
  
  # FPCAder
  Ly1 <- split(fdata$data_noisy[[1]]@X, 1:n_obs)
  fpcaObj1 <- FPCA(Ly = Ly1, Lt = Lt, 
                   optns = list("FVEthreshold" = 0.95, "methodMuCovEst" = "smooth"))
  DPC1 <- FPCAder(fpcaObj1, derOptns = list(method='DPC', bwMu=0.05, bwCov=0.3))
  
  K1 <- min(ncol(FPCA_Dy1$phi), ncol(DPC1$phiDer))
  
  # Refit
  muRes <- muder_fun(yraw = t(fdata$data_noisy[[1]]@X), argvals = argvals)
  muDer <- muRes$muder_hat
  Refit_DPC1 <- tcrossprod(DPC1$phiDer[,1:K1], DPC1$xiDer[,1:K1]) + muDer
  
  # Evaluation
  ISE_list <- ISE_fun(argvals, FPCA_Dy1$phi[,1:K1], DPC1$phiDer[,1:K1])
  RE_list <- RE_fun(FPCA_Dy1$lambda[1:K1], DPC1$lambdaDer[1:K1])
  MSE_list <- MSE_fun(FPCA_Dy1$xiEst[,1:K1], DPC1$xiDer[,1:K1])
  RMISE_list <- RMISE_fun(argvals, t(fdata$data_der[[1]]@X), Refit_DPC1)
  
  list(RE_list = RE_list, ISE_list = ISE_list, MSE_list = MSE_list, RMISE_list = RMISE_list)
}
end_time <- proc.time()
time_case2_fpcader <- end_time - start_time
saveRDS(y1_case2_FPCAder, file = "./SimRes/y1_case2_fpcader.rds")
saveRDS(time_case2_fpcader, file = "./SimRes/time_case2_fpcader.rds")


# Sparse ----
## y1: SP + sparse ----
start_time <- proc.time()
y1_case3_sp <- foreach(rep=1:Nrep, .packages = c("fda", "fdapace", "funData", "mgcv")) %dopar% {
  # Data generation 
  set.seed(rep)
  fdata <- dgp(n_obs, argvals, noise_variance, min_obs, max_obs)
  
  # Ground truth
  Lt <- rep(list(argvals), n_obs)
  LDy1 <- split(fdata$data_der[[1]]@X, 1:n_obs)
  FPCA_Dy1 <- FPCA(Ly = LDy1, Lt = Lt, optns = list("FVEthreshold" = 0.95))
  
  # Single penalty
  y1_dfpcaRes <-  deriv_fpca_fun_sp(yraw = t(fdata$data_sparse[[1]]@X), argvals = argvals, 
                                    d1 = d1, K.p = K.p)
  
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
end_time <- proc.time()
time_case3_sp <- end_time - start_time
saveRDS(y1_case3_sp, file = "./SimRes/y1_case3_sp.rds")
saveRDS(time_case3_sp, file = "./SimRes/time_case3_sp.rds")


## y1: AP + sparse ----
start_time <- proc.time()
y1_case3_ap <- foreach(rep=1:Nrep, .packages = c("fda", "fdapace", "funData", "mgcv")) %dopar% {
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
  RMISE_list <- RMISE_fun(argvals, t(fdata$data_der[[1]]@X), y1_dfpcaRes$Refit)
  
  list(RE_list = RE_list, ISE_list = ISE_list, MSE_list = MSE_list, RMISE_list = RMISE_list)
}
end_time <- proc.time()
time_case3_ap <- end_time - start_time
saveRDS(y1_case3_ap, file = "./SimRes/y1_case3_ap.rds")
saveRDS(time_case3_ap, file = "./SimRes/time_case3_ap.rds")


## y1: FPCAder + sparse ----
start_time <- proc.time()
y1_case3_FPCAder <- foreach(rep=1:Nrep, .packages = c("fda", "fdapace", "funData", "mgcv")) %dopar% {
  # Data generation 
  set.seed(rep)
  fdata <- dgp(n_obs, argvals, noise_variance, min_obs, max_obs)
  
  # Ground truth
  Lt <- rep(list(argvals), n_obs)
  LDy1 <- split(fdata$data_der[[1]]@X, 1:n_obs)
  FPCA_Dy1 <- FPCA(Ly = LDy1, Lt = Lt, optns = list("FVEthreshold" = 0.95))
  
  # FPCAder
  Ly1 <- split(fdata$data_sparse[[1]]@X, 1:n_obs)
  fpcaObj1 <- FPCA(Ly = Ly1, Lt = Lt, 
                   optns = list("FVEthreshold" = 0.95, "methodMuCovEst" = "smooth", nRegGrid = length(argvals)))
  DPC1 <- FPCAder(fpcaObj1, derOptns = list(method='DPC'))
  
  K1 <- min(ncol(FPCA_Dy1$phi), ncol(DPC1$phiDer))
  
  # Refit
  muRes <- muder_fun(yraw = t(fdata$data_sparse[[1]]@X), argvals = argvals)
  muDer <- muRes$muder_hat
  Refit_DPC1 <- tcrossprod(DPC1$phiDer[,1:K1], DPC1$xiDer[,1:K1]) + muDer
  
  # Evaluation
  ISE_list <- ISE_fun(argvals, FPCA_Dy1$phi[,1:K1], DPC1$phiDer[,1:K1])
  RE_list <- RE_fun(FPCA_Dy1$lambda[1:K1], DPC1$lambdaDer[1:K1])
  MSE_list <- MSE_fun(FPCA_Dy1$xiEst[,1:K1], DPC1$xiDer[,1:K1])
  RMISE_list <- RMISE_fun(argvals, t(fdata$data_der[[1]]@X), Refit_DPC1)
  
  list(RE_list = RE_list, ISE_list = ISE_list, MSE_list = MSE_list, RMISE_list = RMISE_list)
}
end_time <- proc.time()
time_case3_fpcader <- end_time - start_time
saveRDS(y1_case3_FPCAder, file = "./SimRes/y1_case3_fpcader.rds")
saveRDS(time_case3_fpcader, file = "./SimRes/time_case3_fpcader.rds")



# No noise ----
## y2: SP + no noise ----
start_time <- proc.time()
y2_case1_sp <- foreach(rep=1:Nrep, .packages = c("fda", "fdapace", "funData", "mgcv")) %dopar% {
  # Data generation 
  set.seed(rep)
  fdata <- dgp(n_obs, argvals, noise_variance, min_obs, max_obs)
  
  # Ground truth
  Lt <- rep(list(argvals), n_obs)
  LDy2 <- split(fdata$data_der[[2]]@X, 1:n_obs)
  FPCA_Dy2 <- FPCA(Ly = LDy2, Lt = Lt, optns = list("FVEthreshold" = 0.95))
  
  # Single penalty
  y2_dfpcaRes <-  deriv_fpca_fun_sp(yraw = t(fdata$data[[2]]@X), 
                                    argvals = argvals, d1 = d1, K.p = K.p)
  # Truncation
  K2 <- min(ncol(FPCA_Dy2$phi), ncol(y2_dfpcaRes$Phi))
  
  # Refit
  Refit_DFPCA <- 
    tcrossprod(y2_dfpcaRes$Phi[,1:K2], y2_dfpcaRes$xi[,1:K2])+y2_dfpcaRes$mu
  
  # Evaluation
  ISE_list <- ISE_fun(argvals, FPCA_Dy2$phi[,1:K2], y2_dfpcaRes$Phi[,1:K2])
  RE_list <- RE_fun(FPCA_Dy2$lambda[1:K2], y2_dfpcaRes$lambda[1:K2])
  MSE_list <- MSE_fun(FPCA_Dy2$xiEst[,1:K2], y2_dfpcaRes$xi[,1:K2])
  RMISE_list <- RMISE_fun(argvals, t(fdata$data_der[[2]]@X), Refit_DFPCA)
  
  list(RE_list = RE_list, ISE_list = ISE_list, MSE_list = MSE_list, RMISE_list = RMISE_list)
}
end_time <- proc.time()
y2_time_case1_sp <- end_time - start_time
saveRDS(y2_case1_sp, file = "./SimRes/y2_case1_sp.rds")
saveRDS(y2_time_case1_sp, file = "./SimRes/y2_time_case1_sp.rds")


## y2: AP + no noise ----
start_time <- proc.time()
y2_case1_ap <- foreach(rep=1:Nrep, .packages = c("fda", "fdapace", "funData", "mgcv")) %dopar% {
  # Data generation 
  set.seed(rep)
  fdata <- dgp(n_obs, argvals, noise_variance, min_obs, max_obs)
  
  # Ground truth
  Lt <- rep(list(argvals), n_obs)
  LDy2 <- split(fdata$data_der[[2]]@X, 1:n_obs)
  FPCA_Dy2 <- FPCA(Ly = LDy2, Lt = Lt, optns = list("FVEthreshold" = 0.95))
  
  # Additive penalty
  y2_dfpcaRes <-  deriv_fpca_fun(yraw = t(fdata$data[[2]]@X),
                                 argvals = argvals, d1 = d1, d2 = d2, K.p = K.p)
  
  K2 <- min(ncol(FPCA_Dy2$phi), ncol(y2_dfpcaRes$Phi))
  
  # Refit
  Refit_DFPCA <- 
    tcrossprod(y2_dfpcaRes$Phi[,1:K2], y2_dfpcaRes$xi[,1:K2])+y2_dfpcaRes$mu
  
  # Evaluation
  ISE_list <- ISE_fun(argvals, FPCA_Dy2$phi[,1:K2], y2_dfpcaRes$Phi[,1:K2])
  RE_list <- RE_fun(FPCA_Dy2$lambda[1:K2], y2_dfpcaRes$lambda[1:K2])
  MSE_list <- MSE_fun(FPCA_Dy2$xiEst[,1:K2], y2_dfpcaRes$xi[,1:K2])
  RMISE_list <- RMISE_fun(argvals, t(fdata$data_der[[2]]@X), Refit_DFPCA)
  
  list(RE_list = RE_list, ISE_list = ISE_list, MSE_list = MSE_list, RMISE_list = RMISE_list)
}
end_time <- proc.time()
y2_time_case1_ap <- end_time - start_time
saveRDS(y2_case1_ap, file = "./SimRes/y2_case1_ap.rds")
saveRDS(y2_time_case1_ap, file = "./SimRes/y2_time_case1_ap.rds")


## y2: FPCAder + no noise ----
start_time <- proc.time()
y2_case1_FPCAder <- foreach(rep=1:Nrep, .packages = c("fda", "fdapace", "funData", "mgcv")) %dopar% {
  # Data generation 
  set.seed(rep)
  fdata <- dgp(n_obs, argvals, noise_variance, min_obs, max_obs)
  
  # Ground truth
  Lt <- rep(list(argvals), n_obs)
  LDy2 <- split(fdata$data_der[[2]]@X, 1:n_obs)
  FPCA_Dy2 <- FPCA(Ly = LDy2, Lt = Lt, optns = list("FVEthreshold" = 0.95))
  
  # FPCAder
  Ly2 <- split(fdata$data[[2]]@X, 1:n_obs)
  fpcaObj2 <- FPCA(Ly = Ly2, Lt = Lt, optns = list("FVEthreshold" = 0.95))
  DPC2 <- FPCAder(fpcaObj2, derOptns = list(method='DPC'))
  
  K2 <- min(ncol(FPCA_Dy2$phi), ncol(DPC2$phiDer))
  
  # Refit
  muRes <- muder_fun(yraw = t(fdata$data[[2]]@X), argvals = argvals)
  muDer <- muRes$muder_hat
  Refit_DPC2 <- tcrossprod(DPC2$phiDer[,1:K2], DPC2$xiDer[,1:K2]) + muDer
  
  # Evaluation
  ISE_list <- ISE_fun(argvals, FPCA_Dy2$phi[,1:K2], DPC2$phiDer[,1:K2])
  RE_list <- RE_fun(FPCA_Dy2$lambda[1:K2], DPC2$lambdaDer[1:K2])
  MSE_list <- MSE_fun(FPCA_Dy2$xiEst[,1:K2], DPC2$xiDer[,1:K2])
  RMISE_list <- RMISE_fun(argvals, t(fdata$data_der[[2]]@X), Refit_DPC2)
  
  list(RE_list = RE_list, ISE_list = ISE_list, MSE_list = MSE_list, RMISE_list = RMISE_list)
}
end_time <- proc.time()
y2_time_case1_fpcader <- end_time - start_time
saveRDS(y2_case1_FPCAder, file = "./SimRes/y2_case1_fpcader.rds")
saveRDS(y2_time_case1_fpcader, file = "./SimRes/y2_time_case1_fpcader.rds")


# Noisy ----
## y2: SP + noisy ----
start_time <- proc.time()
y2_case2_sp <- foreach(rep=1:Nrep, .packages = c("fda", "fdapace", "funData", "mgcv")) %dopar% {
  # Data generation 
  set.seed(rep)
  fdata <- dgp(n_obs, argvals, noise_variance, min_obs, max_obs)
  
  # Ground truth
  Lt <- rep(list(argvals), n_obs)
  LDy2 <- split(fdata$data_der[[2]]@X, 1:n_obs)
  FPCA_Dy2 <- FPCA(Ly = LDy2, Lt = Lt, optns = list("FVEthreshold" = 0.95))
  
  # Single penalty
  y2_dfpcaRes <-  deriv_fpca_fun_sp(yraw = t(fdata$data_noisy[[2]]@X), 
                                    argvals = argvals, d1 = d1, K.p = K.p)
  # Truncation
  K2 <- min(ncol(FPCA_Dy2$phi), ncol(y2_dfpcaRes$Phi))
  
  # Refit
  Refit_DFPCA <- 
    tcrossprod(y2_dfpcaRes$Phi[,1:K2], y2_dfpcaRes$xi[,1:K2])+y2_dfpcaRes$mu
  
  # Evaluation
  ISE_list <- ISE_fun(argvals, FPCA_Dy2$phi[,1:K2], y2_dfpcaRes$Phi[,1:K2])
  RE_list <- RE_fun(FPCA_Dy2$lambda[1:K2], y2_dfpcaRes$lambda[1:K2])
  MSE_list <- MSE_fun(FPCA_Dy2$xiEst[,1:K2], y2_dfpcaRes$xi[,1:K2])
  RMISE_list <- RMISE_fun(argvals, t(fdata$data_der[[2]]@X), Refit_DFPCA)
  
  list(RE_list = RE_list, ISE_list = ISE_list, MSE_list = MSE_list, RMISE_list = RMISE_list)
}
end_time <- proc.time()
y2_time_case2_sp <- end_time - start_time
saveRDS(y2_case2_sp, file = "./SimRes/y2_case2_sp.rds")
saveRDS(y2_time_case2_sp, file = "./SimRes/y2_time_case2_sp.rds")


## y2: AP + noisy ----
start_time <- proc.time()
y2_case2_ap <- foreach(rep=1:Nrep, .packages = c("fda", "fdapace", "funData", "mgcv")) %dopar% {
  # Data generation 
  set.seed(rep)
  fdata <- dgp(n_obs, argvals, noise_variance, min_obs, max_obs)
  
  # Ground truth
  Lt <- rep(list(argvals), n_obs)
  LDy2 <- split(fdata$data_der[[2]]@X, 1:n_obs)
  FPCA_Dy2 <- FPCA(Ly = LDy2, Lt = Lt, optns = list("FVEthreshold" = 0.95))
  
  # Additive penalty
  y2_dfpcaRes <-  deriv_fpca_fun(yraw = t(fdata$data_noisy[[2]]@X),
                                 argvals = argvals, d1 = d1, d2 = d2, K.p = K.p)
  
  K2 <- min(ncol(FPCA_Dy2$phi), ncol(y2_dfpcaRes$Phi))

  # Refit
  Refit_DFPCA <- 
    tcrossprod(y2_dfpcaRes$Phi[,1:K2], y2_dfpcaRes$xi[,1:K2])+y2_dfpcaRes$mu
  
  # Evaluation
  ISE_list <- ISE_fun(argvals, FPCA_Dy2$phi[,1:K2], y2_dfpcaRes$Phi[,1:K2])
  RE_list <- RE_fun(FPCA_Dy2$lambda[1:K2], y2_dfpcaRes$lambda[1:K2])
  MSE_list <- MSE_fun(FPCA_Dy2$xiEst[,1:K2], y2_dfpcaRes$xi[,1:K2])
  RMISE_list <- RMISE_fun(argvals, t(fdata$data_der[[2]]@X), Refit_DFPCA)
  
  list(RE_list = RE_list, ISE_list = ISE_list, MSE_list = MSE_list, RMISE_list = RMISE_list)
}
end_time <- proc.time()
y2_time_case2_ap <- end_time - start_time
saveRDS(y2_case2_ap, file = "./SimRes/y2_case2_ap.rds")
saveRDS(y2_time_case2_ap, file = "./SimRes/y2_time_case2_ap.rds")


## y2: FPCAder + noisy ----
start_time <- proc.time()
y2_case2_FPCAder <- foreach(rep=1:Nrep, .packages = c("fda", "fdapace", "funData", "mgcv")) %dopar% {
  # Data generation 
  set.seed(rep)
  fdata <- dgp(n_obs, argvals, noise_variance, min_obs, max_obs)
  
  # Ground truth
  Lt <- rep(list(argvals), n_obs)
  LDy2 <- split(fdata$data_der[[2]]@X, 1:n_obs)
  FPCA_Dy2 <- FPCA(Ly = LDy2, Lt = Lt, optns = list("FVEthreshold" = 0.95))
  
  # FPCAder
  Ly2 <- split(fdata$data_noisy[[2]]@X, 1:n_obs)
  fpcaObj2 <- FPCA(Ly = Ly2, Lt = Lt, 
                   optns = list("FVEthreshold" = 0.95, "methodMuCovEst" = "smooth"))
  DPC2 <- FPCAder(fpcaObj2, derOptns = list(method='DPC', bwMu=0.05, bwCov=0.3))
  K2 <- min(ncol(FPCA_Dy2$phi), ncol(DPC2$phiDer))
  
  # Refit
  muRes <- muder_fun(yraw = t(fdata$data_noisy[[2]]@X), argvals = argvals)
  muDer <- muRes$muder_hat
  Refit_DPC2 <- tcrossprod(DPC2$phiDer[,1:K2], DPC2$xiDer[,1:K2]) + muDer
  
  # Evaluation
  ISE_list <- ISE_fun(argvals, FPCA_Dy2$phi[,1:K2], DPC2$phiDer[,1:K2])
  RE_list <- RE_fun(FPCA_Dy2$lambda[1:K2], DPC2$lambdaDer[1:K2])
  MSE_list <- MSE_fun(FPCA_Dy2$xiEst[,1:K2], DPC2$xiDer[,1:K2])
  RMISE_list <- RMISE_fun(argvals, t(fdata$data_der[[2]]@X), Refit_DPC2)
  
  list(RE_list = RE_list, ISE_list = ISE_list, MSE_list = MSE_list, RMISE_list = RMISE_list)
}
end_time <- proc.time()
y2_time_case2_fpcader <- end_time - start_time
saveRDS(y2_case2_FPCAder, file = "./SimRes/y2_case2_fpcader.rds")
saveRDS(y2_time_case2_fpcader, file = "./SimRes/y2_time_case2_fpcader.rds")


# Sparse ----
## y2: SP + sparse ----
start_time <- proc.time()
y2_case3_sp <- foreach(rep=1:Nrep, .packages = c("fda", "fdapace", "funData", "mgcv")) %dopar% {
  # Data generation 
  set.seed(rep)
  fdata <- dgp(n_obs, argvals, noise_variance, min_obs, max_obs)
  
  # Ground truth
  Lt <- rep(list(argvals), n_obs)
  LDy2 <- split(fdata$data_der[[2]]@X, 1:n_obs)
  FPCA_Dy2 <- FPCA(Ly = LDy2, Lt = Lt, optns = list("FVEthreshold" = 0.95))
  
  # Single penalty
  y2_dfpcaRes <-  deriv_fpca_fun_sp(yraw = t(fdata$data_sparse[[2]]@X), 
                                    argvals = argvals, d1 = d1, K.p = K.p)
  # Truncation
  K2 <- min(ncol(FPCA_Dy2$phi), ncol(y2_dfpcaRes$Phi))
  
  # Refit
  Refit_DFPCA <- 
    tcrossprod(y2_dfpcaRes$Phi[,1:K2], y2_dfpcaRes$xi[,1:K2])+y2_dfpcaRes$mu
  
  # Evaluation
  ISE_list <- ISE_fun(argvals, FPCA_Dy2$phi[,1:K2], y2_dfpcaRes$Phi[,1:K2])
  RE_list <- RE_fun(FPCA_Dy2$lambda[1:K2], y2_dfpcaRes$lambda[1:K2])
  MSE_list <- MSE_fun(FPCA_Dy2$xiEst[,1:K2], y2_dfpcaRes$xi[,1:K2])
  RMISE_list <- RMISE_fun(argvals, t(fdata$data_der[[2]]@X), Refit_DFPCA)
  
  list(RE_list = RE_list, ISE_list = ISE_list, MSE_list = MSE_list, RMISE_list = RMISE_list)
}
end_time <- proc.time()
y2_time_case3_sp <- end_time - start_time
saveRDS(y2_case3_sp, file = "./SimRes/y2_case3_sp.rds")
saveRDS(y2_time_case3_sp, file = "./SimRes/y2_time_case3_sp.rds")


## y2: AP + sparse ----
start_time <- proc.time()
y2_case3_ap <- foreach(rep=1:Nrep, .packages = c("fda", "fdapace", "funData", "mgcv")) %dopar% {
  # Data generation 
  set.seed(rep)
  fdata <- dgp(n_obs, argvals, noise_variance, min_obs, max_obs)
  
  # Ground truth
  Lt <- rep(list(argvals), n_obs)
  LDy2 <- split(fdata$data_der[[2]]@X, 1:n_obs)
  FPCA_Dy2 <- FPCA(Ly = LDy2, Lt = Lt, optns = list("FVEthreshold" = 0.95))
  
  # Additive penalty
  y2_dfpcaRes <-  deriv_fpca_fun(yraw = t(fdata$data_sparse[[2]]@X),
                                 argvals = argvals, d1 = d1, d2 = d2, K.p = K.p)
  
  K2 <- min(ncol(FPCA_Dy2$phi), ncol(y2_dfpcaRes$Phi))
  
  # Refit
  Refit_DFPCA <- 
    tcrossprod(y2_dfpcaRes$Phi[,1:K2], y2_dfpcaRes$xi[,1:K2])+y2_dfpcaRes$mu
  
  # Evaluation
  ISE_list <- ISE_fun(argvals, FPCA_Dy2$phi[,1:K2], y2_dfpcaRes$Phi[,1:K2])
  RE_list <- RE_fun(FPCA_Dy2$lambda[1:K2], y2_dfpcaRes$lambda[1:K2])
  MSE_list <- MSE_fun(FPCA_Dy2$xiEst[,1:K2], y2_dfpcaRes$xi[,1:K2])
  RMISE_list <- RMISE_fun(argvals, t(fdata$data_der[[2]]@X), Refit_DFPCA)
  
  list(RE_list = RE_list, ISE_list = ISE_list, MSE_list = MSE_list, RMISE_list = RMISE_list)
}
end_time <- proc.time()
y2_time_case3_ap <- end_time - start_time
saveRDS(y2_case3_ap, file = "./SimRes/y2_case3_ap.rds")
saveRDS(y2_time_case3_ap, file = "./SimRes/y2_time_case3_ap.rds")


## y2: FPCAder + sparse ----
start_time <- proc.time()
y2_case3_FPCAder <- foreach(rep=1:Nrep, .packages = c("fda", "fdapace", "funData", "mgcv")) %dopar% {
  # Data generation 
  set.seed(rep)
  fdata <- dgp(n_obs, argvals, noise_variance, min_obs, max_obs)
  
  # Ground truth
  Lt <- rep(list(argvals), n_obs)
  LDy2 <- split(fdata$data_der[[2]]@X, 1:n_obs)
  FPCA_Dy2 <- FPCA(Ly = LDy2, Lt = Lt, optns = list("FVEthreshold" = 0.95))
  
  # FPCAder
  Ly2 <- split(fdata$data_sparse[[2]]@X, 1:n_obs)
  fpcaObj2 <- FPCA(Ly = Ly2, Lt = Lt, 
                   optns = list("FVEthreshold" = 0.95, "methodMuCovEst" = "smooth", nRegGrid = length(argvals)))
  DPC2 <- FPCAder(fpcaObj2, derOptns = list(method='DPC'))
  K2 <- min(ncol(FPCA_Dy2$phi), ncol(DPC2$phiDer))
  
  # Refit
  muRes <- muder_fun(yraw = t(fdata$data_sparse[[2]]@X), argvals = argvals)
  muDer <- muRes$muder_hat
  Refit_DPC2 <- tcrossprod(DPC2$phiDer[,1:K2], DPC2$xiDer[,1:K2]) + muDer
  
  # Evaluation
  ISE_list <- ISE_fun(argvals, FPCA_Dy2$phi[,1:K2], DPC2$phiDer[,1:K2])
  RE_list <- RE_fun(FPCA_Dy2$lambda[1:K2], DPC2$lambdaDer[1:K2])
  MSE_list <- MSE_fun(FPCA_Dy2$xiEst[,1:K2], DPC2$xiDer[,1:K2])
  RMISE_list <- RMISE_fun(argvals, t(fdata$data_der[[2]]@X), Refit_DPC2)
  
  list(RE_list = RE_list, ISE_list = ISE_list, MSE_list = MSE_list, RMISE_list = RMISE_list)
}
end_time <- proc.time()
y2_time_case3_fpcader <- end_time - start_time
saveRDS(y2_case3_FPCAder, file = "./SimRes/y2_case3_fpcader.rds")
saveRDS(y2_time_case3_fpcader, file = "./SimRes/y2_time_case3_fpcader.rds")
