################################################################################
# Multivariate Derivative-Based FPCA Simulation
# Single penalty vs Additive penalty vs FPCAder
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
source("Multivariate FPCAder.R")

# Multivariate: Select the first K eigenfunctions
eigFun_trucate <- function(mfdata, K){
  # mfdata: multivarite eigfun
  # K: truncation number
  n_obs <- funData::nObs(mfdata)
  n_features <- length(mfdata)
  temp <- list()
  for (p in 1:n_features) {
    argvals <- unlist(mfdata[[p]]@argvals)
    eigfun <- mfdata[[p]]@X[1:K,]
    temp[[p]] <- funData(argvals = argvals, X = eigfun)
  }
  res <- do.call(multiFunData, temp)
  return(res)
}

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
uniExpansions = list(list(type = "uFPCA"), list(type = "uFPCA"))


# Simulation ----
Nrep <- 500
cl <- parallel::makeCluster(6)
doParallel::registerDoParallel(cl)
RE_list <- ISE_list <- RMISE_list <- MSE_list <- list()

# No noise ----
## SP ----
start_time <- proc.time()
SP_nonoise <- foreach(rep=1:Nrep, 
                      .packages = c("fda", "fdapace", "funData", "mgcv", "MFPCA")) %dopar% {
      # Data generation 
      set.seed(rep)
      fdata_all <- dgp(n_obs, argvals, noise_variance, min_obs, max_obs)
      n_features <- length(fdata_all$data)
      # Ground truth
      mfpca_truth <- MFPCA(fdata_all$data_der, M=2, uniExpansions = uniExpansions)
      
      # SP
      dmfpcaRes <- deriv_mfpca_fun_sp(fdata = fdata_all$data, d1 = d1, K.p = K.p)
      eigfun_new <- eigFun_trucate(dmfpcaRes$eigenfunctions, 2)    
      
      ISE_list <- multi_ISE_fun(mfpca_truth$functions, eigfun_new)
      RE_list <- multi_RE_fun(mfpca_truth$values, dmfpcaRes$eigenvalues[1:2])
      MSE_list <- multi_MSE_fun(mfpca_truth$scores, dmfpcaRes$scores[,1:2])
      
      Refit <- list()
      for (p in 1:n_features) {
        temp <- dmfpcaRes$scores[, 1:2] %*% eigfun_new[[p]]@X
        temp <- sweep(temp, 2, dmfpcaRes$mu[[p]], "+")
        Refit[[p]] <- funData(argvals, temp)
      }
      RMISE_list <- c(multi_RMISE_fun(fdata_all$data_der, Refit))
      
      list(RE_list = RE_list, ISE_list = ISE_list, 
           MSE_list = MSE_list, RMISE_list = RMISE_list)
}
end_time <- proc.time()
time_case1_sp <- end_time - start_time
saveRDS(SP_nonoise, file = "./SimRes Multi/SP_nonoise.rds")
saveRDS(time_case1_sp, file = "./SimRes Multi/time_case1_sp.rds")

## AP ----
start_time <- proc.time()
AP_nonoise <- foreach(rep=1:Nrep, 
                      .packages = c("fda", "fdapace", "funData", "mgcv", "MFPCA")) %dopar% {
    # Data generation 
    set.seed(rep)
    fdata_all <- dgp(n_obs, argvals, noise_variance, min_obs, max_obs)
    n_features <- length(fdata_all$data)
    # Ground truth
    mfpca_truth <- MFPCA(fdata_all$data_der, M=2, uniExpansions = uniExpansions)
    
    # AP
    dmfpcaRes <- deriv_mfpca_fun(fdata = fdata_all$data, d1 = d1, d2 = d2, K.p = K.p)
    eigfun_new <- eigFun_trucate(dmfpcaRes$eigenfunctions, 2)    
    
    ISE_list <- multi_ISE_fun(mfpca_truth$functions, eigfun_new)
    RE_list <- multi_RE_fun(mfpca_truth$values, dmfpcaRes$eigenvalues[1:2])
    MSE_list <- multi_MSE_fun(mfpca_truth$scores, dmfpcaRes$scores[,1:2])
    
    Refit <- list()
    for (p in 1:n_features) {
      temp <- dmfpcaRes$scores[, 1:2] %*% eigfun_new[[p]]@X
      temp <- sweep(temp, 2, dmfpcaRes$mu[[p]], "+")
      Refit[[p]] <- funData(argvals, temp)
    }
    RMISE_list <- c(multi_RMISE_fun(fdata_all$data_der, Refit))
    
    list(RE_list = RE_list, ISE_list = ISE_list, 
         MSE_list = MSE_list, RMISE_list = RMISE_list)
}
end_time <- proc.time()
time_case1_ap <- end_time - start_time
saveRDS(AP_nonoise, file = "./SimRes Multi/AP_nonoise.rds")
saveRDS(time_case1_ap, file = "./SimRes Multi/time_case1_ap.rds")


## FPCAder ----
start_time <- proc.time()
fpcader_nonoise <- foreach(rep=1:Nrep, 
                      .packages = c("fda", "fdapace", "funData", "mgcv", "MFPCA")) %dopar% {
    # Data generation 
    set.seed(rep)
    fdata_all <- dgp(n_obs, argvals, noise_variance, min_obs, max_obs)
    n_features <- length(fdata_all$data)
    # Ground truth
    mfpca_truth <- MFPCA(fdata_all$data_der, M=2, uniExpansions = uniExpansions)
    
    # FPCAder
    FPCAderRes <- FPCAder_multi(fdata = fdata_all$data, methodMuCovEst = "cross-sectional")
    eigfun_new <- eigFun_trucate(FPCAderRes$eigenfunctions, 2)    
    
    ISE_list <- multi_ISE_fun(mfpca_truth$functions, eigfun_new)
    RE_list <- multi_RE_fun(mfpca_truth$values, FPCAderRes$eigenvalues[1:2])
    MSE_list <- multi_MSE_fun(mfpca_truth$scores, FPCAderRes$scores[,1:2])
    
    Refit <- list()
    for (p in 1:n_features) {
      temp <- FPCAderRes$scores[, 1:2] %*% eigfun_new[[p]]@X
      temp <- sweep(temp, 2, FPCAderRes$mu[[p]], "+")
      Refit[[p]] <- funData(argvals, temp)
    }
    RMISE_list <- c(multi_RMISE_fun(fdata_all$data_der, Refit))
    
    list(RE_list = RE_list, ISE_list = ISE_list, 
         MSE_list = MSE_list, RMISE_list = RMISE_list)
  }
end_time <- proc.time()
time_case1_fpcader <- end_time - start_time
saveRDS(fpcader_nonoise, file = "./SimRes Multi/fpcader_nonoise.rds")
saveRDS(time_case1_fpcader, file = "./SimRes Multi/time_case1_fpcader.rds")


# Noisy ----
## SP ----
start_time <- proc.time()
SP_noisy <- foreach(rep=1:Nrep, 
                      .packages = c("fda", "fdapace", "funData", "mgcv", "MFPCA")) %dopar% {
    # Data generation 
    set.seed(rep)
    fdata_all <- dgp(n_obs, argvals, noise_variance, min_obs, max_obs)
    n_features <- length(fdata_all$data)
    # Ground truth
    mfpca_truth <- MFPCA(fdata_all$data_der, M=2, uniExpansions = uniExpansions)
    
    # SP
    dmfpcaRes <- deriv_mfpca_fun_sp(fdata = fdata_all$data_noisy, d1 = d1, K.p = K.p)
    eigfun_new <- eigFun_trucate(dmfpcaRes$eigenfunctions, 2)    
    
    ISE_list <- multi_ISE_fun(mfpca_truth$functions, eigfun_new)
    RE_list <- multi_RE_fun(mfpca_truth$values, dmfpcaRes$eigenvalues[1:2])
    MSE_list <- multi_MSE_fun(mfpca_truth$scores, dmfpcaRes$scores[,1:2])
    
    Refit <- list()
    for (p in 1:n_features) {
      temp <- dmfpcaRes$scores[, 1:2] %*% eigfun_new[[p]]@X
      temp <- sweep(temp, 2, dmfpcaRes$mu[[p]], "+")
      Refit[[p]] <- funData(argvals, temp)
    }
    RMISE_list <- c(multi_RMISE_fun(fdata_all$data_der, Refit))
    
    list(RE_list = RE_list, ISE_list = ISE_list, 
         MSE_list = MSE_list, RMISE_list = RMISE_list)
}
end_time <- proc.time()
time_case2_sp <- end_time - start_time
saveRDS(SP_noisy, file = "./SimRes Multi/SP_noisy.rds")
saveRDS(time_case2_sp, file = "./SimRes Multi/time_case2_sp.rds")

## AP ----
start_time <- proc.time()
AP_noisy <- foreach(rep=1:Nrep, 
                      .packages = c("fda", "fdapace", "funData", "mgcv", "MFPCA")) %dopar% {
    # Data generation 
    set.seed(rep)
    fdata_all <- dgp(n_obs, argvals, noise_variance, min_obs, max_obs)
    n_features <- length(fdata_all$data)
    # Ground truth
    mfpca_truth <- MFPCA(fdata_all$data_der, M=2, uniExpansions = uniExpansions)
    
    # AP
    dmfpcaRes <- deriv_mfpca_fun(fdata = fdata_all$data_noisy, d1 = d1, d2 = d2, K.p = K.p)
    eigfun_new <- eigFun_trucate(dmfpcaRes$eigenfunctions, 2)    
    
    ISE_list <- multi_ISE_fun(mfpca_truth$functions, eigfun_new)
    RE_list <- multi_RE_fun(mfpca_truth$values, dmfpcaRes$eigenvalues[1:2])
    MSE_list <- multi_MSE_fun(mfpca_truth$scores, dmfpcaRes$scores[,1:2])
    
    Refit <- list()
    for (p in 1:n_features) {
      temp <- dmfpcaRes$scores[, 1:2] %*% eigfun_new[[p]]@X
      temp <- sweep(temp, 2, dmfpcaRes$mu[[p]], "+")
      Refit[[p]] <- funData(argvals, temp)
    }
    RMISE_list <- c(multi_RMISE_fun(fdata_all$data_der, Refit))
    
    list(RE_list = RE_list, ISE_list = ISE_list, 
         MSE_list = MSE_list, RMISE_list = RMISE_list)
  }
end_time <- proc.time()
time_case2_ap <- end_time - start_time
saveRDS(AP_noisy, file = "./SimRes Multi/AP_noisy.rds")
saveRDS(time_case2_ap, file = "./SimRes Multi/time_case2_ap.rds")

## FPCAder ----
start_time <- proc.time()
fpcader_noisy <- foreach(rep=1:Nrep, 
                           .packages = c("fda", "fdapace", "funData", "mgcv", "MFPCA")) %dopar% {
   # Data generation 
   set.seed(rep)
   fdata_all <- dgp(n_obs, argvals, noise_variance, min_obs, max_obs)
   n_features <- length(fdata_all$data)
   # Ground truth
   mfpca_truth <- MFPCA(fdata_all$data_der, M=2, uniExpansions = uniExpansions)
   
   # FPCAder
   FPCAderRes <- FPCAder_multi(fdata = fdata_all$data_noisy, methodMuCovEst = "smooth", bwMu=0.05, bwCov=0.3)
   eigfun_new <- eigFun_trucate(FPCAderRes$eigenfunctions, 2)    
   
   ISE_list <- multi_ISE_fun(mfpca_truth$functions, eigfun_new)
   RE_list <- multi_RE_fun(mfpca_truth$values, FPCAderRes$eigenvalues[1:2])
   MSE_list <- multi_MSE_fun(mfpca_truth$scores, FPCAderRes$scores[,1:2])
   
   Refit <- list()
   for (p in 1:n_features) {
     temp <- FPCAderRes$scores[, 1:2] %*% eigfun_new[[p]]@X
     temp <- sweep(temp, 2, FPCAderRes$mu[[p]], "+")
     Refit[[p]] <- funData(argvals, temp)
   }
   RMISE_list <- c(multi_RMISE_fun(fdata_all$data_der, Refit))
   
   list(RE_list = RE_list, ISE_list = ISE_list, 
        MSE_list = MSE_list, RMISE_list = RMISE_list)
}
end_time <- proc.time()
time_case2_fpcader <- end_time - start_time
saveRDS(fpcader_noisy, file = "./SimRes Multi/fpcader_noisy.rds")
saveRDS(time_case2_fpcader, file = "./SimRes Multi/time_case2_fpcader.rds")


# Sparse ----
## SP ----
start_time <- proc.time()
SP_sparse <- foreach(rep=1:Nrep, 
                    .packages = c("fda", "fdapace", "funData", "mgcv", "MFPCA")) %dopar% {
    # Data generation 
    set.seed(rep)
    fdata_all <- dgp(n_obs, argvals, noise_variance, min_obs, max_obs)
    n_features <- length(fdata_all$data)
    # Ground truth
    mfpca_truth <- MFPCA(fdata_all$data_der, M=2, uniExpansions = uniExpansions)
    
    # SP
    dmfpcaRes <- deriv_mfpca_fun_sp(fdata = fdata_all$data_sparse, d1 = d1, K.p = K.p)
    eigfun_new <- eigFun_trucate(dmfpcaRes$eigenfunctions, 2)    
    
    ISE_list <- multi_ISE_fun(mfpca_truth$functions, eigfun_new)
    RE_list <- multi_RE_fun(mfpca_truth$values, dmfpcaRes$eigenvalues[1:2])
    MSE_list <- multi_MSE_fun(mfpca_truth$scores, dmfpcaRes$scores[,1:2])
    
    Refit <- list()
    for (p in 1:n_features) {
      temp <- dmfpcaRes$scores[, 1:2] %*% eigfun_new[[p]]@X
      temp <- sweep(temp, 2, dmfpcaRes$mu[[p]], "+")
      Refit[[p]] <- funData(argvals, temp)
    }
    RMISE_list <- c(multi_RMISE_fun(fdata_all$data_der, Refit))
    
    list(RE_list = RE_list, ISE_list = ISE_list, 
         MSE_list = MSE_list, RMISE_list = RMISE_list)
  }
end_time <- proc.time()
time_case3_sp <- end_time - start_time
saveRDS(SP_sparse, file = "./SimRes Multi/SP_sparse.rds")
saveRDS(time_case3_sp, file = "./SimRes Multi/time_case3_sp.rds")

## AP ----
start_time <- proc.time()
AP_sparse <- foreach(rep=1:Nrep, 
                    .packages = c("fda", "fdapace", "funData", "mgcv", "MFPCA")) %dopar% {
    # Data generation 
    set.seed(rep)
    fdata_all <- dgp(n_obs, argvals, noise_variance, min_obs, max_obs)
    n_features <- length(fdata_all$data)
    # Ground truth
    mfpca_truth <- MFPCA(fdata_all$data_der, M=2, uniExpansions = uniExpansions)
    
    # AP
    dmfpcaRes <- deriv_mfpca_fun(fdata = fdata_all$data_sparse, d1 = d1, d2 = d2, K.p = K.p)
    eigfun_new <- eigFun_trucate(dmfpcaRes$eigenfunctions, 2)    
    
    ISE_list <- multi_ISE_fun(mfpca_truth$functions, eigfun_new)
    RE_list <- multi_RE_fun(mfpca_truth$values, dmfpcaRes$eigenvalues[1:2])
    MSE_list <- multi_MSE_fun(mfpca_truth$scores, dmfpcaRes$scores[,1:2])
    
    Refit <- list()
    for (p in 1:n_features) {
      temp <- dmfpcaRes$scores[, 1:2] %*% eigfun_new[[p]]@X
      temp <- sweep(temp, 2, dmfpcaRes$mu[[p]], "+")
      Refit[[p]] <- funData(argvals, temp)
    }
    RMISE_list <- c(multi_RMISE_fun(fdata_all$data_der, Refit))
    
    list(RE_list = RE_list, ISE_list = ISE_list, 
         MSE_list = MSE_list, RMISE_list = RMISE_list)
  }
end_time <- proc.time()
time_case3_ap <- end_time - start_time
saveRDS(AP_sparse, file = "./SimRes Multi/AP_sparse.rds")
saveRDS(time_case3_ap, file = "./SimRes Multi/time_case3_ap.rds")


## FPCAder ----
start_time <- proc.time()
fpcader_sparse <- foreach(rep=1:Nrep, 
                         .packages = c("fda", "fdapace", "funData", "mgcv", "MFPCA")) %dopar% {
   # Data generation 
   set.seed(rep)
   fdata_all <- dgp(n_obs, argvals, noise_variance, min_obs, max_obs)
   n_features <- length(fdata_all$data)
   # Ground truth
   mfpca_truth <- MFPCA(fdata_all$data_der, M=2, uniExpansions = uniExpansions)
   
   # FPCAder
   FPCAderRes <- FPCAder_multi(fdata = fdata_all$data_sparse, methodMuCovEst = "smooth")
   eigfun_new <- eigFun_trucate(FPCAderRes$eigenfunctions, 2)    
   
   ISE_list <- multi_ISE_fun(mfpca_truth$functions, eigfun_new)
   RE_list <- multi_RE_fun(mfpca_truth$values, FPCAderRes$eigenvalues[1:2])
   MSE_list <- multi_MSE_fun(mfpca_truth$scores, FPCAderRes$scores[,1:2])
   
   Refit <- list()
   for (p in 1:n_features) {
     temp <- FPCAderRes$scores[, 1:2] %*% eigfun_new[[p]]@X
     temp <- sweep(temp, 2, FPCAderRes$mu[[p]], "+")
     Refit[[p]] <- funData(argvals, temp)
   }
   RMISE_list <- c(multi_RMISE_fun(fdata_all$data_der, Refit))
   
   list(RE_list = RE_list, ISE_list = ISE_list, 
        MSE_list = MSE_list, RMISE_list = RMISE_list)
 }
end_time <- proc.time()
time_case3_fpcader <- end_time - start_time
saveRDS(fpcader_sparse, file = "./SimRes Multi/fpcader_sparse.rds")
saveRDS(time_case3_fpcader, file = "./SimRes Multi/time_case3_fpcader.rds")


