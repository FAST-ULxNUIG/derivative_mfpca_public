################################################################################
# Multivariate FPCAder
# Extend FPCAder to Multivariate Cases following Happ's method
################################################################################

FPCAder_multi <- function(fdata, methodMuCovEst = "cross-sectional", bwMu=0.05, bwCov=0.1){
  # methodMuCovEst: smooth for noisy and sparse cases

  n_obs <- funData::nObs(fdata)
  n_features <- length(fdata)
  
  # Computation of the univariate derivative-based results
  eigfun_uni_list <- list()
  scores_uni_list <- list()
  mu_uni_list <- list()
  
  for (p in 1:n_features) {
    argvals <- unlist(fdata[[p]]@argvals)
    muRes <- muder_fun(yraw = t(fdata[[p]]@X), argvals = argvals)
    mu_uni_list[[p]] <- muRes$muder_hat
    
    Lt <- rep(list(argvals), n_obs)
    Ly <- split(fdata[[p]]@X, 1:n_obs)
    fpcaObj <- FPCA(Ly = Ly, Lt = Lt, 
                    optns = list("FVEthreshold" = 0.95, "methodMuCovEst" = methodMuCovEst, 
                                 nRegGrid = length(argvals)))
    DPC <- FPCAder(fpcaObj, derOptns = list(method='DPC', bwMu=bwMu, bwCov=bwCov))
    eigfun_uni_list[[p]] <- DPC$phiDer
    scores_uni_list[[p]] <- DPC$xiDer
  }
  
  # Computation of the multivariate derivative-based results
  npc <- sapply(scores_uni_list, function(x) ncol(x))
  scores <- do.call(cbind, scores_uni_list) # univariate scores
  scores_cov <- cov(scores)
  eigen_obj <- eigen(scores_cov, symmetric = T)
  
  # Computation of the multivariate DFPCs
  eigenfunctions <- list()
  for (p in 1:n_features) {
    if (p == 1) {
      idx <- 1:npc[p]
    } else{
      idx <- 1:npc[p] + sum(npc[1:(p - 1)])
    }
    eigenfunctions[[p]] <- eigfun_uni_list[[p]] %*% eigen_obj$vectors[idx,]
  }
  
  # Computation of the multivariate derivative scores
  scores <- scores %*% eigen_obj$vectors
  
  # Variation explained (>= 95%)
  positiveInd <- eigen_obj[["values"]] >= 0
  eigval <- sort(eigen_obj[["values"]][positiveInd], decreasing = TRUE)
  FVE <- cumsum(eigval)/sum(eigval)
  K <- min(which(FVE > 0.95)) # the number of eigenfuns
  
  # Truncation
  eigenvalues <- eigen_obj$values[1:K]
  eigenfunctions <- lapply(
    1:n_features,
    function(p) funData::funData(argvals, t(eigenfunctions[[p]][, 1:K]))
  )
  scores <- scores[, 1:K]
  
  # Refit
  Refit <- list()
  for (p in 1:n_features) {
    argvals <- fdata[[p]]@argvals[[1]]
    cent <- scores %*% eigenfunctions[[p]]@X
    values <- sweep(cent, 2, mu_uni_list[[p]], "+")
    Refit[[p]] <- funData(argvals, values)
  }
  
  return(list(
    mu = mu_uni_list,
    eigenvalues = eigenvalues,
    eigenfunctions = funData::multiFunData(eigenfunctions),
    scores = scores,
    FVE = FVE[1:K],
    Refit = funData::multiFunData(Refit)
  ))
}
