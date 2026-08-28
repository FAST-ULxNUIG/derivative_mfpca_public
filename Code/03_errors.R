################################################################################
# Metrics to evaluate the error between the ground truth and the estimation 
# (univariate + multivariate functional data)
################################################################################

# Univaraite cases  ----
# ISE for eigenfuns
ISE_fun <- function(argvals, ref, new){
  # ref/new: L*K eigenfunctions
  if (is.vector(new)) {
    new_fd <- flipFuns(funData(argvals, t(ref)), funData(argvals, t(new)))
    new <- t(new_fd@X)
    ise <- trapzRcpp(argvals, (ref - new) ^ 2)
  } else {
    K <- ncol(new)
    new_fd <- flipFuns(funData(argvals, t(ref[, 1:K])), funData(argvals, t(new[, 1:K])))
    new <- t(new_fd@X)
    ise <- c()
    for (k in 1:K) {
      error <- (ref[, k] - new[, k]) ^ 2
      ise[k] <- trapzRcpp(argvals, error)
    }
  }
  return(ise)
}


# RE for the eigenvalues 
RE_fun <- function(eigenvalues, estimates){
  # eigenvalues: true eigenvalues
  # estimates: estimation of the eigenvalues
  n_components <- length(eigenvalues)
  re <- sapply(1:n_components, function(k) {
    abs(eigenvalues[k] - estimates[k]) / eigenvalues[k]
  })
  return(re)
}

# MSE for the scores
# functions to flip scores
score_flipFun <- function(true_score, est_score){
  if (is.vector(est_score)) {
    signs <- sign(sum(est_score * true_score))
    flip_score <- signs*est_score
  } else {
    signs <- sign(colSums(est_score * true_score))
    flip_score <- sweep(est_score, 2, signs, `*`)
  }
  return(flip_score)
}


# functions to calculate MSE for scores
MSE_fun <- function(scores, estimates){
  if (is.vector(estimates)) {
    est_score <- score_flipFun(scores, estimates)
    MSE_res <- mean((scores-est_score)^2)/var(scores) 
  } else {
  K <- ncol(scores)
  est_score <- score_flipFun(scores, estimates)
  MSE_res <- sapply(1:K, function(k) mean((scores[,k]-est_score[,k])^2)/var(scores[,k]))
  }
  return(MSE_res)
}

# RMISE for the derivatives
RMISE_fun <- function(argvals, curves, estimates){
  # curves: true derivatives Lxn
  # estimates: estimates of the derivatives Lxn
  n_obs <- ncol(curves)
  error_norm <- sapply(1:n_obs, function(i) {
    tmp <- (curves[,i]-estimates[,i])^2
    fdapace::trapzRcpp(argvals, tmp)
  })
  actual_norm <- sapply(1:n_obs, function(i) {
    fdapace::trapzRcpp(argvals, curves[,i]^2)
  })
  rmise <- mean(error_norm/actual_norm)
  return(rmise)
}


# Multivariate cases ----
# RMISE for the derivatives
multi_RMISE_fun <- function(curves, estimates){
  # curves: true derivatives (class multiFunData)
  # estimates: estimates of the derivatives (class multiFunData)
  # Note: They should be sampled on the same grid.
  n_features <- length(curves)
  n_obs <- nObs(curves)
  
  rmise <- matrix(0, nrow = 1, ncol = n_features)
  for (p in 1:n_features) {
    arg <- curves[[p]]@argvals[[1]]
    error_norm <- mean(sapply(1:n_obs, function(i) {
      tmp <- (curves[[p]]@X[i, ] - estimates[[p]]@X[i, ])^2
      fdapace::trapzRcpp(arg, tmp)
    }))
    actual_norm <- mean(sapply(1:n_obs, function(i) {
      fdapace::trapzRcpp(arg, curves[[p]]@X[i, ]^2)
    }))
    rmise[1, p] <- error_norm / actual_norm
  }
  return(rmise)
}


# ISE for the eigenfunctions
multi_ISE_fun <- function(curves, estimates){
  # curves: true eigenfunctions (class multiFunData)
  # estimates: estimates of the eigenfunctions (class multiFunData)
  # Note: They should be sampled on the same grid.
  n_features <- length(curves)
  n_components <- dim(curves[[n_features]]@X)[1]
  
  # The eigenfunctions are defined up to a sign.
  estimates <- flipFuns(curves, estimates)
  
  ise <- rep(0, n_components)
  for (k in 1:n_components) {
    ise[k] <- sum(sapply(1:n_features, function(p) {
      arg <- curves[[p]]@argvals[[1]]
      tmp <- (curves[[p]]@X[k, ] - estimates[[p]]@X[k, ])^2
      fdapace::trapzRcpp(arg, tmp)
    }))
  }
  return(ise)
}


# RE for the eigenvalues
multi_RE_fun <- function(eigenvalues, estimates){
  # eigenvalues: true eigenvalues
  # estimates: estimation of the eigenvalues
  n_components <- length(eigenvalues)
  re <- sapply(1:n_components, function(k) {
    abs(eigenvalues[k] - estimates[k]) / eigenvalues[k]
  })
  return(re)
}

# MSE for the scores
# functions to flip scores
multi_score_flipFun <- function(true_score, est_score){
  signs <- sign(colSums(est_score * true_score))
  flip_score <- sweep(est_score, 2, signs, `*`)
  return(flip_score)
}
# functions to calculate MSE for scores
multi_MSE_fun <- function(scores, estimates){
  K <- ncol(scores)
  est_score <- multi_score_flipFun(scores, estimates)
  MSE_res <- sapply(1:K, function(k) mean((scores[,k]-est_score[,k])^2)/var(scores[,k]))
  return(MSE_res)
}
