################################################################################
# Define the multivariate functional data generating process
################################################################################

init_dgp <- function() {
  # Mean vector
  mvn_mean <- c(0, 0.5, 3.75)
  
  # Covariance matrix
  sds <- diag(c(1, 0.14, 0.7))
  rho <- 0.2
  rho_mat <- matrix(rho, nrow = 3, ncol = 3)
  diag(rho_mat) <- 1
  mvn_var <- sds %*% rho_mat %*% sds
  
  equations <- list(
    expression(a + 1 / (c * t / 5 + 2 * b * exp(-16 * t ^ 2))), 
    expression(a - cos((t * c/4) * (2 * t - pi)) + 2 * exp(-16 * b * t^2))
  )
  derivatives <- lapply(equations, D, name = 't')
  return(list(
    'mvn_mean' = mvn_mean,
    'mvn_var' = mvn_var,
    'equations' = equations,
    'derivatives' = derivatives
  ))
}

dgp <- function(n_obs, argvals, noise_variance, min_obs, max_obs) {
  # n_obs: number of observations (int)
  # argvals: sampling points (vector)
  # noise_variance: variance of the noise (float)
  # min_obs: minimum number of points (scalar)
  # max_obs: maximum number of points (scalar)
  init <- init_dgp()
  equations <- init$equations
  derivatives <- init$derivatives
  
  coefs <- MASS::mvrnorm(
    n = n_obs, mu = init$mvn_mean, Sigma = sqrt(init$mvn_var)
  )
  while (any(
    coefs[, 2] < 0.1 | coefs[, 2] > 1 | coefs[, 3] < 1 | coefs[, 3] > 6.5
  )) {
    cond1 <- coefs[, 2] < 0.1 | coefs[, 2] > 1 
    cond2 <- coefs[, 3] < 1 | coefs[, 3] > 6.5
    idx <- which(cond1 | cond2)
    coefs[idx, ] <- MASS::mvrnorm(
      length(idx), mu = init$mvn_mean, Sigma = sqrt(init$mvn_var)
    )
  }
  
  y1 <- y2 <- matrix(0, nrow = length(argvals), ncol = n_obs)
  Dy1 <- Dy2 <- matrix(0, nrow = length(argvals), ncol = n_obs)
  for (i in 1:n_obs) {
    t <- argvals
    a <- coefs[i, 1]
    b <- coefs[i, 2] 
    c <- coefs[i, 3] 
    y1[, i] <- eval(equations[[1]], list(a = a, b = b, c = c, t = t))
    y2[, i] <- eval(equations[[2]], list(a = a, b = b, c = c, t = t))
    
    Dy1[, i] <- eval(derivatives[[1]], list(a = a, b = b, c = c, t = t))
    Dy2[, i] <- eval(derivatives[[2]], list(a = a, b = b, c = c, t = t))
  }
  
  data <- funData::multiFunData(
    funData::funData(t, t(y1)),
    funData::funData(t, t(y2))
  )
  data_der <- funData::multiFunData(
    funData::funData(t, t(Dy1)),
    funData::funData(t, t(Dy2))
  )
  # Add noise to functional data
  data_noisy <- funData::addError(data, sd = noise_variance)
  
  # Sparsity
  minObs <- rep(min_obs, length(equations))
  maxObs <- rep(max_obs, length(equations))
  data_sparse <- funData::sparsify(data_noisy, minObs, maxObs)
  
  return(list(
    'data' = data,
    'data_der' = data_der,
    'data_noisy' = data_noisy,
    'data_sparse' = data_sparse
  ))
}