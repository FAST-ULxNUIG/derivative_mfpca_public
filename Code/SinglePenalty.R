################################################################################
# Single penalty 
################################################################################

# Pre-setting for basic matrices ----
pre_setting_sp <- function(argvals, d1, K.p, l=1, p=3){
  # x: a vector containing observed locations on the functional domain
  # d1: difference order for P1 
  # d2: difference order for P2
  # l: derivative order; default 1
  # K.p: the number of knots;
  # p: degree; cubic splines; default 3
  knots <- seq(-p, K.p + p, length = K.p + 1 + 2 * p) / K.p
  knots <- knots * (max(argvals) - min(argvals)) + min(argvals)
  
  B <- spline.des(knots = knots, x = argvals, ord = p + 1, 
                  outer.ok = TRUE, sparse = TRUE)$design
  Bder <- spline.des(knots = knots, x = argvals, ord = p + 1, derivs = l, 
                     outer.ok = TRUE, sparse = TRUE)$design
  G <- crossprod(B) #G=t(B)%*%B
  eigG <- eigen(G, symmetric = TRUE)
  G_invhalf <- eigG$vectors %*% diag(1/sqrt(eigG$values)) %*% t(eigG$vectors)
  nbasis <- ncol(B)
  D1 <- diff(diag(1, nrow=nbasis), diff=d1)
  P1 <- crossprod(D1)
  
  H <- G_invhalf %*% P1 %*% G_invhalf
  H <- (H + t(H)) / 2
  eigH <- eigen(H, symmetric = TRUE)
  
  U <- eigH$vectors
  s  <- pmax(eigH$values, 0)
  if(min(s)<=0.0000001) {
    s <- s + 0.000001;
  }
  
  A0 <- G_invhalf%*%U
  As <- B %*% G_invhalf %*% U
  return(list(B=as.matrix(B), Bder=as.matrix(Bder), G_invhalf=G_invhalf, 
              P1=P1, U=U, s=s, A0=A0, As=as.matrix(As)))
}




# gcv ----
gcv_fun_sp <- function(Y, As, s, totalmiss,
                       lower, upper, search.length){
  # Y: observations of dimension nxL
  # As: from pre_setting_sp
  # s: from pre_setting_sp
  L <- ncol(Y)
  
  Ytilde <- crossprod(As, t(Y))
  C_diag <- rowSums(Ytilde^2, na.rm = TRUE)
  Y_square <- sum(Y^2, na.rm = TRUE)
  Ytilde_square <- sum(Ytilde^2, na.rm = TRUE)
  
  gcv_temp <- function(x) {
    lambda <- exp(x)
    lambda_s <- (lambda * s)^2/(1 + lambda * s)^2
    gcv <- sum(C_diag * lambda_s) - Ytilde_square + Y_square
    trace <- sum(1/(1 + lambda * s))
    gcv <- gcv/(1 - trace/L/(1 - totalmiss))^2
    return(gcv)
  }
  
  Lambda <- seq(lower, upper, length = search.length)
  Length <- length(Lambda)
  Gcv <- rep(0, Length)
  for (i in 1:Length) Gcv[i] <- gcv_temp(Lambda[i])
  i0 <- which.min(Gcv)
  lambda <- exp(Lambda[i0])
  
  res <- list("lambda"=lambda, "gcv"=Gcv[i0])
  return(res)
}


# Theta_fun ----
Theta_fun_sp <- function(Y, argvals, d1, K.p, totalmiss, 
                         lower = -20, upper = 20, 
                         search.length = 100){
  # Y: demeaned observations with dimension nxL
  # argvals: a vector containing observed locations on the functional domain
  # d1/d2: difference order 
  # K.p: the number of knots;
  # totalmiss: 0 for dense data
  # lower, upper: select lambda based on GCV
  L <- ncol(Y)
  n <- nrow(Y)
  
  List <- pre_setting_sp(argvals = argvals, d1 = d1, K.p = K.p)
  B <- List$B
  P1 <- List$P1
  G_invhalf <- List$G_invhalf
  A0 <- List$A0
  As <- List$As
  s <- List$s
  
  
  # GCV
  fit <- gcv_fun_sp(Y = Y, As = As, s = s, 
                    totalmiss = totalmiss, lower = lower, 
                    upper = upper, search.length = search.length)
  
  rho <- fit$lambda
  gcv <- fit$gcv
  
  Sigmas <- diag(1/(1+rho*s))
  Ytilde <- crossprod(As, t(Y)) #Ytilde=t(AS)%*%Y
  YS <- Sigmas%*%Ytilde
  
  # temp0 = YS%*%t(YS)/n = A%*%Sigma%*%A^{T}
  temp0 <- YS %*% t(YS)/n
  # Theta (of dimension cxc)
  Theta <- A0%*%temp0%*%t(A0)
  
  # when data are noisy
  Kraw <- crossprod(Y)/n # Kraw
  Ksmooth <- B %*% Theta %*% t(B)
  gap <- diag(Kraw)-diag(Ksmooth)
  sigma2 <- fdapace:::trapzRcpp(argvals, pmax(gap, 0))
  if (isTRUE(sigma2 != 0)) {
    temp1 <- diag(sigma2/(1+rho*s)^2)
    Theta <- A0%*%(temp0-temp1)%*%t(A0)
  }
  return(list(Theta=Theta, A0=A0, Sigmas=Sigmas, sigma2=sigma2,
              rho=rho, gcv=gcv))
}


# Estimation of eigfuns and eigvals ----
# Weights function
quadWeights <- function(argvals, method = "trapezoidal"){
  ret <- switch(method,
                trapezoidal = {D <- length(argvals)
                1/2*c(argvals[2] - argvals[1], argvals[3:D] -argvals[1:(D-2)], argvals[D] - argvals[D-1])},
                midpoint = c(0,diff(argvals)),
                stop("function quadWeights: choose either trapezoidal or midpoint quadrature rule"))
  return(ret)
}

# function to estimate derivative-based eigfun and eigval
dfpca_eigen_fun_sp <- function(argvals, Bder, Theta){
  # argvals: a vector containing observed locations on the functional domain
  # Bder: Lxc matrix from pre_setting function
  # Theta: cxc matrix from Theta_fun function
  # weights
  w <- quadWeights(argvals = argvals)
  W <- diag(w)
  # Gint
  Gint <- t(Bder) %*% W %*% Bder
  eigGint <- eigen(Gint, symmetric = TRUE)
  Gint_half <- eigGint$vectors %*% diag(sqrt(eigGint$values)) %*% t(eigGint$vectors)
  Gint_invhalf <- eigGint$vectors %*% diag(sqrt(1/eigGint$values)) %*% t(eigGint$vectors)
  # Eigenfuns and eigenvals
  M <- Gint_half %*% Theta %*% Gint_half
  eigM <- eigen(M)
  eigvals <- eigM$values[eigM$values>0]
  per <- cumsum(eigvals)/sum(eigvals)
  U <- eigM$vectors
  Phi <- Bder %*% Gint_invhalf %*% U
  K <- min(which(per > 0.95)) # the number of eigenfuns
  FVE <- cumsum(eigvals[1:K])/sum(eigvals)
  return(list(eigvals = eigvals[1:K], eigfuns=Phi[,1:K], FVE=FVE))
}


# Derivative-based FPCA results ----
deriv_fpca_fun_sp <- function(yraw, argvals, d1, K.p){
  # yraw: Lxn observations
  # argvals: a vector containing observed locations on the functional domain
  # d1/d2: difference penalty order
  # K.p: the number of knots;
  
  # Pre-setting
  n <- ncol(yraw)
  L <- nrow(yraw)
  List <- pre_setting_sp(argvals = argvals, d1 = d1, K.p = K.p)
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
  
  # Mean centred
  x_all <- rep(argvals, times = n)
  y_all <- as.vector(yraw)
  
  mu_hat <- c(JOPS::psNormal(x = x_all, y = y_all, xl = min(argvals), xr = max(argvals),
                             xgrid = argvals)$ygrid)
  muder_hat <- JOPS::psNormal_Deriv(x = x_all, y = y_all, xl = min(argvals), xr = max(argvals),
                                    xgrid = argvals)$d_pred
  Ycentered <- sweep(yraw, 1, mu_hat, "-")
  
  # Theta
  ThetaRes <- Theta_fun_sp(Y = t(Ycentered), argvals = argvals, 
                           d1 = d1, K.p = K.p, totalmiss = totalmiss)
  Theta_mat <- ThetaRes$Theta
  A0 <- ThetaRes$A0
  Sigmas <- ThetaRes$Sigmas
  
  # Estimate derivative-based eigfun and eigval
  EigRes <- dfpca_eigen_fun_sp(argvals = argvals, Bder = Bder, Theta = Theta_mat)
  Phi <- EigRes$eigfuns
  lambda <- EigRes$eigvals
  # Estimate derivative-based scores
  # Smoothed curves
  A0Sigmas <- A0 %*% Sigmas %*% t(A0)
  Yhat <- B %*% A0Sigmas %*% t(B) %*% Ycentered
  Yder_hat <- Bder %*% A0Sigmas %*% t(B) %*% Ycentered
  
  # 'raw' derivative covariance 
  Kder_raw <- Yder_hat %*% t(Yder_hat) /n 
  # 'smoothed' derivative covariance
  Kder_smooth <- Bder %*% ThetaRes$Theta %*% t(Bder)
  gap <- diag(Kder_raw) - diag(Kder_smooth)
  sigma2 <- fdapace:::trapzRcpp(argvals, pmax(gap, 0))
  
  # Derivative-based scores
  if (length(lambda) == 1) {
    #note that Phi and lambda are derivatived-based
    Phi <- matrix(Phi, ncol=1)
    tmp <- crossprod(Phi) + sigma2/lambda
  } else{
    tmp <- crossprod(Phi) + sigma2*diag(1/lambda)
  }
  xi <- solve(tmp) %*% t(Phi) %*% Yder_hat
  
  # Derivative-based refitted curves
  Refit <- Phi %*% xi + muder_hat
  
  return(list(mu = muder_hat, Phi = Phi, lambda = lambda, sigma2=sigma2,
              FVE = EigRes$FVE, xi = t(xi), Refit = Refit,
              rho = ThetaRes$rho, gcv = ThetaRes$gcv))
}


# Multivariate derivative FPCA ----
deriv_mfpca_fun_sp <- function(fdata, d1, K.p){
  # fdata: multivariate functional data
  # K.p: the number of knots;
  
  n_obs <- funData::nObs(fdata)
  n_features <- length(fdata)
  
  # Computation of the univariate derivative-based results
  eigfun_uni_list <- list()
  scores_uni_list <- list()
  mu_uni_list <- list()
  for (p in 1:n_features) {
    argvals <- fdata[[p]]@argvals[[1]]
    yraw <- t(fdata[[p]]@X)
    dfpcaRes <- 
      deriv_fpca_fun_sp(yraw = yraw, argvals = argvals, 
                        d1 = d1, K.p = K.p)
    
    eigfun_uni_list[[p]] <- dfpcaRes$Phi
    scores_uni_list[[p]] <- dfpcaRes$xi
    mu_uni_list[[p]] <- dfpcaRes$mu
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

