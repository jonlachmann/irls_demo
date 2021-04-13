# Title     : IRLS with subsampling
# Objective : IRLS algorithm with WLS quantile optimization and stochastic subsampling
# Created by: jonlachmann
# Created on: 2021-04-13

# X is the covariates matrix
# y is the dependent variable
# family is the glm family
# quant is the top quantile of observations to use for WLS, (0-1]
# subs is the proportion of observations to use for subsampling
# maxit is the maximum number of iterations
# tol is the convergence tolerance
# cooling is the cooling schedule parameters, i.e. iterations to stay constant for, initial value, decay value
# expl is the settings to avoid exploding deviance in subsampling, i.e. iterations to not check, proportional increase signifying explosion

irls <- function (X, y, family, quant=1, subs=1, maxit=100, tol=1e-7, cooling = c(3,0.9,0.9), expl = c(3,1.5)) {
  temp <- cooling[2]
  nobs <- nrow(X)
  nvars <- ncol(X)
  weights <- rep(1,nobs)
  eval(family$initialize)

  # Get initial eta, mu and deviance
  eta <- family$linkfun(mustart)
  if (subs != 1) {
    # Get subsample
    sub_size <- nobs*subs
    subsi <- sample.int(nobs, sub_size, replace=F)

    eta <- eta[subsi]
    mu <- family$linkinv(eta)
    dev <- sum(family$dev.resids(y[subsi], mu, weights[subsi]))/subs
  } else {
    mu <- family$linkinv(eta)
    dev <- sum(family$dev.resids(y, mu, weights))
  }

  devhist <- matrix(NA,maxit,1)
  betahist <- matrix(NA,maxit,nvars)
  devhist[1] <- dev
  iter <- 1
  conv <- F
  while (!conv & iter < maxit) {
    # Get var_mu
    var_mu <- family$variance(mu)
    # Get mu eta
    mu_eta <- family$mu.eta(eta)
    # Get z
    if (subs != 1) z <- eta + (y[subsi] - mu) / mu_eta
    else z <- eta + (y - mu) / mu_eta
    # Get w
    w <- sqrt((mu_eta^2)/var_mu)

    if (iter > 1) betaold <- beta
    # Do WLS, use quantile if set
    if (quant != 1) {
      quanti <- (w>=quantile(w, (1-quant)))
      if (subs != 1) X_s <- X[subsi,][quanti,]
      else X_s <- X[quanti,]
      w_s <- w[quanti]
      z_s <- z[quanti]
      beta <- .Call(stats:::C_Cdqrls, X_s * as.numeric(w_s), z_s * w_s, 1e-7, check=FALSE)$coefficients
    } else if (subs != 1) {
      beta <- .Call(stats:::C_Cdqrls, X[subsi,] * as.numeric(w), z * w, 1e-7, check=FALSE)$coefficients
    } else {
      beta <- .Call(stats:::C_Cdqrls, X * as.numeric(w), z * w, 1e-7, check=FALSE)$coefficients
    }

    # Do new subsampling
    if (subs != 1) subsi <- sample.int(nobs, sub_size, replace=F)

    # Get eta
    if (subs != 1) eta <- X[subsi,] %*% beta
    else eta <- X %*% beta

    # Get mu
    mu <- family$linkinv(eta)

    # Calculate deviance
    devold <- dev
    if (subs != 1) dev <- sum(family$dev.resids(y[subsi], mu, weights[subsi]))/subs
    else dev <- sum(family$dev.resids(y, mu, weights))

    # Do cooling schedule
    if (iter > cooling[1]+2) {
      beta <- temp * beta + (1-temp) * betaold
      temp <- temp * cooling[3]
    }

    # Check convergence
    if ((abs(dev - devold)) / (0.1 + abs(dev)) < tol) {
      conv <- T
    }
    # Check for exploding deviance in subsampling
    if (iter > expl[1] && (dev-devold) / abs(devold) > expl[2]) {
      print(paste0("Exploding deviance at iteration ", iter))
      break
    }

    betahist[iter,] <- beta
    iter <- iter + 1
    devhist[iter] <- dev
  }
  return(list(coefficients=beta, iter=iter, loglik=-devhist[iter]/2, devhist=devhist[1:iter], betahist=betahist[1:(iter-1),]))
}