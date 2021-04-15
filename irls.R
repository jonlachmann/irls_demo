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

  # Get initial eta, mu and deviance
  if (subs != 1) {
    # Get subsample
    sub_size <- nobs*subs
    subsi <- sample.int(nobs, sub_size, replace=F)

    eta <- family$linkfun((y[subsi]+0.5)/2)
    mu <- family$linkinv(eta)
    dev <- sum(family$dev.resids(y[subsi], family$linkinv(rowSums(X[subsi,,drop=F])), 1))/subs
  } else {
    eta <- family$linkfun((y+0.5)/2)
    mu <- family$linkinv(eta)
    dev <- sum(family$dev.resids(y, family$linkinv(rowSums(X)), 1))
  }

  devhist <- matrix(NA,maxit,1)
  rankhist <- matrix(NA,maxit,1)
  betahist <- matrix(NA,maxit,nvars)
  devhist[1] <- dev
  betahist[1,] <- 0
  explosions <- 0
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
      if (subs != 1) X_s <- X[subsi,,drop=F][quanti,,drop=F]
      else X_s <- X[quanti,,drop=F]
      w_s <- w[quanti]
      z_s <- z[quanti]
      fit <- .Call(stats:::C_Cdqrls, X_s * as.numeric(w_s), z_s * w_s, 1e-7, check=FALSE)
      beta <- fit$coefficients
    } else if (subs != 1) {
      fit <- .Call(stats:::C_Cdqrls, X[subsi,,drop=F] * as.numeric(w), z * w, 1e-7, check=FALSE)
      beta <- fit$coefficients
    } else {
      fit <- .Call(stats:::C_Cdqrls, X * as.numeric(w), z * w, 1e-7, check=FALSE)
      beta <- fit$coefficients
    }
    rankhist[iter] <- fit$rank

    # Do new subsampling
    if (subs != 1) subsi <- sample.int(nobs, sub_size, replace=F)

    # Get eta
    if (subs != 1) eta <- X[subsi,,drop=F] %*% beta
    else eta <- X %*% beta

    # Get mu
    mu <- family$linkinv(eta)

    # Calculate deviance
    devold <- dev
    if (subs != 1) dev <- sum(family$dev.resids(y[subsi], mu, 1))/subs
    else dev <- sum(family$dev.resids(y, mu, 1))

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
      # Reset beta and start counting explosions
      beta <- betaold
      explosions <- explosions + 1
      if (explosions > 5) {
        print(paste0("5 exploding deviances at iteration ", iter))
        break
      }
    } else explosions <- 0

    iter <- iter + 1
    betahist[iter,] <- beta
    devhist[iter] <- dev
  }
  return(list(coefficients=beta, iter=iter, loglik=-devhist[iter]/2, rank=rankhist[iter-1], devhist=devhist[1:iter], betahist=betahist[1:iter, ,drop=F], rankhist=rankhist))
}