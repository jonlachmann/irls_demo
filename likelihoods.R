# Title     : TODO
# Objective : TODO
# Created by: jonlachmann
# Created on: 2021-04-14

logistic.loglik.aic.sub <- function (y, x, model, complex, params) {
  mod <- irls(as.matrix(x[,model]), y, binomial(), 1, params$subs, maxit=75, cooling = c(3,0.9,0.95), expl=c(3,1.5))
  ret1 <- mod$loglik - mod$rank
  dev2 <- get_deviance(mod$betahist[which.min(mod$devhist),], as.matrix(x[,model]), y, binomial())
  ret2 <- -(dev2/2) - mod$rank
  return(list(sub=ret1, full=ret2))
}

logistic.loglik.aic <- function (y, x, model, complex, params) {
  suppressWarnings({mod <- glm.fit(as.matrix(x[,model]), y, family=binomial())})
  ret <- -(mod$deviance/2) - mod$rank
  return(ret)
}