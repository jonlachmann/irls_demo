# Title     : TODO
# Objective : TODO
# Created by: jonlachmann
# Created on: 2021-04-13

nvars <- 15
nobs <- 10^6

library(mvtnorm)
# Generate data
{
  set.seed(1911)
  covmat <- matrix(rnorm(nvars^2, runif(nvars^2, -5,5), runif(nvars^2, 0, 5)), nvars)
  covmat <- covmat %*% t(covmat)

  million_x <- cbind(1, matrix(rmvnorm(nobs, runif(nvars, -5,5), covmat), nobs))
  covars <- sample.int(nvars, 8)
  betas <- runif(9, -10, 10)
  million_y_g <- million[,c(1,covars)] %*% betas + rnorm(nobs, 0, 3)
  million_y_l <- rbinom(nobs, 1, (1/(1+exp(-million_y_g))))
}

system.time(glmmm <- glm.fit(million_x[1:1000000,], million_y_l[1:1000000], family=binomial()))
system.time(glmmm <- glm.fit(million_x[1:1000000,], million_y_l[1:1000000], family=binomial(), start=sub_mod$betahist[121,]))

system.time(quant_mod <- irls(million_x, million_y_l, binomial(), 0.1, 1, cooling = c(10,0.8,0.9), expl=c(100,3)))
system.time(sub_mod <- irls(million_x, million_y_l, binomial(), 1, 0.003, maxit=200, cooling = c(4,0.9,0.95), expl=c(3,3)))
both_mod <- irls(million_x, million_y_l, binomial(), 0.5, 0.1, cooling = c(4,0.7,0.7), expl=c(3,1.5))

glm.fit()

min(sub_mod$devhist)
plot(sub_mod$devhist, type="l")
multiplot(sub_mod$betahist)
plot(quant_mod$devhist, type="l")
get_deviance(sub_mod$betahist[121,], million_x, million_y_l, binomial())

