# Title     : Simulation study
# Objective : Explore the performance of subsampling IRLS with a simulation study
# Created by: jonlachmann
# Created on: 2021-04-13

# Set up simulations study parameters
nvars <- 15
nobs <- 10^6
full_model_count <- 2^nvars

library(mvtnorm)
# Generate the data to use
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

# Calculate the full model set using regular glm (SLOW!)
full_10K <- vector("list", full_model_count)
progress <- 0
for (i in 1:full_model_count) {
  modelvector <- as.logical(c(T,intToBits(i)[1:15]))
  loglik <- logistic.loglik.aic(million_y_l[1:10000], million_x[1:10000,], modelvector, NULL, NULL)
  full_10K[[i]] <- list(prob=NA, model=modelvector[-1], crit=loglik, alpha=NA)
  if (i %% floor(full_model_count/40) == 0) progress <- print.progressbar(progress, 40)
}

# Calculate the full model set using 1% at each iteration and two methods
full_10K_sub_1 <- vector("list", full_model_count)
full_10K_sub2_1 <- vector("list", full_model_count)
progress <- 0
for (i in 1:full_model_count) {
  modelvector <- as.logical(c(T,intToBits(i)[1:15]))
  logliks <- logistic.loglik.aic.sub(million_y_l[1:10000], million_x[1:10000,], modelvector, NULL, list(subs = 0.01))
  full_10K_sub_1[[i]] <- list(prob=NA, model=modelvector[-1], crit=logliks$sub, alpha=NA)
  full_10K_sub2_1[[i]] <- list(prob=NA, model=modelvector[-1], crit=logliks$full, alpha=NA)
  if (i %% floor(full_model_count/40) == 0) progress <- print.progressbar(progress, 40)
}

# Calculate the full model set using 2% at each iteration and two methods
full_10K_sub_2 <- vector("list", full_model_count)
full_10K_sub2_2 <- vector("list", full_model_count)
progress <- 0
for (i in 1:full_model_count) {
  modelvector <- as.logical(c(T,intToBits(i)[1:15]))
  logliks <- logistic.loglik.aic.sub(million_y_l[1:10000], million_x[1:10000,], modelvector, NULL, list(subs = 0.02))
  full_10K_sub_2[[i]] <- list(prob=NA, model=modelvector[-1], crit=logliks$sub, alpha=NA)
  full_10K_sub2_2[[i]] <- list(prob=NA, model=modelvector[-1], crit=logliks$full, alpha=NA)
  if (i %% floor(full_model_count/40) == 0) progress <- print.progressbar(progress, 40)
}

# Calculate the full model set using 10% at each iteration and two methods
full_10K_sub_10 <- vector("list", full_model_count)
full_10K_sub2_10 <- vector("list", full_model_count)
progress <- 0
for (i in 1:full_model_count) {
  modelvector <- as.logical(c(T,intToBits(i)[1:15]))
  logliks <- logistic.loglik.aic.sub(million_y_l[1:10000], million_x[1:10000,], modelvector, NULL, list(subs = 0.1))
  full_10K_sub_10[[i]] <- list(prob=NA, model=modelvector[-1], crit=logliks$sub, alpha=NA)
  full_10K_sub2_10[[i]] <- list(prob=NA, model=modelvector[-1], crit=logliks$full, alpha=NA)
  if (i %% floor(full_model_count/40) == 0) progress <- print.progressbar(progress, 40)
}

margprobs_10K_true <- marginal.probs.renorm(full_10K)
margprobs_10K_v1 <- marginal.probs.renorm(full_10K_sub_2)
margprobs_10K_v2 <- marginal.probs.renorm(full_10K_sub2_2)

barplot(margprobs_10K_true)
barplot(margprobs_10K_v1)
barplot(margprobs_10K_v2)



v1_2_mat <- matrix(unlist(full_10K_sub_2), ncol=18, byrow=T)
v2_2_mat <- matrix(unlist(full_10K_sub2_2), ncol=18, byrow=T)

hist(v1_2_mat[,17], breaks=150)
hist(v2_2_mat[,17], breaks=150)

compare2 <- matrix(NA, full_model_count, 3)
for (i in 1:full_model_count){
  try({compare2[i,1] <- full_10K[[i]]$crit})
  print(i)
}

multiplot(compare2[(compare2[,1] > -3000 & compare2[,1] < -2950),c(1,3)])

best50 <- compare2[(compare2[,1] > -114.54), c(1,3)]

best50exp <- exp(best50)

best50expmarg1 <- best50[,1] / sum(best50[,1])
best50expmarg2 <- best50[,2] / sum(best50[,2])
multiplot(cbind(best50expmarg1, best50expmarg2))

(compare2[,1] < -3000 & compare2[,1] > -2700)

compare2[,2] <- v1_2_mat[,17]
compare2[,3] <- v2_2_mat[,17]
cor(compare2)

ddiff <- compare2[,1] - compare2[,3]

compare2

hist(ddiff, breaks=1500, xlim=c(0,200), freq=F)

plot(ddiff[1:20000], type="l")

sort(compare2, )

mean(compare2[,1]) - mean(compare2[,3])

comparesort <- compare2[order(compare2[,1]),]

multiplot(comparesort[32200:32768,c(1,3)])

library(Rmpfr)


compare2mpfr1 <- mpfr(compare2[,1], 256)
compare2mpfr3 <- mpfr(compare2[,3], 256)

exps <- exp(compare2mpfr3)
sum(exps)

exps[32767]

exps <- cbind(exps, exp(compare2mpfr3-max(compare2mpfr3)+100))



