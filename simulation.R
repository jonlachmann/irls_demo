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
system.time(sub_mod2 <- irls(million_x, million_y_l, binomial(), 1, 0.01, maxit=300, cooling = c(4,0.9,0.99), expl=c(3,1.5)))
both_mod <- irls(million_x, million_y_l, binomial(), 0.5, 0.1, cooling = c(4,0.7,0.7), expl=c(3,1.5))

glm.fit()

min(sub_mod$devhist)
plot(sub_mod$devhist[10:300], type="l")
par(mfrow=c(4,4))
for (i in 1:16) {
  multiplot(cbind(sub_mod2$betahist[,i], true_betas[,i]))
}

library(GMJMCMC)

true_betas <- matrix(rep(glmmm$coefficients, 299), 299, byrow=T)

plot(quant_mod$devhist, type="l")
get_deviance(sub_mod$betahist[11,], million_x, million_y_l, binomial())
mean(sub_mod$devhist[100:200])

logistic.loglik.aic <- function (y, x, model, complex, params) {
  suppressWarnings({mod <- glm.fit(as.matrix(x[,model]), y, family=binomial())})
  ret <- -(mod$deviance/2) - mod$rank
  return(ret)
}

full_model_count <- 2^15

# Calculate the full model set
full_10K <- vector("list", full_model_count)
progress <- 0
for (i in 1:full_model_count) {
  modelvector <- as.logical(c(T,intToBits(i)[1:15]))
  loglik <- logistic.loglik.aic(million_y_l[1:10000], million_x[1:10000,], modelvector, NULL, NULL)
  full_10K[[i]] <- list(prob=NA, model=modelvector[-1], crit=loglik, alpha=NA)
  if (i %% floor(full_model_count/40) == 0) progress <- print.progressbar(progress, 40)
}

system.time(glmmm <- glm.fit(million_x[1:10000,], million_y_l[1:10000], family=binomial()))
system.time(sub_mod <- irls(million_x[1:10000,], million_y_l[1:10000], binomial(), 1, 0.1, maxit=300, cooling = c(3,0.9,0.95), expl=c(3,3)))
min(sub_mod$devhist)
get_deviance(sub_mod$betahist[which.min(sub_mod$devhist),], million_x[1:10000,], million_y_l[1:10000], binomial())

logistic.loglik.aic.sub <- function (y, x, model, complex, params) {
  suppressWarnings({mod <- irls(as.matrix(x[,model]), y, binomial(), 1, 0.1, maxit=300, cooling = c(3,0.9,0.95), expl=c(3,3))})
  ret1 <- mod$loglik - mod$rank
  dev2 <- get_deviance(mod$betahist[which.min(mod$devhist),], as.matrix(x[,model]), y, binomial())
  ret2 <- -(dev2/2) - mod$rank
  return(list(sub=ret1, full=ret2))
}

full_10K_sub <- vector("list", full_model_count)
full_10K_sub_full <- vector("list", full_model_count)
library(parallel)
full_10K_res <- mclapply(X=1:full_model_count, FUN=function(i) {
  modelvector <- as.logical(c(T,intToBits(i)[1:15]))
  try({
    logliks <- logistic.loglik.aic.sub(million_y_l[1:10000], million_x[1:10000,], modelvector, NULL, NULL)
    return(list(sub=list(prob=NA, model=modelvector[-1], crit=logliks$sub, alpha=NA),
              sub_full=list(prob=NA, model=modelvector[-1], crit=logliks$full, alpha=NA)))
  })
})

comparison <- matrix(NA,full_model_count, 3)
for (i in 1:full_model_count){
  try({comparison[i,1] <- full_10K[[i]]$crit})
  try({comparison[i,2] <- full_10K_res[[i]]$sub$crit})
  try({comparison[i,3] <- full_10K_res[[i]]$sub_full$crit})
  print(i)
}
full_10K[[1]]$crit
full_10K_res[[1]]$sub$crit

comparison <- comparison[!is.na(comparison[,2]),]

comparison_diffs <- matrix(NA, nrow(comparison), 2)
comparison_diffs[,1] <- comparison[,1] - comparison[,2]
comparison_diffs[,2] <- comparison[,1] - comparison[,3]

cor(comparison)

quantiles <- (comparison>=quantile(comparison[,1], (0.99)))

bestmods <- matrix(comparison[quantiles], ncol=3, byrow=T)

rmsediff1_m <- sqrt(mean((modelprobs[,2]-modelprobs[,1])^2))
rmsediff2_m <- sqrt(mean((modelprobs[,3]-modelprobs[,1])^2))

modelprobs <- bestmods / colSums(bestmods)
modelprobs_all <- comparison / colSums(comparison)

rmsediff1_m_all <- sqrt(mean((modelprobs_all[,2]-modelprobs_all[,1])^2))
rmsediff2_m_all <- sqrt(mean((modelprobs_all[,3]-modelprobs_all[,1])^2))


multiplot(matrix(comparison[quantiles], ncol=3, byrow=T))

matrix(comparison[quantiles], ncol=3, byrow=T)

rmsediff1 <- sqrt(mean((comparison_diffs[,1])^2))
rmsediff2 <- sqrt(mean((comparison_diffs[,2])^2))

multiplot(comparison_diffs[1:500,], ylim=c(-150,200))

multiplot(comparison[1:500,])


truee <- matrix(unlist(full_10K), ncol=18, byrow=T)
v1_mat <- matrix(unlist(sub_models_v1), ncol=18, byrow=T)

hist(truee[,17], breaks=50)
hist(v1_mat[,17], breaks=50)

for (i in 1:full_model_count) {
  if (abs(sub_models_v1[[i]]$crit -)
}

margprobs_10K_true <- marginal.probs.renorm(full_10K)
margprobs_10K_v1 <- marginal.probs.renorm(sub_models_v1)
margprobs_10K_v2 <- marginal.probs.renorm(sub_models_v2)

par(mfrow=c(1,1))
barplot(margprobs_10K_true, main="True marginal inclusion probabilities")
barplot(margprobs_10K_v1, main="Estimated marginal inclusion probabilities")

sub_models_v1 <- vector("list", full_model_count)
sub_models_v2 <- vector("list", full_model_count)
for (i in 1:full_model_count) {
  if (length(full_10K_res[[i]]) == 2) {
    sub_models_v1[[i]] <- full_10K_res[[i]][[1]]
    sub_models_v2[[i]] <- full_10K_res[[i]][[2]]
  }
}
for (i in 1:full_model_count) {
  if (length(sub_models_v1[[i]]) != 4) sub_models_v1[[i]] <- NULL
  if (length(sub_models_v2[[i]]) != 4) sub_models_v1[[i]] <- NULL
}



save(full_10K_res, file="full10Kres")

for (i in 1:10) {
  modelvector <- as.logical(c(T,intToBits(1)[1:15]))
  logliks <- logistic.loglik.aic.sub(million_y_l[1:10000], million_x[1:10000,], modelvector, NULL, NULL)
  full_10K_sub[[i]] <- list(prob=NA, model=modelvector[-1], crit=logliks$sub, alpha=NA)
  full_10K_sub_full[[i]] <- list(prob=NA, model=modelvector[-1], crit=logliks$full, alpha=NA)

  if (i %% floor(full_model_count/40) == 0) progress <- print.progressbar(progress, 40)
}

