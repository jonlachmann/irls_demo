# Title     : TODO
# Objective : TODO
# Created by: jonlachmann
# Created on: 2021-04-15

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

system.time(glmmm <- glm.fit(million_x[1:10000,], million_y_l[1:10000], family=binomial()))
system.time(sub_mod <- irls(million_x[1:10000,], million_y_l[1:10000], binomial(), 1, 0.1, maxit=300, cooling = c(3,0.9,0.95), expl=c(3,3)))
min(sub_mod$devhist)
get_deviance(sub_mod$betahist[which.min(sub_mod$devhist),], million_x[1:10000,], million_y_l[1:10000], binomial())


comparison <- matrix(NA,full_model_count, 3)
for (i in 1:full_model_count){
  try({comparison[i,1] <- full_10K[[i]]$crit})
  try({comparison[i,2] <- full_10K_res[[i]]$sub$crit})
  try({comparison[i,3] <- full_10K_res[[i]]$sub_full$crit})
  print(i)
}



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

modddd <- as.logical(c(T,intToBits(3456)[1:15]))

system.time(glmmodell <- glm.fit(as.matrix(million_x[1:10000,]), million_y_l[1:10000], family=binomial()))

system.time(irlsmo <- irls(million_x[1:10000,], million_y_l[1:10000], binomial(), 1, 0.1))

multiplot(irlsmo$betahist[1:15,])
abline(h=glmmodell$coefficients)

plot(irlsmo$devhist[50:100])
abline(h=glmmodell$deviance)
