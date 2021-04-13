# Title     : TODO
# Objective : TODO
# Created by: jonlachmann
# Created on: 2021-04-13

{
  library(RCurl)
  ### Download simulated logistic data as per example 2
  logistic_x <- read.csv(header=F, text=getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/supplementaries/Mode%20Jumping%20MCMC/supplementary/examples/Simulated%20Logistic%20Data%20With%20Multiple%20Modes%20(Example%203)/sim3-X.txt"))
  logistic_y <- read.csv(header=F, text=getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/supplementaries/Mode%20Jumping%20MCMC/supplementary/examples/Simulated%20Logistic%20Data%20With%20Multiple%20Modes%20(Example%203)/sim3-Y.txt"))
  colnames(logistic_y) <- "Y"
  # Modify data
  logistic_x$V2<-(logistic_x$V10+logistic_x$V14)*logistic_x$V9
  logistic_x$V5<-(logistic_x$V11+logistic_x$V15)*logistic_x$V12

  logistic_x <- as.matrix(cbind(1,logistic_x))
  logistic_y <- as.matrix(logistic_y)
}


mod_count <- 500

mliks <- matrix(NA, mod_count, 4)
for (i in 1:mod_count) {
  model <- as.logical(c(T,intToBits(i*100)[1:20]))
  log_x_mod <- logistic_x[,model]
  glm_mod <- glm.fit(log_x_mod, logistic_y, family=binomial())
  quant_mod <- irls(log_x_mod, logistic_y, binomial(), 0.5, 1, cooling = c(10,0.8,0.8), expl=c(100,3))
  sub_mod <- irls(log_x_mod, logistic_y, binomial(), 1, 0.1, cooling = c(4,0.7,0.7), expl=c(3,1.5))
  both_mod <- irls(log_x_mod, logistic_y, binomial(), 0.5, 0.1, cooling = c(4,0.7,0.7), expl=c(3,1.5))
  mliks[i,1] <- -glm_mod$deviance/2
  mliks[i,2] <- quant_mod$loglik
  mliks[i,3] <- sub_mod$loglik
  mliks[i,4] <- both_mod$loglik
  print(i)
}

multiplot(mliks[,c(1,2)], ylim=c(-2000,-1000))

sim.map