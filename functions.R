# Title     : General functions
# Objective : Functions to be used when analysing the IRLS algorithm
# Created by: jonlachmann
# Created on: 2021-04-13

# Calculate the full deviance for a model
get_deviance <- function(beta, X, y, family) {
  mu <- family$linkinv(X %*% beta)
  sum(family$dev.resids(y, mu, rep(1,nrow(X))))
}

# Plot many columns in a matrix, log scale can be enabled too
multiplot <- function (mat, logscale=F, ylim=c(min(mat), max(mat))) {
  if (logscale) {
    mat[mat > 0] <- log(mat[mat > 0])
    mat[mat < 0] <- -log(-mat[mat < 0])
  }
  rbcol <- rainbow(ncol(mat))
  plot(mat[,1], type="l", ylim=ylim, col=rbcol[1])
  for (i in 2:ncol(mat)) lines(mat[,i], col=rbcol[i])
}

# Print a progress bar while iterating over a population
print.progressbar <- function (progress, size=40) {
  cat("\r", "|")
  for (p in 1:size-1) {
    if (progress >= p) cat("=")
    else cat(" ")
  }
  cat("|")
  return(progress+1)
}

# Function for calculating feature importance through renormalized model estimates
marginal.probs.renorm <- function (models) {
  model.size <- length(models[[1]]$model)
  models.matrix <- matrix(unlist(models), ncol=model.size+3, byrow=T)
  models.matrix <- models.matrix[(!duplicated(models.matrix[,2:(model.size+1)], dim=1)),]
  max_mlik <- max(models.matrix[,(model.size+2)])+57.71633
  crit.sum <- sum(exp(mpfr(models.matrix[,(model.size+2)]-max_mlik,128)))
  probs <- matrix(NA,1,model.size)
  for (i in 2:(model.size+1)) probs[i-1] <- asNumeric(sum(exp(mpfr(models.matrix[as.logical(models.matrix[,i]),(model.size+2)]-max_mlik,128)))/crit.sum)
  return(probs)
}