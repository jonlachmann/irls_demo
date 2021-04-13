# Title     : General functions
# Objective : Functions to be used when analysing the IRLS algorithm
# Created by: jonlachmann
# Created on: 2021-04-13

# Calculate the full deviance for a model
get_deviance <- function(beta, X, y, family) {
  mu <- family$linkinv(X %*% beta)
  sum(family_b$dev.resids(y, mu, rep(1,nrow(X))))
}

# Plot many columns in a matrix, log scale can be enabled too
multiplot <- function (mat, logscale=F) {
  if (logscale) {
    mat[mat > 0] <- log(mat[mat > 0])
    mat[mat < 0] <- -log(-mat[mat < 0])
  }
  plot(mat[,1], type="l", ylim=c(min(mat),max(mat)))
  for (i in 2:ncol(mat)) lines(mat[,i], col=colors()[i*5])
}