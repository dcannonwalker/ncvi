# Transpose of differential wrt mvn mean from poisson likelihood,
# where the second fixed effect parameter has a mixture prior of
# Gaussian and point mass at zero.
d_mixture_mean_poisson <- function(data, pars) {
  C1 <- data$C1
  C0 <- data$C0
  y <- data$y
  mu <- pars$phi$mu
  Sigma <- pars$phi$Sigma
  pi <- pars$pi
  # Wrap in 'c()' to drop matrix to vector
  A1 <- c(C1 %*% mu + 0.5 * diag(C1 %*% Sigma %*% t(C1)))
  A0 <- c(C0 %*% mu + 0.5 * diag(C0 %*% Sigma %*% t(C0)))
  (1 - pi) * t(C1) %*% (y - exp(A1)) + (pi) * t(C0) %*% (y - exp(A0))
}

# Vec inverse of differential wrt vec of mvn covariance from poisson likelihood,
# where the second fixed effect parameter has a mixture prior of
# Gaussian and point mass at zero.
d_mixture_cov_poisson <- function(data, pars) {
  C1 <- data$C1
  C0 <- data$C0
  mu <- pars$phi$mu
  Sigma <- pars$phi$Sigma
  pi <- pars$pi
  # Wrap in 'c()' to drop matrix to vector
  # Should use 'drop()' instead?
  A1 <- c(C1 %*% mu + 0.5 * diag(C1 %*% Sigma %*% t(C1)))
  A0 <- c(C0 %*% mu + 0.5 * diag(C0 %*% Sigma %*% t(C0)))

  (1 - pi) * t(C1) %*% diag(exp(A1)) %*% C1 +
    pi * t(C0) %*% diag(exp(A0)) %*% C0
}
