# Transpose of differential wrt mvn mean from poisson likelihood,
# where the second fixed effect parameter has a mixture prior of
# Gaussian and point mass at zero.
d_mvn_mean_mixture <- function(data, pars) {

  y <- data$y
  P <- data$P
  C  <- data$C
  if(is.null(data$S)) S <- 0
  else S <- data$S

  C0 <- C
  C0[, P] <- 0
  mu <- pars$phi$mu
  Sigma <- pars$phi$Sigma
  pi <- pars$pi
  # Wrap in 'c()' to drop matrix to vector
  A1 <- c(C %*% mu + S + 0.5 * diag(C %*% Sigma %*% t(C)))
  A0 <- c(C0 %*% mu + S + 0.5 * diag(C0 %*% Sigma %*% t(C0)))

  if (sum(is.infinite(A1)) > 0 | sum(is.infinite(A0)) > 0) {
    message(paste0("A1 or A0 is infite: ", A0, A1,
                   " and pi is: ", pi), "\n")
  }

  d_likelihood <-
    (1 - pi) * t(C) %*% (y - exp(A1)) +
    (pi) * t(C0) %*% (y - exp(A0))

  Tau <- pars$theta$Tau
  M <- pars$theta$M
  Tau %*% (M - mu) + d_likelihood
}

# Vec inverse of differential wrt vec of mvn covariance from poisson likelihood,
# where the second fixed effect parameter has a mixture prior of
# Gaussian and point mass at zero.
d_mvn_cov_mixture <- function(data, pars) {
  y <- data$y
  P <- data$P
  C  <- data$C
  if(is.null(data$S)) S <- 0
  else S <- data$S

  C0 <- C
  C0[, P] <- 0
  mu <- pars$phi$mu
  Sigma <- pars$phi$Sigma
  pi <- pars$pi
  # Wrap in 'c()' to drop matrix to vector
  # Should use 'drop()' instead?
  A1 <- c(C %*% mu + S + 0.5 * diag(C %*% Sigma %*% t(C)))
  A0 <- c(C0 %*% mu + S + 0.5 * diag(C0 %*% Sigma %*% t(C0)))
  if(!is.null(data$S)) {
    A1 <- A1 + data$S
    A0 <- A0 + data$S
  }

  d_likelihood <-
    (1 - pi) * t(C) %*% diag(exp(A1)) %*% C +
    pi * t(C0) %*% diag(exp(A0)) %*% C0

  Tau <- pars$theta$Tau
  -0.5 * (Tau + d_likelihood)
}
