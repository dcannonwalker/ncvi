update_phi <- function(data, pars, differentials) {

  phi <- pars$phi
  theta <- pars$theta
  pi <- pars$pi
  y <- data$y
  C <- data$C
  G <- data$G

  for (i in seq(1, G)) {

    phi[[i]] <- nc_update_mvn(data = list(y = y[[i]], C = C),
                              pars = list(phi = phi[[i]], theta = theta),
                              differentials = differentials)

  }

  # return phi
  phi

}

# Transpose of differential wrt mvn mean from poisson likelihood,
# where the second fixed effect parameter has a mixture prior of
# Gaussian and point mass at zero.
d_mixture_mean_poisson <- function(data, pars) {
  C <- data$C
  Cstar <- C # pros/cons move Cstar into 'data'?
  Cstar[, 2] <- 0
  y <- data$y
  mu <- pars$phi$mu
  Sigma <- pars$phi$Sigma
  pi <- pars$phi$pi
  # Wrap in 'c()' to drop matrix to vector
  A <- c(C %*% mu + 0.5 * diag(C %*% Sigma %*% t(C)))
  Astar <- c(Cstar %*% mu + 0.5 * diag(Cstar %*% Sigma %*% t(Cstar)))

  pi * t(C) %*% (y - exp(A)) + (1 - pi) * t(Cstar) %*% (y - exp(Astar))
}

# Vec inverse of differential wrt vec of mvn covariance from poisson likelihood,
# where the second fixed effect parameter has a mixture prior of
# Gaussian and point mass at zero.
d_mixture_cov_poisson <- function(data, pars) {
  C <- data$C
  Cstar <- C
  Cstar[, 2] <- 0
  mu <- pars$phi$mu
  Sigma <- pars$phi$Sigma
  pi <- pars$phi$pi
  # Wrap in 'c()' to drop matrix to vector
  # Should use 'drop()' instead?
  A <- c(C %*% mu + 0.5 * diag(C %*% Sigma %*% t(C)))
  Astar <- c(Cstar %*% mu + 0.5 * diag(Cstar %*% Sigma %*% t(Cstar)))

  pi * t(C) %*% diag(exp(A)) %*% C +
    (1 - pi) * t(Cstar) %*% diag(exp(Astar)) %*% Cstar
}

