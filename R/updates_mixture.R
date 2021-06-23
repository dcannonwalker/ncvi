update_pars <- function(data, pars, args) {

  pars$phi <- update_phi(data, pars, args$differentials)

  pars$theta <- update_theta(data, pars, args$priors)

  pars$pi <- update_pi(data, pars)

  pars
}

# can augment with precision updates

update_theta <- function(data, pars, priors) {
  theta <- pars$theta
  U <- data$U
  pars_mu0 <- conjugate_update_mvn_REH(data, pars)
  theta$M <- c(pars_mu0$M, rep(0, U))
  theta$R <- pars_mu0$R

  theta

}

update_phi <- function(data, pars, differentials) {

  phi <- pars$phi
  theta <- pars$theta
  pi <- pars$pi
  y <- data$y
  C <- data$C
  C1 <- data$C1
  C0 <- data$C0
  G <- data$G

  for (i in seq(1, G)) {

    phi[[i]] <- nc_update_mvn(data = list(y = y[[i]], C = C, C1 = C1, C0 = C0),
                              pars = list(phi = phi[[i]], theta = theta,
                                          pi = pi[[i]]),
                              differentials = differentials)

  }

  # return phi
  phi

}

update_pi <- function(data, pars) {

  phi <- pars$phi
  theta <- pars$theta
  pi <- pars$pi
  y <- data$y
  C <- data$C
  C1 <- data$C1
  C0 <- data$C0
  G <- data$G

  for (i in seq(1, G)) {

    pi[[i]] <- update_bern(data = list(y = y[[i]], C = C, C1 = C1, C0 = C0),
                           pars = list(phi = phi[[i]], theta = theta))

  }

  # return pi
  pi
}

update_bern <- function(data, pars) {

  phi <- pars$phi
  mu <- phi$mu
  Sigma <- phi$Sigma
  theta <- pars$theta
  y <- data$y
  C1 <- data$C1
  C0 <- data$C0
  pi0 <- theta$pi0
  A0 <- log(pi0) + t(y) %*% C0 %*% mu -
    sum(exp(C0 %*% mu  +
            0.5 * diag(C0 %*% Sigma %*% t(C0))))
  A1 <- log(1 - pi0) + t(y) %*% C1 %*% mu -
    sum(exp(C1 %*% mu  + 0.5 * diag(C1 %*% Sigma %*% t(C1))))

  pi <- (exp(A1 - A0) + 1)^-1

  #return variational mean
  c(pi)

}

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
  pi * t(C0) %*% (y - exp(A0)) + (1 - pi) * t(C1) %*% (y - exp(A1))
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

  pi * t(C0) %*% diag(exp(A0)) %*% C0 +
    (1 - pi) * t(C1) %*% diag(exp(A1)) %*% C1
}

