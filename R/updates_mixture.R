update_pars_mixture <- function(data, pars, args) {

  pars$phi <- update_phi_mixture(data, pars)

  pars$theta <- update_theta_mixture(data, pars, args$priors)

  pars$pi <- update_pi(data, pars)

  pars
}

# can augment with precision updates

update_theta_mixture <- function(data, pars, priors) {
  theta <- pars$theta
  U <- data$U
  pars_mu0 <- conjugate_update_mvn_REH(data, pars)
  theta$M <- c(pars_mu0$M, rep(0, U))
  theta$R <- pars_mu0$R
  if (!is.null(priors)) {
    list_beta <- update_precision_beta_mixture(data, pars, priors)
    theta$list_beta <- list_beta
    if (U != 0) {
      list_u <- update_precision_u_mixture(data, pars, priors)
      theta$precision_u <- c(list_u$mean)
      theta$list_u <- list_u
    }
    theta$precision_beta <- c(list_beta$mean)
    if (length(list_beta$mean) == 1) {
      theta$Tau <- diag(c(rep(theta$precision_beta, data$P),
                          rep(theta$precision_u, data$U)))
    }

    else if (length(list_beta$mean) == P) {
      theta$Tau <- diag(c(theta$precision_beta),
                        rep(theta$precision_u, data$U))
    }

  }

  theta

}

update_precision_u_mixture <- function(data, pars, priors) {
  phi <- pars$phi
  M <- pars$theta$M
  R <- pars$theta$R
  G <- data$G
  P <- data$P
  U <- data$U
  a_u <- priors$a_u
  b_u <- priors$b_u
  collection_mu_i2 <- lapply(phi,
                             function(x) x$mu[(P + 1):(P + U)]^2)
  sum_mu_i2 <- Reduce("+", collection_mu_i2)
  collection_diag_Sigma_i <- lapply(phi,
                                    function(x) diag(x$Sigma)[(P + 1):(P + U)])
  sum_diag_Sigma_i <- Reduce("+", collection_diag_Sigma_i)

  a <- G * U / 2 + a_u
  b <- b_u + (sum(sum_mu_i2) + sum(sum_diag_Sigma_i)) / 2

  list(mean = a / b, a = a, b = b)

}

update_precision_beta_mixture <- function(data, pars, priors) {
  phi <- pars$phi
  G <- data$G
  P <- data$P
  M <- pars$theta$M[1:P]
  R <- pars$theta$R
  a_beta <- priors$a_beta
  b_beta <- priors$b_beta
  if (length(pars$theta$precision_beta) == 1) {
    collection_mu_i <- lapply(phi, function(x) x$mu[1:P])
    collection_mu_i2 <- lapply(collection_mu_i, function(x) x^2)
    sum_mu_i <- Reduce("+", collection_mu_i)
    sum_mu_i2 <- Reduce("+", collection_mu_i2)
    collection_diag_Sigma_i <- lapply(phi, function(x) diag(x$Sigma)[1:P])
    sum_diag_Sigma_i <- Reduce("+", collection_diag_Sigma_i)

    a <- G * P / 2 + a_beta
    b <- b_beta + (sum(sum_mu_i2) + sum(sum_diag_Sigma_i) -
                     2 * t(M) %*% sum_mu_i + G * sum(M^2) + G * sum(diag(R))) / 2

    list(mean = a / b, a = a, b = b)
  }

  else if (length(pars$theta$precision_beta) == P) {
    collection_mu_i <- lapply(phi, function(x) x$mu[1:P])
    collection_mu_i2 <- lapply(collection_mu_i, function(x) x^2)
    sum_mu_i <- Reduce("+", collection_mu_i)
    sum_mu_i2 <- Reduce("+", collection_mu_i2)
    collection_diag_Sigma_i <- lapply(phi, function(x) diag(x$Sigma)[1:P])
    sum_diag_Sigma_i <- Reduce("+", collection_diag_Sigma_i)

    a <- rep(G / 2 + a_beta, P)
    b <- b_beta + (sum_mu_i2 + sum_diag_Sigma_i -
                     2 * M * sum_mu_i + G * M^2 + G * diag(R)) / 2
  }

  a <- G * P / 2 + a_beta
  b <- b_beta + (sum(sum_mu_i2) + sum(sum_diag_Sigma_i) -
                   2 * t(M) %*% sum_mu_i + G * sum(M^2) + G * sum(diag(R))) / 2

  list(mean = a / b, a = a, b = b)
}

update_phi_mixture <- function(data, pars) {

  phi <- pars$phi
  theta <- pars$theta
  pi <- pars$pi
  y <- data$y
  C <- data$C
  G <- data$G
  P <- data$P

  for (i in seq(1, G)) {

    phi[[i]] <- nc_update_mvn(data = list(y = y[[i]], C = C, P = P),
                              pars = list(phi = phi[[i]], theta = theta,
                                          pi = pi[[i]]),
                              differentials = list(
                                Sigma = d_mvn_cov_mixture,
                                mu = d_mvn_mean_mixture
                              ))

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
  G <- data$G
  P <- data$P

  for (i in seq(1, G)) {

    pi[[i]] <- update_pi_i(data = list(y = y[[i]], C = C, P = P),
                           pars = list(phi = phi[[i]], theta = theta),
                           i = i)

  }

  # return pi
  pi
}

update_pi_i <- function(data, pars, i) {

  mu <- pars$phi$mu
  Sigma <- pars$phi$Sigma
  y <- data$y
  P <- data$P
  C <- data$C
  C0 <- C
  C0[, P] <- 0

  pi0 <- pars$theta$pi0

  A0 <- C0 %*% mu  +
    0.5 * diag(C0 %*% Sigma %*% t(C0))
  A1 <- C %*% mu  + 0.5 * diag(C %*% Sigma %*% t(C))
  B0 <- log(pi0) + t(y) %*% C0 %*% mu -
    sum(exp(A0))
  B1 <- log(1 - pi0) + t(y) %*% C %*% mu -
    sum(exp(A1))

  pi <- c((exp(B1 - B0) + 1)^-1)
  if (is.nan(pi)) {
    message("pi is nan", i)
  }

  #return variational mean
  max(min(pi, 1 - 1e-10), 1e-8)

}



