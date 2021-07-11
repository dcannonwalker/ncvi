update_phi_REH <- function(data, pars) {

  phi <- pars$phi
  theta <- pars$theta
  y <- data$y
  C <- data$C
  G <- data$G

  for (i in seq(1, G)) {

    phi[[i]] <- nc_update_mvn(data = list(y = y[[i]], C = C),
                              pars = list(phi = phi[[i]],
                                          theta = theta),
                              differentials = list(
                                Sigma = function(data, pars) {
                                  d_mvn_cov(data, pars,
                                            d_mvn_cov_reh)
                                },
                                mu = function(data, pars) {
                                  d_mvn_mean(data, pars,
                                             d_mvn_mean_reh)
                                }
                              ))

  }

  # return phi
  phi

}

update_theta_REH <- function(data, pars, priors) {
  theta <- pars$theta
  U <- data$U
  pars_mu0 <- conjugate_update_mvn_REH(data, pars)
  #theta$Tau <- conjugate_update_gamma_precision_matrix()
  theta$M <- c(pars_mu0$M, rep(0, U))
  theta$R <- pars_mu0$R
  if (!is.null(priors)) {
    list_beta <- update_precision_beta_REH(data, pars, priors)
    list_u <- update_precision_u_REH(data, pars, priors)
    theta$precision_u <- c(list_u$mean)
    theta$precision_beta <- c(list_beta$mean)
    theta$Tau <- diag(c(rep(theta$precision_beta, data$P),
                        rep(theta$precision_u, data$U)))
    theta$list_beta <- list_beta
    theta$list_u <- list_u
  }

  theta

}

update_pars_REH <- function(data, pars, args) {
  pars$phi <- update_phi_REH(data, pars)
  pars$theta <- update_theta_REH(data, pars, args$priors)

  pars
}

update_precision_beta_REH <- function(data, pars, priors) {
  phi <- pars$phi
  G <- data$G
  P <- data$P
  M <- pars$theta$M[1:P]
  R <- pars$theta$R
  a_beta <- priors$a_beta
  b_beta <- priors$b_beta
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

update_precision_u_REH <- function(data, pars, priors) {
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
