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
                                            d_mvn_cov_poisson)
                                },
                                mu = function(data, pars) {
                                  d_mvn_mean(data, pars,
                                             d_mvn_mean_poisson)
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

elbo_extra_REH <- function(data, pars, priors = NULL) {

  phi <- pars$phi
  M <- pars$theta$M
  R <- pars$theta$R
  Tau <- pars$theta$Tau
  U <- data$U
  G <- data$G
  precision_mu0 <- pars$theta$precision_mu0
  list_beta <- pars$theta$list_beta
  list_u <- pars$theta$list_u


  collection_mu_i <- lapply(phi, function(x) x$mu)

  sum_mu_i <- Reduce("+", collection_mu_i)

  extra <- -(G * sum((M^2 + c(diag(R), rep(0, U)) -
               2 * (sum_mu_i * M)) *
              diag(Tau))) / 2 -
    precision_mu0 * sum(M^2) / 2

  if (!is.null(priors)) {
    extra <- extra + (priors$a_beta - list_beta$a) *
      (digamma(list_beta$a) - log(list_beta$b)) +
      (priors$b_beta - list_beta$b) * list_beta$mean +
      (priors$a_u - list_u$a) *
      (digamma(list_u$a) - log(list_u$b)) +
      (priors$b_u - list_u$b) * list_u$mean
  }
  extra
}

## this is essentially generalized from elbo_i with the
## replacement of
elbo_i_REH <- function(data, pars) {

  y <- data$y
  C <- data$C
  mu <- pars$phi$mu
  Sigma <- pars$phi$Sigma
  Tau <- pars$theta$Tau

  A <- c(C %*% mu + 0.5 * diag(C %*% Sigma %*% t(C)))

  y %*% C %*% mu - sum(exp(A)) -
    sum((mu^2 + diag(Sigma)) * diag(Tau)) / 2 +
    0.5 * log(det(Sigma))
}

#' @export
elbo_REH <- function(data, pars, old_elbo = NULL, priors = NULL) {

  phi <- pars$phi
  theta <- pars$theta
  C <- data$C
  y <- data$y
  G <- data$G
  pars_list <- list()
  for (i in seq(1, length(data$y))) {
    pars_list[[i]] <- list(data = list(y = data$y[[i]], C = C),
                           pars = list(phi = phi[[i]],
                                       theta = theta))

  }
  if (!is.null(priors)) {
    extra <- elbo_extra_REH(data, pars, priors)
  }
  else extra <- elbo_extra_REH(data, pars)
  elbo_list <- lapply(pars_list, function(p) {
    elbo_i_REH(data = p$data,
               pars = p$pars)
  })
  elbo_vector <- unlist(elbo_list)
  elbo_vector[G+1] <- extra

  elbo <- sum(elbo_vector)
  if (!is.null(old_elbo)) {
    list(delta_vector = elbo_vector - old_elbo$elbo_vector,
         delta = sum(elbo_vector - old_elbo$elbo_vector),
         elbo_vector = elbo_vector,
         elbo = elbo)
  }
  else list(elbo_vector = elbo_vector, elbo = elbo)

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
