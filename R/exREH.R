update_phi_REH <- update_phi_new

update_theta_REH <- function(data, pars, priors) {
  theta <- pars$theta
  U <- data$U
  pars_mu0 <- conjugate_update_mvn_REH(data, pars)
  theta$M <- c(pars_mu0$M, rep(0, U))
  theta$R <- pars_mu0$R

  theta

}

update_pars_REH <- function(data, pars, args) {
  pars$phi <- update_phi_REH(data, pars, args$differentials)
  pars$theta <- update_theta_REH(data, pars, args$priors)

  pars
}

elbo_extra_REH <- function(data, pars) {

  phi <- pars$phi
  M <- pars$theta$M
  R <- pars$theta$R
  Tau <- pars$theta$Tau
  U <- data$U
  G <- data$G
  precision_mu0 <- pars$theta$precision_mu0


  collection_mu_i <- lapply(phi, function(x) x$mu)

  sum_mu_i <- Reduce("+", collection_mu_i)

  -(G * sum((M^2 + c(diag(R), rep(0, U)) -
               2 * (sum_mu_i * M)) *
              diag(Tau))) / 2 -
    precision_mu0 * sum(M^2) / 2




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

elbo_REH <- function(data, pars) {

  phi <- pars$phi
  theta <- pars$theta
  C <- data$C
  y <- data$y
  G <- data$G

  L <- elbo_extra_REH(data, pars)

  for (i in seq(1,G)) {
    L <- L + elbo_i_REH(data = list(y = y[[i]], C = C),
                        pars = list(phi = phi[[i]], theta = theta))
  }

  L

}
