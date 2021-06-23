## To do: generalize
## To do: fix internal precision pars names

#' Contribution to ELBO of terms not included in group-wise ELBO
#'
#' @inheritParams elbo_hierarchical
#' @export
elbo_extra <- function(data, pars) {

  phi <- pars$phi
  r <- pars$theta$precision_beta
  t <- pars$theta$precision_mu0
  M <- pars$theta$M
  G <- data$G

  collection_mu_i <- lapply(phi, function(x) x$mu)

  sum_mu_i <- Reduce("+", collection_mu_i)
  # check that the precision pars are in the correct places?
  -0.5 * (t + (G * r)) * sum(M^2) + r * sum_mu_i %*% M
}

#' Contribution to EBLO from a single group
#'
#' @inheritParams elbo_hierarchical
#' @export
elbo_i <- function(data, pars) {

  y <- data$y
  C <- data$C
  mu <- pars$phi$mu
  Sigma <- pars$phi$Sigma
  sig2 <- 1 / pars$theta$precision_beta

  A <- c(C %*% mu + 0.5 * diag(C %*% Sigma %*% t(C)))


  # returns:
  y %*% C %*% mu - sum(exp(A)) -
    (sum(mu^2) + sum(diag(Sigma))) / (2 * sig2) +
    0.5 * log(det(Sigma))

}

#' Complete ELBO function
#'
#' @param data List of `y, C, G`
#' @param pars List of `phi, theta`
#' @export
elbo_hierarchical <- function(data, pars) {

  phi <- pars$phi
  theta <- pars$theta
  C <- data$C
  y <- data$y
  G <- data$G

  L <- elbo_extra(data, pars)

  for (i in seq(1, G)) {
    L <- L + elbo_i(data = list(y = y[[i]], C = C),
                   pars = list(phi = phi[[i]], theta = theta))
  }

  L

}
