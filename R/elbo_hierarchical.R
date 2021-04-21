## To do: generalize

#' Contribution to ELBO of terms not included in group-wise ELBO
#'
#' @inheritParams elbo
elbo_extra <- function(data, pars) {

  phi <- pars$phi

  r <- pars$theta$r

  t <- pars$theta$t

  G <- data$G

  collection_mu_i <- lapply(phi, function(x) x$mu)

  sum_mu_i <- Reduce("+", collection_mu_i)

  -0.5 * (r + (G * t)) * sum(M^2) + t * sum_mu_i %*% M
}

#' Contribution to EBLO from a single group
#'
#' @inheritParams elbo
elbo_i <- function(data, pars) {

  y <- data$y
  C <- data$C
  mu <- pars$phi$mu
  Sigma <- pars$phi$Sigma
  sig2 <- 1 / pars$theta$r

  A <- c(C %*% mu + 0.5 * diag(C %*% Sigma %*% t(C)))


  # returns:
  y %*% C %*% mu - sum(exp(A)) -
    (sum(mu^2) + sum(diag(Sigma))) / (2 * sig2) +
    0.5 * log(det(Sigma))

}

#' Complete ELBO function
#'
#' @param data
#' @param pars
elbo_hierarchical <- function(data, pars) {

  phi <- pars$phi
  theta <- pars$theta
  t <- pars$theta$t
  C <- data$C
  y <- data$y
  G <- data$G

  L <- elbo_extra(data, pars)

  for (i in seq(1,G)) {
    L <- L + elboi(data = list(y = y[[i]], C = C),
                   pars = list(phi = phi[[i]], theta = theta))
  }

  L

}
