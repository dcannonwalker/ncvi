## Differentials for multivariate normal distribution.

#' Vec inverse of the differential with respect to vec of covariance
#'
#' @param data List of observed or known variables.
#'   If `d_likelihood()` requires only a subset of `data` param
#'   given to `ncvi()`, `data` should be only this subset.
#' @param pars List of parameters, including `theta$Tau`.
#'   If `d_likelihood()` requires only a subset of the `pars` param
#'   given to `ncvi()`, `pars` should be only this subset.
#' @param d_likelihood Function to provide the contribution of the
#'   data likelihood to the differential. Should take params `data`
#'   and `pars` as given to this function.
#' @export
d_mvn_cov <- function(data, pars, d_likelihood) {
  Tau <- pars$theta$Tau
  -0.5 * (Tau + d_likelihood(data, pars))
}

#' The contribution of the likelihood to the
#' differential w.r.t. vec of the mvn covariance,
#' for Poisson-distributed data
#'
#' @param data List of observed or known variables, including 'C'.
#' @param pars List of parameters, including `phi$mu`, `phi$Sigma`,
#' @export
d_mvn_cov_poisson <- function(data, pars) {
  C <- data$C
  mu <- pars$phi$mu
  Sigma <- pars$phi$Sigma
  # Wrap in 'c()' to drop matrix to vector
  # Should use 'drop()' instead?
  A <- c(C %*% mu + 0.5 * diag(C %*% Sigma %*% t(C)))

  t(C) %*% diag(exp(A)) %*% C
}

#'
#' @param data List of observed or known variables, including 'C', 'y'.
#'   If `d_likelihood()` requires only a subset of `data` param
#'   given to `ncvi()`, `data` should be only this subset.
#' @param pars List of parameters, including `theta$Tau`, `theta$M`,
#'   `phi$mu`, `phi$Sigma`.
#'   If `d_likelihood()` requires only a subset of the `pars` param
#'   given to `ncvi()`, `pars` should be only this subset.
#' @param d_likelihood Function to provide the contribution of the
#'   data likelihood to the differential. Should take params `data`
#'   and `pars` as given to this function.
#' @export
d_mvn_mean <- function(data, pars, d_likelihood) {
  Tau <- pars$theta$Tau
  M <- pars$theta$M
  mu <- pars$phi$mu

  Tau %*% (M - mu) + d_likelihood(data, pars)
}

#' The contribution of the likelihood to the
#' differential w.r.t. mvn mean,
#' for Poisson-distributed data
#'
#' @param data List of observed or known variables, including 'C'.
#' @param pars List of parameters, including `phi$mu`, `phi$Sigma`,
#' @export
d_mvn_mean_poisson <- function(data, pars) {
  C <- data$C
  y <- data$y
  mu <- pars$phi$mu
  Sigma <- pars$phi$Sigma
  A <- c(C %*% mu + 0.5 * diag(C %*% Sigma %*% t(C)))

  t(C) %*% (y - exp(A))
}
