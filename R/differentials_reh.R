#' The contribution of the likelihood to the
#' differential w.r.t. vec of the mvn covariance,
#' for Poisson-distributed data
#'
#' @param data List of observed or known variables, including 'C'.
#' @param pars List of parameters, including `phi$mu`, `phi$Sigma`,
#' @export
d_mvn_cov_reh <- function(data, pars) {
  C <- data$C
  mu <- pars$phi$mu
  Sigma <- pars$phi$Sigma
  # Wrap in 'c()' to drop matrix to vector
  # Should use 'drop()' instead?
  A <- c(C %*% mu + 0.5 * diag(C %*% Sigma %*% t(C)))

  t(C) %*% diag(exp(A)) %*% C
}

#' The contribution of the likelihood to the
#' differential w.r.t. mvn mean,
#' for Poisson-distributed data
#'
#' @param data List of observed or known variables, including 'C'.
#' @param pars List of parameters, including `phi$mu`, `phi$Sigma`,
#' @export
d_mvn_mean_reh <- function(data, pars) {
  C <- data$C
  y <- data$y
  mu <- pars$phi$mu
  Sigma <- pars$phi$Sigma
  A <- c(C %*% mu + 0.5 * diag(C %*% Sigma %*% t(C)))

  t(C) %*% (y - exp(A))
}
