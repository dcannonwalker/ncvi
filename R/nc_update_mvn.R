#' Non-conjugate update for the parameters of
#' multivariate normal distribution.
#'
#' @param data List of observed and known variables passed down from the
#' ncvi wrapper function.
#' @param pars List of unknown variables passed down from ncvi function.
#' @param Sigma Current variational covariance matrix for the MVN.
#' @param mu Current variational mean vector for the MVN.
#' @param differentials List of differentials functions, should have
#' 'mu' and 'Sigma'.
#' @export

nc_update_mvn <- function(data, pars, Sigma, mu,  differentials) {
  Sigma <- MASS::ginv(
    -2*differential_Sigma(data, pars)
    )
  mu <- mu + Sigma%*%differentials$mu(data, pars)
  list(mu = mu, Sigma = Sigma)
}
