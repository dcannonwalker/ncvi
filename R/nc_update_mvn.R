#' Non-conjugate update for the parameters of
#' multivariate normal distribution.
#'
#' @param data List of observed and known variables passed down from the
#' `fit_ncvi()` wrapper function
#' @param pars List of unknown variables,
#'   including `phi$mu`
#' @param differentials List of differentials functions, should have
#' 'mu' and 'Sigma'
#' @export
nc_update_mvn <- function(data, pars, differentials, i) {
  Sigma <-
    MASS::ginv(
      -2 * differentials$Sigma(data, pars)
    )
  pars$phi$Sigma <- Sigma
  mu <- c(pars$phi$mu + Sigma %*% differentials$mu(data, pars))

  list(mu = mu, Sigma = Sigma)
}
