## Differentials for multivariate normal distribution

#' Vec inverse of the differential with respect to vec of mvn covariance
#'
#' @param data List of observed or known variables;
#'   if `d_likelihood()` requires only a subset of `data` param
#'   given to `fit_ncvi()`, `data` should be only this subset
#' @param pars List of parameters, including `theta$Tau`;
#'   if `d_likelihood()` requires only a subset of the `pars` param
#'   given to `fit_ncvi()`, `pars` should be only this subset
#' @param d_likelihood Function to provide the contribution of the
#'   data likelihood to the differential; should take params `data`
#'   and `pars` as given to this function
#' @export
d_mvn_cov <- function(data, pars, d_likelihood) {
  Tau <- pars$theta$Tau
  -0.5 * (Tau + d_likelihood(data, pars))
}


#' Transpose of the differential w.r.t. mvn mean
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




