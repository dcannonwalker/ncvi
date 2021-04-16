#' Non-conjugate variational inference
#'
#' General form of wrapper function for a
#' variational inference algorithm that may include
#' both conjugate and non-conjugate conditional distributions.
#' @param data List of observed or known variables used by the
#' differentials or update functions.
#' @param init List of initial values for unknown variables
#' updated by the algorithm.
#' @param differentials List of functions for differentials used in
#' non-conjugate updates of 'phi' parameters.
#' @param update_phi Function to carry out non-conjugate updates of
#' 'phi' parameter. Should take args 'data', 'pars', and 'differentials'.
#' @param update_theta Function to carry out conjugate updates of 'theta'
#' parameter. Should take args 'data' and 'pars'.
#' @param elbo Function to calculate the ELBO for the model.
#' @param options List of options.
#' Should include positive real 'elbo_delta', the threshold
#' change in ELBO to terminate the algorithm.
#' @export

ncvi <- function(data, init, differentials,
                  update_phi, update_theta,
                  elbo, options){
  L = elbo(data, init)
  pars = init
  iter = 0
  delta = options$elbo_delta + 1
  while(delta > options$elbo_delta){
    pars$phi <- update_phi(data, pars, differentials)
    pars$theta <- update_theta(data, pars)
    delta <- abs(L - (L <- elbo(data,pars)))
    iter <- iter + 1
    print(c(L, iter))
  }
  return(pars)
}
