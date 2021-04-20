## To do: consider change to `update_pars()` function in place of
## the two individual update functions.

#' Non-conjugate variational inference
#'
#' General form of wrapper function for a
#' variational inference algorithm that may include
#' both conjugate and non-conjugate conditional distributions.
#'
#' @param data List of observed or known variables used by the
#'   differentials or update functions
#' @param init List of initial values for unknown variables
#'   updated by the algorithm; should have 'phi' and 'theta'
#' @param update_phi Function to carry out non-conjugate updates of
#'   'phi' parameter
#' @param update_theta Function to carry out conjugate updates of 'theta'
#'   parameter
#' @param elbo Function to calculate the ELBO for the model
#' @param options List of options.
#'   Should include positive real 'elbo_delta', the threshold
#'   change in ELBO to terminate the algorithm
#' @param ... Additional arguments to be passed to `update_phi()`,
#'   e.g. `differentials`
#' @export

fit_ncvi <- function(data, init,
                 update_phi,
                 update_theta,
                 elbo,
                 options = list(max_iter = 100,
                                elbo_delta = 0.0001,
                                verbose = T),
                 ...) {

  args <- list(...)
  L = elbo(data, init)
  pars = init
  iter = 0
  delta = options$elbo_delta + 1


  while (iter < options$max_iter &&
         delta > options$elbo_delta) {

    pars$phi <- update_phi(data, pars, args$differentials)

    pars$theta <- update_theta(data, pars, args$priors)

    delta <- abs(L - (L <- elbo(data,pars)))

    iter <- iter + 1

    if (options$verbose == T) print(c(L, iter))

  }

  pars

}
