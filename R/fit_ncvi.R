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
#' @param update_pars Function to carry out the update to
#'   `pars` at each step
#' @param elbo Function to calculate the ELBO for the model
#' @param options List of options.
#'   Should include positive real 'elbo_delta', the threshold
#'   change in ELBO to terminate the algorithm
#' @param ... Additional arguments to be passed to `update_pars()`,
#'   e.g. `differentials`
#' @export

fit_ncvi <- function(data, init,
                 update_pars,
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

    pars <- update_pars(data, pars, args)

    delta <- abs(L - (L <- elbo(data,pars)))

    iter <- iter + 1

    if (options$verbose == T) print(c(L, iter))

  }

  list(pars = pars,
       data = data,
       elbo = L)

}
