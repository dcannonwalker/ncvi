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
#'   updated by the algorithm
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
                                verbose = T,
                                fixed_iter = F),
                 ...) {

  args <- list(...)
  pars = init
  iter = 0
  delta = options$elbo_delta + 1
  L = NULL

  if (options$fixed_iter == F) {
    while (iter < options$max_iter &&
           abs(delta) > options$elbo_delta) {

      pars_old <- pars

      pars <- update_pars(data, pars, args)

      elbo_out <- elbo(data, pars_new = pars, pars_old = pars_old)

      iter <- iter + 1

      L <- elbo_out$L

      delta <- elbo_out$delta

      if (is.na(delta)) {
        print("delta is NA, continuing")
        delta <- options$elbo_delta + 1
      }

      else if (is.nan(delta)) {
        print("delta is NaN, continuing")
        delta <- options$elbo_delta + 1
      }

      else if(is.infinite(delta)) {
        print("delta is infinite, continuing")
        delta <- options$elbo_delta + 1
      }

      if (options$verbose == T) print(c(delta, iter, L))

    }

  }
  else if (options$fixed_iter == T) {
    while (iter < options$max_iter) {

      pars <- update_pars(data, pars, args)

      iter <- iter + 1

    }
  }

  list(pars = pars,
       data = data,
       L = L)

}
