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
                                    verbose = TRUE,
                                    fixed_iter = FALSE,
                                    short_out = FALSE),
                     ...) {

  args <- list(...)
  pars = init
  if(is.null(pars$problem_index)) pars$problem_index <- c()
  iter = 0
  L = elbo(data, pars, priors = args$priors)
  L$delta <- options$elbo_delta + 1
  maxL_pars <- pars
  if (!is.nan(L$elbo) &
      !is.infinite(L$elbo) &
      !is.na(L$elbo) &
      !is.null(L$elbo)) {
    maxL <- L$elbo
  }
  if (options$fixed_iter == FALSE) {
    while (iter < options$max_iter &&
           abs(L$delta) > options$elbo_delta) {

      pars <- update_pars(data, pars, args)

      iter <- iter + 1

      L <- elbo(data, pars, old_elbo = L, priors = args$priors)

      delta <- L$delta

      if (is.na(delta)) {
        print("delta is NA, continuing")
        L$delta <- options$elbo_delta + 1
      }

      else if (is.nan(delta)) {
        print("delta is NaN, continuing")
        L$delta <- options$elbo_delta + 1
      }

      else if(is.infinite(delta)) {
        print("delta is infinite, continuing")
        L$delta <- options$elbo_delta + 1
      }

      else if(!is.nan(L$elbo) &
              !is.infinite(L$elbo) &
              !is.na(L$elbo) &
              !is.null(L$elbo)) {
        if (L$elbo > maxL) {
          maxL_pars <- pars
          maxL <- L$elbo
        }
      }


      if (options$verbose == TRUE) print(data.frame(iter = iter,
                                                    elbo = L$elbo,
                                                    delta = L$delta,
                                                    max_elbo = maxL))

    }

  }
  else if (options$fixed_iter == TRUE) {
    while (iter < options$max_iter) {

      pars <- update_pars(data, pars, args)

      iter <- iter + 1

    }
  }

  if(options$short_out) {
    list(beta = t(sapply(pars$phi, function(p) p$mu[1:4])),
         pi = pars$pi)
  }

  else {
    list(
      beta = t(sapply(pars$phi, function(p) p$mu[1:4])),
      pi = pars$pi,
      pars = pars,
      maxL_pars = maxL_pars,
      data = data,
      L = L)
  }

}
