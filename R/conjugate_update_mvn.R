## conjugate updates for mvn hierarchical mean
## To do: generalize
## To do: better/shorter function names?

#' Conjugate update for the parameters of multivariate normal
#' distribution for the hierarchical mean of
#' multivariate normal child nodes
#'
#' @param pars List with `phi` and `theta`
#' @export
conjugate_update_mvn_hierarchical_mean <- function(pars) {
  # update the variance
  phi <- pars$phi
  Tau <- pars$theta$Tau
  t <- pars$theta$t
  P <- nrow(Tau)

  R <- MASS::ginv(

    # Tau is the variational expectation of the prior precision
    # of the Beta_i,
    # i.e. E(inverse(Sigma_0))
    # t is the variational expectation of the precision param
    # of mu_0,
    # i.e. E(1/sig2)

    Tau + (t * diag(P))

  )
  # update the mean
  # collect and sum the variational means of the Beta_i

  collection_mu_i <- lapply(phi, function(x) x$mu)

  sum_mu_i <- Reduce("+", collection_mu_i)

  M <- R %*% sum(Tau %*% sum_mu_i)

}
