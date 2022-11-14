## To do: this probably doesn't need to be its own function
##  can probably just incorporate into the update_theta() function
##  for this example
#' Conjugate update mu0 in the REH example
#'
#' @param data List with `G, P`
#' @param pars List with `phi` and `theta`
#' @export
conjugate_update_mvn_REH <- function(data, pars) {

  phi <- pars$phi
  P <- data$P
  G <- data$G
  U <- data$U
  Tau <- pars$theta$Tau[1:P, 1:P]
  precision_mu0 <- pars$theta$precision_mu0

  # update the variance

  R <- MASS::ginv(

    # Tau is the variational expectation of the prior precision
    # of the Beta_i,
    # i.e. E(inverse(Sigma_0))
    # t is the variational expectation of the precision param
    # of mu_0,
    # i.e. E(1/sig2)

    G*Tau + (precision_mu0 * diag(P))

  )
  # update the mean
  # collect and sum the variational means of the Beta_i

  collection_mu_i <- lapply(phi, function(x) x$mu[1:P])

  sum_mu_i <- Reduce("+", collection_mu_i)

  M <- drop(R %*% (Tau %*% sum_mu_i))

  # return R and M

  list(R = R, M = M)

}
