## conjugate updates for mvn hierarchical mean
## To do: generalize
## To do: better/shorter function names?
## To do: move explanatory comments out of the function body

#' Conjugate update for the parameters of multivariate normal
#' distribution for the hierarchical mean of
#' multivariate normal child nodes
#'
#' @param pars List with `phi` and `theta`
#' @export
conjugate_update_mvn_hierarchical_mean <- function(data, pars) {

  phi <- pars$phi
  Tau <- pars$theta$Tau
  precision_mu0 <- pars$theta$precision_mu0
  K <- nrow(Tau)
  G <- data$G

  # update the variance

  R <- MASS::ginv(

    # Tau is the variational expectation of the prior precision
    # of the Beta_i,
    # i.e. E(inverse(Sigma_0))
    # precision_mu0 is the variational expectation of the precision param
    # of mu_0,
    # i.e. E(1/sig2)

    G*Tau + (precision_mu0 * diag(K))

  )
  # update the mean
  # collect and sum the variational means of the Beta_i

  collection_mu_i <- lapply(phi, function(x) x$mu)

  sum_mu_i <- Reduce("+", collection_mu_i)

  M <- drop(R %*% (Tau %*% sum_mu_i))

  # return R and M

  list(R = R, M = M)

}

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
