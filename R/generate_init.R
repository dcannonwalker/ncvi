#' Generate a set of inits for hierarchical model, with random effects
#'
#' @param G
#' @param N
#' @param P
#' @param U
#' @param precision_mu0
#' @param precision_beta
#' @param precision_u
#' @param mu_adjust
generate_inits <- function(G, N, P, U, precision_mu0, precision_beta,
                           precision_u, mu_adjust = 0) {
  init_phi <- list()

  for (i in seq(1, G)) {

    init_phi[[i]] <- list(
      mu = rep(mu_adjust, (P + U)),
      Sigma = diag((P + U)))

  }

  init_theta <- list(
    M = rep(0, (P + U)),
    R = diag(P),
    Tau = diag(c(rep(precision_beta, P), rep(precision_u, U))),
    precision_beta = precision_beta,
    precision_u = precision_u,
    precision_mu0 = precision_mu0
  )

  list(
    phi = init_phi,
    theta = init_theta
  )
}
