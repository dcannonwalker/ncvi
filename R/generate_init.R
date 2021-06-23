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
generate_inits_REH <- function(G, N, P, U, precision_mu0, precision_beta,
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

generate_inits_mixtureFixedOnly <- function(settings) {
  P <- settings$P
  N <- settings$N
  G <- settings$G
  X <- settings$X
  M <- settings$M
  precision_mu0 <- settings$precision_mu0
  precision_beta <- settings$precision_beta
  pi0 <- settings$pi0

  mu0 <- settings$mu0
  if (is.null(mu0)) mu0 <- M

  inits_phi <- list()
  inits_pi <- list()

  for (i in seq(1, G)) {
    inits_phi[[i]] <- list(mu = mu0, Sigma = diag(1 / precision_beta, P))
    inits_pi[[i]] <- pi0
  }

  inits_theta <- list(
    M = mu0,
    R = diag(1 / precision_mu0, P),
    Tau = diag(rep(precision_beta, P)),
    precision_beta = precision_beta,
    precision_mu0 = precision_mu0,
    pi0 = pi0
  )

  init <- list(phi = inits_phi,
               theta = inits_theta,
               pi = inits_pi)


}

generate_inits_mixture <- function(settings) {
  P <- settings$P
  U <- settings$U
  N <- settings$N
  G <- settings$G
  M <- settings$M

  mu0 <- settings$mu0
  if (is.null(mu0)) mu0 <- M

  precision_mu0 <- settings$precision_mu0
  precision_beta <- settings$precision_beta
  precision_u <- settings$precision_u
  pi0 <- settings$pi0

  inits_phi <- list()
  inits_pi <- list()

  for (i in seq(1, G)) {
    inits_phi[[i]] <- list(mu = c(mu0, rep(0, U)),
                           Sigma = diag(c(rep(1 / precision_beta, P),
                                        rep(1 / precision_u, U))))
    inits_pi[[i]] <- pi0
  }

  inits_theta <- list(
    R = diag(1 / precision_mu0, P),
    precision_beta = precision_beta,
    precision_mu0 = precision_mu0,
    pi0 = pi0,
    M = c(mu0, rep(0, U)),
    Tau = diag(c(rep(precision_beta, P), rep(precision_u, U))),
    precision_u = precision_u
  )

  init <- list(phi = inits_phi,
               theta = inits_theta,
               pi = inits_pi)


}
