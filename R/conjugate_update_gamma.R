## updates for Gamma-distributed precision parameters

#' Conjugate update for the parameters of ind. Gamma distributions
#' with ind. mvn child nodes that have a mvn prior mean; prior precision
#' of child nodes
#'
#' @param pars List with `phi, theta`
#' @param priors List with `a, b`
#' @param data List with `S, H`
#' @export
conjugate_update_gamma_precision_matrix <- function(data, pars, priors) {

  S <- data$S
  H <- data$H
  phi <- pars$phi
  M <- S %*% pars$theta$M
  R <- S %*% pars$theta$R

  collection_mu_i <- lapply(phi, function(x) x$mu)

  sum_mu_i <- Reduce("+", collection_mu_i)

  mu <- S %*% sum_mu_i

  collection_diag_Sigma_i <- lapply(phi, function(x) diag(x$Sigma))

  sum_diag_Sigma_i <- S %*% Reduce("+", collection_diag_Sigma_i)

  a <- priors$a + colSums(S)/2

  b <- priors$b + 0.5 * (
    sum(mu^2 + sum + M^2 + (mu * M) + diag(R))
    )

  Tau <- diag(H %*% (a / b))

  # return a, b, Tau
  list(a = a, b = b, Tau = Tau)

}

conjugate_update_gamma_precision <- function() {



  list(t = t)

}

