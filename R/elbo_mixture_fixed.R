elbo_i_mixture <- function(data, pars) {
  y <- data$y
  C1 <- data$C1
  C0 <- data$C0
  Tau <- pars$theta$Tau # expectation of prior precision for phi
  Sigma0 <- MASS::ginv(Tau)
  M <- pars$theta$M
  R <- pars$theta$R
  pi0 <- pars$theta$pi0
  pi <- pars$pi
  mu <- pars$phi$mu
  Sigma <- pars$phi$Sigma

  #   A <- c(C %*% mu + 0.5 * diag(C %*% Sigma %*% t(C)))
  #   Astar <- c(Cstar %*% mu + 0.5 * diag(Cstar %*% Sigma %*% t(Cstar)))
  A1 <- c(C1 %*% mu + 0.5 * diag(C1 %*% Sigma %*% t(C1)))
  A0 <- c(C0 %*% mu + 0.5 * diag(C0 %*% Sigma %*% t(C0)))

  # contribution of ln p(y | beta, z)
  elbo_i <- (
    pi * (t(y) %*% C0 %*% mu - sum(exp(A0))) +
    (1 - pi) * (t(y) %*% C1 %*% mu - sum(exp(A1))) +

  # contribution of ln p(beta | mu)
    -log(det(Sigma0)) / 2 - sum((mu^2  + diag(Sigma)) * diag(Tau)) / 2 +
    mu %*% Tau %*% M - sum((M^2 + diag(R)) * diag(Tau)) / 2 +

  # contribution of ln p(Z | pi0) - ln q(Z)

    contribution_pi_i(pi = pi, pi0 = pi0)

  # contribution of ln q(beta)
    -log(det(Sigma)) / 2

  )


  elbo_i

}

contribution_pi_i <- function(pi, pi0) {
  if (pi == 0) {
    log(1 - pi0)
  }
  else if (pi == 1) {
    log(pi0)
  }
  else if (pi > 0 && pi < 1) {
    pi * log(pi0) + (1 - pi) * log(1 - pi0) -
      pi * log(pi) - (1 - pi) * log(1 - pi)
  }
}
elbo_extra_mixture <- function(data, pars) {

  M <- pars$theta$M
  R <- pars$theta$R
  precision_mu0 <- pars$theta$precision_mu0

  -precision_mu0 * (sum(M^2) + sum(diag(R))) / 2 -
                      det(R) / 2

}

elbo_mixture_list <- function(data, pars) {

  y <- data$y
  C1 <- data$C1
  C0 <- data$C0
  G <- data$G
  theta <- pars$theta
  pi <- pars$pi
  phi <- pars$phi

  extra <- elbo_extra_mixture(data, pars)

  elbo_i_list <- list()
  for (i in seq(1, G)) {
    elbo_i_list[[i]] <- elbo_i_mixture(data = list(y = y[[i]], C1 = C1,
                                                   C0 = C0),
                                       pars = list(phi = phi[[i]],
                                                   pi = pi[[i]], theta = theta)
    )
  }

  list(extra = extra, elbo_i_list = elbo_i_list)
}

elbo_mixture <- function(data, pars, old_elbo = NULL) {
  elbo_list <- elbo_mixture_list(data, pars)
  elbo_vector <- unlist(elbo_list$elbo_i_list)
  elbo_vector[length(elbo_vector) + 1] <- elbo_list$extra

  elbo <- sum(elbo_vector)
  if (!is.null(old_elbo)) {
    list(delta_vector = elbo_vector - old_elbo$elbo_vector,
         delta = sum(elbo_vector - old_elbo$elbo_vector),
         elbo_vector = elbo_vector,
         elbo = elbo)
  }
  else list(elbo_vector = elbo_vector, elbo = elbo)

}
