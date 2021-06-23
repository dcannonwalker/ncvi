# # elbo functions
# elbo_mixture <- function(data, pars) {
#
#   y <- data$y
#   G <- data$G
#   C <- data$C
#   phi <- pars$phi
#   theta <- pars$theta
#   pi <- pars$pi
#
#   L <- elbo_extra_mixture(data, pars)
#
#   for (i in seq(1, G)) {
#     L <- L + elbo_i_mixture(data = list(y = y[[i]], C = C),
#                             pars = list(phi = phi[[i]], theta = theta, pi = pi[[i]]))
#   }
#
#   L
#
# }
#
# elbo_mixture_list <- function(data, pars) {
#
#   y <- data$y
#   G <- data$G
#   C <- data$C
#   phi <- pars$phi
#   theta <- pars$theta
#   pi <- pars$pi
#
#   elbo_list <- list(
#     extra = elbo_extra_mixture(data, pars)
#   )
#
#   for (i in seq(1, G)) {
#     elbo_list$groupwise[[i]] <-
#       elbo_i_mixture(data = list(y = y[[i]], C = C),
#                      pars = list(phi = phi[[i]],
#                                  theta = theta,
#                                  pi = pi[[i]]))
#   }
#
#   elbo_list
#
# }
#
# elbo_i_mixture <- function(data, pars) {
#
#   y <- data$y
#   C <- data$C
#   Cstar <- C
#   Cstar[, 2] <- 0
#   mu <- pars$phi$mu
#   Sigma <- pars$phi$Sigma
#   Tau <- pars$theta$Tau
#   M <- pars$theta$M
#   pi0 <- pars$theta$pi0
#   pi <- pars$pi
#
#   A <- c(C %*% mu + 0.5 * diag(C %*% Sigma %*% t(C)))
#   Astar <- c(Cstar %*% mu + 0.5 * diag(Cstar %*% Sigma %*% t(Cstar)))
#
#   if (pi != 0 && pi != 1) {
#     (1 - pi) * (sum(exp(Astar)) - y %*% Cstar %*% mu) +
#       pi * (sum(exp(A)) - y %*% C %*% mu) -
#       sum((mu^2 + diag(Sigma) - 2 * mu * M) * diag(Tau)) / 2 +
#       0.5 * log(det(Sigma)) +
#       (1 - pi) * log(1 - pi0) + pi * log(pi0) -
#       (1 - pi) * log(1 - pi) - pi * log(pi)
#   }
#   else if (pi == 0) {
#     (sum(exp(Astar)) - y %*% Cstar %*% mu) +
#       sum((mu^2 + diag(Sigma) - 2 * mu * M) * diag(Tau)) / 2 +
#       0.5 * log(det(Sigma)) +
#       log(1 - pi0)
#   }
#   else if (pi == 1) {
#     (sum(exp(A)) - y %*% C %*% mu) -
#       sum((mu^2 + diag(Sigma) - 2 * mu * M) * diag(Tau)) / 2 +
#       0.5 * log(det(Sigma)) +
#       log(pi0)
#   }
#
# }
#
# elbo_extra_mixture <- function(data, pars) {
#
#   phi <- pars$phi
#   M <- pars$theta$M
#   R <- pars$theta$R
#   Tau <- pars$theta$Tau
#   U <- data$U
#   G <- data$G
#   precision_mu0 <- pars$theta$precision_mu0
#
#   -0.5 * (G * sum((M^2 + c(diag(R), rep(0, U))) * diag(Tau))) -
#     0.5 * precision_mu0 * sum(M^2)
#
# }
