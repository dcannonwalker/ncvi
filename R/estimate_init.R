estimate_init <- function(data, theta) {
  init_phi <- list()
  init_pi <- c()
  estmu <- matrix(nrow = data$G, ncol = data$P)
  for (i in seq(1, data$G)) {
    init_phi[[i]] <- list(mu = NULL, Sigma = NULL)
    est <- est_init_i(data$y[[i]], data$C, data$P)
    init_phi[[i]]$mu <- est$mu
    estmu[i, ] <- init_phi[[i]]$mu
    init_pi[i] <- theta$pi0
  }
  test <- apply(estmu, 2, sd)
  for (i in seq(1, data$G)) init_phi[[i]]$Sigma <- diag(test)
  theta$M <- apply(estmu, 2, mean)
  list(phi = init_phi, theta = theta, pi = init_pi)
}

est_init_i <- function(y, X, P) {
  fit <- glm(y ~ 0 + X, family = "poisson")
  list(
    p = summary(fit)$coefficients[P, 4],
    mu = coef(fit)
  )

}

estimate_init_mu0 <- function(data) {

}
