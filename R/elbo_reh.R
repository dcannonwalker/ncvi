# REH elbo

elbo_extra_REH <- function(data, pars, priors = NULL) {

  phi <- pars$phi
  M <- pars$theta$M
  R <- pars$theta$R
  Tau <- pars$theta$Tau
  U <- data$U
  G <- data$G
  precision_mu0 <- pars$theta$precision_mu0
  list_beta <- pars$theta$list_beta
  list_u <- pars$theta$list_u


  collection_mu_i <- lapply(phi, function(x) x$mu)

  sum_mu_i <- Reduce("+", collection_mu_i)

  extra <- -(G * sum((M^2 + c(diag(R), rep(0, U)) -
                        2 * (sum_mu_i * M)) *
                       diag(Tau))) / 2 -
    precision_mu0 * sum(M^2) / 2

  if (!is.null(priors)) {
    extra <- extra + (priors$a_beta - list_beta$a) *
      (digamma(list_beta$a) - log(list_beta$b)) +
      (priors$b_beta - list_beta$b) * list_beta$mean +
      (priors$a_u - list_u$a) *
      (digamma(list_u$a) - log(list_u$b)) +
      (priors$b_u - list_u$b) * list_u$mean
  }
  extra
}

## this is essentially generalized from elbo_i with the
## replacement of
elbo_i_REH <- function(data, pars) {

  y <- data$y
  C <- data$C
  mu <- pars$phi$mu
  Sigma <- pars$phi$Sigma
  Tau <- pars$theta$Tau

  A <- c(C %*% mu + 0.5 * diag(C %*% Sigma %*% t(C)))

  y %*% C %*% mu - sum(exp(A)) -
    sum((mu^2 + diag(Sigma)) * diag(Tau)) / 2 +
    0.5 * log(det(Sigma))
}

#' @export
elbo_REH <- function(data, pars, old_elbo = NULL, priors = NULL) {

  phi <- pars$phi
  theta <- pars$theta
  C <- data$C
  y <- data$y
  G <- data$G
  pars_list <- list()
  for (i in seq(1, length(data$y))) {
    pars_list[[i]] <- list(data = list(y = data$y[[i]], C = C),
                           pars = list(phi = phi[[i]],
                                       theta = theta))

  }
  if (!is.null(priors)) {
    extra <- elbo_extra_REH(data, pars, priors)
  }
  else extra <- elbo_extra_REH(data, pars)
  elbo_list <- lapply(pars_list, function(p) {
    elbo_i_REH(data = p$data,
               pars = p$pars)
  })
  elbo_vector <- unlist(elbo_list)
  elbo_vector[G+1] <- extra

  elbo <- sum(elbo_vector)
  if (!is.null(old_elbo)) {
    list(delta_vector = elbo_vector - old_elbo$elbo_vector,
         delta = sum(elbo_vector - old_elbo$elbo_vector),
         elbo_vector = elbo_vector,
         elbo = elbo)
  }
  else list(elbo_vector = elbo_vector, elbo = elbo)

}
