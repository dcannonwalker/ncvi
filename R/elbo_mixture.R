elbo_extra_mixture <- function(data, pars, priors = NULL) {

  problem_index <- pars$problem_index
  phi <- pars$phi
  if(length(problem_index) > 0) phi <- phi[-problem_index]
  M <- pars$theta$M
  R <- pars$theta$R
  Tau <- pars$theta$Tau
  U <- data$U
  G <- data$G
  precision_mu0 <- pars$theta$precision_mu0
  list_beta <- pars$theta$list_beta
  if (U != 0) list_u <- pars$theta$list_u


  collection_mu_i <- lapply(phi, function(x) x$mu)

  sum_mu_i <- Reduce("+", collection_mu_i)

  extra <- -(G * sum((M^2 + c(diag(R), rep(0, U)) -
                        2 * (sum_mu_i * M)) *
                       diag(Tau))) / 2 -
    precision_mu0 * sum(M^2) / 2

  if (!is.null(priors)) {
    if (length(pars$theta$precision_beta) == 1) {
      extra <- extra + (priors$a_beta - list_beta$a) *
        (digamma(list_beta$a) - log(list_beta$b)) +
        (priors$b_beta - list_beta$b) * list_beta$mean
    }

    else if (length(pars$theta$precision_beta) == data$P) {
      extra <- extra + sum((priors$a_beta - list_beta$a) *
        (digamma(list_beta$a) - log(list_beta$b)) +
        (priors$b_beta - list_beta$b) * list_beta$mean)
    }


    if (U != 0) {
        extra <- extra + (priors$a_u - list_u$a) *
          (digamma(list_u$a) - log(list_u$b)) +
          (priors$b_u - list_u$b) * list_u$mean
        }

  }
  extra
}

## this is essentially generalized from elbo_i with the
## replacement of
elbo_i_mixture <- function(data, pars) {

  if(is.null(data) | is.null(pars)) return(0)
  y <- data$y
  C <- data$C
  P <- data$P
  C0 <- C
  C0[, P] <- 0
  mu <- pars$phi$mu
  Sigma <- pars$phi$Sigma
  Tau <- pars$theta$Tau
  pi <- pars$pi
  pi0 <- pars$theta$pi0

  A1 <- c(C %*% mu + 0.5 * diag(C %*% Sigma %*% t(C)))
  A0 <- c(C0 %*% mu + 0.5 * diag(C0 %*% Sigma %*% t(C0)))

  (1 - pi) * (y %*% C %*% mu - sum(exp(A1))) -
    sum((mu^2 + diag(Sigma)) * diag(Tau)) / 2 +
    0.5 * log(det(Sigma)) +
    pi * (y %*% C0 %*% mu - sum(exp(A0))) -
    sum((mu^2 + diag(Sigma)) * diag(Tau)) / 2 +
    0.5 * log(det(Sigma)) +
    contribution_pi_i(pi = pi, pi0 = pi0)
}

contribution_pi_i <- function(pi, pi0) {
  if (is.nan(pi)) {
    message("pi is NaN")
    0.5
  }
  else if (pi == 0) {
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

#' @export
elbo_mixture <- function(data, pars, old_elbo = NULL, priors = NULL) {

  # exclude problem phi_i
  problem_index <- pars$problem_index
  index <- seq(1, data$G)
  if(length(problem_index) > 0) {
    index <- index[-problem_index]
  }

  phi <- pars$phi
  theta <- pars$theta
  pi <- pars$pi
  C <- data$C
  y <- data$y
  G <- data$G
  P <- data$P
  pars_list <- list()
  for (i in index) {
    pars_list[[i]] <- list(data = list(y = data$y[[i]], C = C, P = P),
                           pars = list(phi = phi[[i]],
                                       pi = pi[[i]],
                                       theta = theta))

  }
  if (!is.null(priors)) {
    extra <- elbo_extra_mixture(data, pars, priors)
  }
  else extra <- elbo_extra_mixture(data, pars)
  elbo_list <- lapply(pars_list, function(p) {
    elbo_i_mixture(data = p$data,
               pars = p$pars)
  })
  elbo_vector <- unlist(elbo_list)
  elbo_vector[G+1] <- extra


  if (sum(is.nan(elbo_vector)) == 0 & sum(is.infinite(elbo_vector)) == 0) {
    elbo <- sum(elbo_vector)
    if (!is.null(old_elbo)) {
      list(delta_vector = elbo_vector - old_elbo$elbo_vector,
           delta = sum(elbo_vector - old_elbo$elbo_vector),
           elbo_vector = elbo_vector,
           elbo = elbo)
    }
    else list(elbo_vector = elbo_vector, elbo = elbo)
  }
  else if (sum(is.nan(elbo_vector)) > 0 | sum(is.infinite(elbo_vector) > 0)) {
    which_nan = which(is.nan(elbo_vector))
    which_inf = which(is.infinite(elbo_vector))
    if (!is.null(old_elbo)) {
      filter <- !is.nan(elbo_vector) & !is.nan(old_elbo$elbo_vector) &
        !is.infinite(elbo_vector) & !is.infinite(old_elbo$elbo_vector)
      delta_vector <- elbo_vector[filter] - old_elbo$elbo_vector[filter]
      list(elbo_vector = elbo_vector,
           elbo = sum(elbo_vector[filter]),
           delta_vector = delta_vector,
           delta = sum(delta_vector),
           which_nan = which_nan,
           which_inf = which_inf)
    }
    else {
      filter <- !is.nan(elbo_vector)
      list(elbo_vector = elbo_vector, elbo = sum(elbo_vector[filter]),
              which_nan = which_nan,
           which_inf = which_inf)
    }
  }


}

