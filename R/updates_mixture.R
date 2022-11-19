update_pars_mixture <- function(data, pars, args) {

  problem_index <- pars$problem_index
  if(is.null(problem_index)) problem_index <- c()

  pars$phi <- update_phi_mixture(data, pars, problem_index)

  problem_index <- pars$phi$problem_index
  pars$phi <- pars$phi$phi

  pars$theta <- update_theta_mixture(data, pars, args$priors,
                                     problem_index)

  pars$pi <- update_pi(data, pars, problem_index)

  pars$problem_index <- problem_index

  pars
}

# can augment with precision updates

update_theta_mixture <- function(data, pars, priors, problem_index) {
  theta <- pars$theta
  U <- data$U

  # exclude problem phi_i
  if(length(problem_index) > 0) {
    pars$phi <- pars$phi[-problem_index]
    data$G <- length(pars$phi)
    pars$pi <- pars$pi[-problem_index]

  }

  pars_mu0 <- conjugate_update_mvn_REH(data, pars)
  theta$M <- c(pars_mu0$M, rep(0, U))
  theta$R <- pars_mu0$R
  # adjust for removed genes


  # condition on which priors where given
  if (!is.null(priors)) {
    list_beta <- update_precision_beta_mixture(data, pars, priors)
    theta$list_beta <- list_beta
    if (U != 0) {
      list_u <- update_precision_u_mixture(data, pars, priors)
      theta$precision_u <- c(list_u$mean)
      theta$list_u <- list_u
    }
    theta$precision_beta <- c(list_beta$mean)
    if (length(list_beta$mean) == 1) {
      theta$Tau <- diag(c(rep(theta$precision_beta, data$P),
                          rep(theta$precision_u, data$U)))
    }

    else if (length(list_beta$mean) == data$P) {
      theta$Tau <- diag(c(theta$precision_beta,
                          rep(theta$precision_u, data$U)))
    }

    if (!is.null(priors$list_pi0)) {
      list_pi0 <- update_pi0(data, pars, priors)
      theta$pi0 <- c(list_pi0$mean)
      theta$list_pi0 <- list_pi0
    }

  }

  theta

}

update_precision_u_mixture <- function(data, pars, priors) {
  phi <- pars$phi
  M <- pars$theta$M
  R <- pars$theta$R
  G <- data$G
  P <- data$P
  U <- data$U
  a_u <- priors$a_u
  b_u <- priors$b_u
  collection_mu_i2 <- lapply(phi,
                             function(x) x$mu[(P + 1):(P + U)]^2)
  sum_mu_i2 <- Reduce("+", collection_mu_i2)
  collection_diag_Sigma_i <- lapply(phi,
                                    function(x) diag(x$Sigma)[(P + 1):(P + U)])
  sum_diag_Sigma_i <- Reduce("+", collection_diag_Sigma_i)

  a <- G * U / 2 + a_u
  b <- b_u + (sum(sum_mu_i2) + sum(sum_diag_Sigma_i)) / 2

  list(mean = a / b, a = a, b = b)

}

update_precision_beta_mixture <- function(data, pars, priors) {
  phi <- pars$phi
  G <- data$G
  P <- data$P
  M <- pars$theta$M[1:P]
  R <- pars$theta$R
  a_beta <- priors$a_beta
  b_beta <- priors$b_beta
  if (length(pars$theta$precision_beta) == 1) {
    collection_mu_i <- lapply(phi, function(x) x$mu[1:P])
    collection_mu_i2 <- lapply(collection_mu_i, function(x) x^2)
    sum_mu_i <- Reduce("+", collection_mu_i)
    sum_mu_i2 <- Reduce("+", collection_mu_i2)
    collection_diag_Sigma_i <- lapply(phi, function(x) diag(x$Sigma)[1:P])
    sum_diag_Sigma_i <- Reduce("+", collection_diag_Sigma_i)

    a <- G * P / 2 + a_beta
    b <- b_beta + (sum(sum_mu_i2) + sum(sum_diag_Sigma_i) -
                     2 * t(M) %*% sum_mu_i + G * sum(M^2) +
                     G * sum(diag(R))) / 2

    list(mean = a / b, a = a, b = b)
  }

  else if (length(pars$theta$precision_beta) == P) {
    collection_mu_i <- lapply(phi, function(x) x$mu[1:P])
    collection_mu_i2 <- lapply(collection_mu_i, function(x) x^2)
    sum_mu_i <- Reduce("+", collection_mu_i)
    sum_mu_i2 <- Reduce("+", collection_mu_i2)
    collection_diag_Sigma_i <- lapply(phi, function(x) diag(x$Sigma)[1:P])
    sum_diag_Sigma_i <- Reduce("+", collection_diag_Sigma_i)

    a <- G / 2 + a_beta
    b <- b_beta + (sum_mu_i2 + sum_diag_Sigma_i -
                     2 * M * sum_mu_i + G * M^2 + G * diag(R)) / 2
  }

  list(mean = a / b, a = a, b = b)
}

update_phi_mixture <- function(data, pars, problem_index) {

  phi <- pars$phi
  theta <- pars$theta
  pi <- pars$pi
  y <- data$y
  C <- data$C
  G <- data$G
  P <- data$P
  if(is.null(data$S)) S <- 0
  else S <- data$S

  index <- seq(1, G)
  if(length(problem_index) > 0) {
    index <- index[-problem_index]
  }
  for (i in index) {
    if(length(phi[[i]]) > 1) {
      problem_phii <- TRUE
      tryCatch(
        error = function(cnd) {
          message(paste0("Problem for phi ", i))
        },
        expr = {
          phi[[i]] <- nc_update_mvn(data = list(y = y[[i]], S = S,
                                                C = C, P = P),
                                    pars = list(phi = phi[[i]],
                                                theta = theta,
                                                pi = pi[[i]]),
                                    differentials = list(
                                      Sigma = d_mvn_cov_mixture,
                                      mu = d_mvn_mean_mixture
                                    ), i = i)
          problem_phii <- FALSE
        }, finally = {
          if(problem_phii) problem_index <- unique(c(problem_index, i))
        }
      )
    }
  }

  # return phi
  list(phi = phi, problem_index = problem_index)

}

update_pi <- function(data, pars, problem_index) {

  phi <- pars$phi
  theta <- pars$theta
  pi <- pars$pi
  y <- data$y
  C <- data$C
  G <- data$G
  P <- data$P
  if(is.null(data$S)) S <- 0
  else S <- data$S

  index <- seq(1, G)
  if(length(problem_index) > 0) {
    index <- index[-problem_index]
  }
  for (i in index) {
    pi[[i]] <- update_pi_i(data = list(y = y[[i]], C = C, P = P,
                                       S = S,
                                       min_pi = data$min_pi),
                           pars = list(phi = phi[[i]], theta = theta),
                           i = i)

  }

  # return pi
  pi
}

update_pi_i <- function(data, pars, i) {

  min_pi <- data$min_pi
  mu <- pars$phi$mu
  Sigma <- pars$phi$Sigma
  y <- data$y
  P <- data$P
  C <- data$C
  if(is.null(data$S)) S <- 0
  else S <- data$S

  C0 <- C
  C0[, P] <- 0

  pi0 <- pars$theta$pi0

  A0 <- C0 %*% mu  + S +
    0.5 * diag(C0 %*% Sigma %*% t(C0))
  A1 <- C %*% mu  + S +
    0.5 * diag(C %*% Sigma %*% t(C))
  B0 <- log(pi0) + t(y) %*% (C0 %*% mu + S) -
    sum(exp(A0))
  B1 <- log(1 - pi0) + t(y) %*% (C %*% mu + S)-
    sum(exp(A1))

  if(is.nan(B1 - B0)) {
    A0 <- Rmpfr::mpfr(C0 %*% mu + S +
                        0.5 * diag(C0 %*% Sigma %*% t(C0)), 10)
    A1 <- Rmpfr::mpfr(C %*% mu + S +
                        0.5 * diag(C %*% Sigma %*% t(C)), 10)
    B0 <- log(pi0) + t(y) %*% (C0 %*% mu + S) -
      sum(exp(A0))
    B1 <- log(1 - pi0) + t(y) %*% (C %*% mu + S)-
      sum(exp(A1))
  }

  pi <- c((exp(B1 - B0) + 1)^-1)
  if (is.nan(pi)) {
    message("pi is nan", i)
    return(1)
  }

  #return variational mean
  Rmpfr::asNumeric(max(min(pi, 1 - min_pi), min_pi))

}

update_pi0 <- function(data, pars, priors) {
  pi <- pars$pi
  list_pi0 <- priors$list_pi0
  a <- sum(pi) + list_pi0$a
  b <- sum(1 - pi) + list_pi0$b
  mean = a / (a + b)
  list(a = a, b = b, mean = mean)
}

