#' Simulate a set of data to fit with
#' the REH mixture prior model using `fit_ncvi()` or `stan()`
#'
#' @export
sim_data_mixture <- function(settings) {
  if (is.null(settings$seed)) {
    seed <- sample(1:10000, size = 1)
    set.seed(seed)
  }
  else {
    seed <- settings$seed
    set.seed(seed)
  }

  P <- settings$P # number fixed effects / treatments
  N <- settings$N # number of samples / observations per group or gene
  U <- settings$U # number of random effects
  G <- settings$G # number of groups or genes
  X <- settings$X # f.e. design
  Z <- settings$Z # r.e. design
  ## known parameters

  # set values if fixing these precision parameters
  if (is.null(settings$priors)) {
    precision_mu0 <- settings$precision_mu0 # prior precision for mu0
    precision_beta <- settings$precision_beta # prior precision for beta
    precision_u <- NULL

    # will expect NULL if U = 0
    if (U != 0) precision_u <- settings$precision_u
  }


  # set values if estimating these precision parameters
  else if (!is.null(settings$priors)) {
    priors <- settings$priors
    a_beta <- priors$a_beta
    b_beta <- priors$b_beta
    precision_mu0 <- settings$precision_mu0

    # if using a set value for precision_beta
    if (is.null(settings$precision_beta)) {
      precision_beta <- a_beta / b_beta
      if (!is.null(settings$varComps)) {
        if (settings$varComps == T) {
          precision_beta = rep(a_beta / b_beta, P)
        }
      }
    }

    else if (!is.null(settings$precision_beta)) {
      precision_beta = settings$precision_beta
    }

    if (!is.null(settings$sim_precision)) {

      if(settings$sim_precision == "model") {
        precision_beta <- rgamma(P, shape = a_beta, rate = b_beta)
      }

      else if (settings$sim_precision == "flat") {
        precision_beta <- sim_precision_beta(a_beta, b_beta, type = "flat")
      }
    }


    list_beta <- list(mean = precision_beta,
                      a = a_beta,
                      b = b_beta)
    precision_u <- NULL

    if (U != 0) {
      a_u <- priors$a_u
      b_u <- priors$b_u
      if (!is.null(settings$sim_precision)) {
        if (settings$sim_precision == "model") {
          precision_u <- rgamma(1, shape = a_u, rate = b_u)
        }
      }

      else precision_u <- a_u / b_u

      list_u <- list(mean = precision_u,
                     a = a_u,
                     b = b_u)
    }
  }

  # set value if pi0 is fixed / given
  if (is.null(settings$priors$pi0)) pi0 <- settings$pi0
  else pi0 <- settings$priors$pi0

  ## simulate mu0, unless mu0 is provided
  if (is.null(settings$mu0)) {
    mu0 <- mvtnorm::rmvnorm(n = 1,
                            mean = rep(0, P),
                            sigma = diag(1/precision_mu0, P))
  }

  else mu0 <- settings$mu0

  if (length(precision_beta) == 1) {
    beta <- mvtnorm::rmvnorm(n = G,
                             mean = mu0,
                             sigma = diag(1/precision_beta, P))
  }

  else if (length(precision_beta) == P) {
    beta <- mvtnorm::rmvnorm(n = G,
                             mean = mu0,
                             sigma = diag(1/precision_beta))
  }

  D <- rbinom(n = G, 1, pi0)

  beta[, P] <- (1 - D) * beta[, P]
  u <- NULL

  C <- X
  phi <- beta

  if (U != 0) {
    u <- mvtnorm::rmvnorm(n = G,
                          mean = rep(0, U),
                          sigma = diag(1/precision_u, U))
    C <- cbind(X, Z)
    phi <- cbind(beta, u)
  }

  eta <- c(apply(phi, 1, function(p) C %*% p))
  if (!is.null(settings$checkSmall)) {
    if (settings$checkSmall == T) {
      eta[eta < 0] <- 0
    }
  }
  if (is.null(settings$family)) {
    y_vector = rpois(N * G, lambda = exp(eta))
  }
  else if (settings$family == "poisson" | is.null(settings$family)) {
    y_vector = rpois(N * G, lambda = exp(eta))
  }
  else if (settings$family == "nbinom") {
    if (!is.null(settings$nbinom$dispersion)) {
      dispersion <- settings$nbinom$dispersion
      y_vector = rnbinom(N * G, mu = exp(eta),
                         size = 1 / dispersion)
    }
    else if (settings$nbinom$est) {
      dispersion <- sapply(beta[, 1], function(b) dispersion_est(b))
      y_vector = rnbinom(N * G, mu = exp(eta),
                         size = 1 / rep(dispersion, each = N))
    }
    else if (!is.null(settings$nbinom$priors)) {
      shape <- settings$nbinom$priors$shape
      scale <- settings$nbinom$priors$scale
      dispersion <- rep(rgamma(G, shape = shape, scale = scale), each = N)
      y_vector = rnbinom(N * G, mu = exp(eta),
                         size = 1 / dispersion)

    }

  }
  if (settings$adjust$decide == T) {
    y_vector <- y_vector + settings$adjust$offset
  }
  y <- vector2list(y = y_vector,
                   G = G,
                   N = N)



  if (length(precision_beta) == 1) {
    true_theta = list(
      M = c(mu0, rep(0, U)),
      R = diag(1 / precision_mu0, P),
      Tau = diag(c(rep(precision_beta, P), rep(precision_u, U))),
      precision_beta = precision_beta,
      precision_u = precision_u,
      precision_mu0 = precision_mu0,
      pi0 = pi0
    )
  }

  else if (length(precision_beta) == P) {
    true_theta = list(
      M = c(mu0, rep(0, U)),
      R = diag(1 / precision_mu0, P),
      Tau = diag(c(precision_beta, rep(precision_u, U))),
      precision_beta = precision_beta,
      precision_u = precision_u,
      precision_mu0 = precision_mu0,
      pi0 = pi0
    )
  }

  if(settings$adjust$decide == T) {
    true_theta$M[1] = true_theta$M[1] + log(settings$adjust$offset)
  }

  if (!is.null(priors)) {
    true_theta$list_beta <- list_beta
    if (U != 0) true_theta$list_u <- list_u
  }

  temp_Sigma <- MASS::ginv(true_theta$Tau)
  true_phi <- list()
  for (i in seq(1, G)) {
    true_phi[[i]] <- list(mu = c(beta[i, ], u[i, ]),
                     Sigma = temp_Sigma)
  }

  truepars_init <- list(theta = true_theta,
                        phi = true_phi,
                        pi = D)
  truepars_ncvi <- list(mu0 = mu0,
                        beta = beta,
                        u = u,
                        precision_mu0 = precision_mu0,
                        precision_beta = precision_beta,
                        precision_u = precision_u,
                        D = D,
                        pi0 = pi0)
  if (!is.null(settings$family)) {
    if (settings$family == "nbinom") {
      truepars_ncvi$dispersion <- dispersion
    }
  }


  data_ncvi <- list(y = y,
                   C = C,
                   G = G,
                   P = P,
                   U = U,
                   min_pi = settings$min_pi)



  init_phi <- list()
  init_pi <- list()

  for (i in seq(1, G)) {

    init_phi[[i]] <- list(
      mu = rep(0, (P + U)),
      Sigma = diag((P + U)))

    init_pi[[i]] <- pi0

  }

  init_theta <- list(
    M = rep(0, (P + U)),
    R = diag(P),
    Tau = diag(c(rep(precision_beta, P), rep(precision_u, U))),
    precision_beta = precision_beta,
    precision_u = precision_u,
    precision_mu0 = precision_mu0,
    pi0 = pi0
  )
  if (!is.null(priors)) {
    init_theta$list_beta <- list_beta
    if (U != 0) init_theta$list_u <- list_u
  }

  init_ncvi <- list(
    phi = init_phi,
    theta = init_theta,
    pi = init_pi
  )

  group <- rep(1:G, each = N)





  preparation <- rep(c(0, 1), each = N / 2)
  treatment <- rep(rep(c(0, 1), each = N / 4), 2)
  sample <- factor(rep(1:(N / 2), 2))
  design <- model.matrix(~treatment + preparation + preparation:treatment +
                           sample)
  design <- design[, c(1:(3 + N / 2 - 2), ncol(design))]
  y_df <- y |> list2DF() |> data.table::transpose()
  other_data <- list(counts = y_df,
                     group = treatment,
                     design = design)

  list(data = data_ncvi, init = init_ncvi,
       other_data = other_data,
       seed = seed,
       truepars = truepars_ncvi,
       true_init = truepars_init,
       D = D,
       priors = priors,
       adjust = settings$adjust)
}
