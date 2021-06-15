#' Simulate a set of data to fit with
#' the REH model using `fit_ncvi()` or `stan()`
#'
#' @export
sim_data_REH <- function(settings) {
  if (is.null(settings$seed)) set.seed(sample(1:10000, size = 1))
  else set.seed(settings$seed)

  P <- settings$P # number fixed effects / treatments
  N <- settings$N # number of samples / observations per group or gene
  U <- settings$U # number of random effects
  G <- settings$G # number of groups or genes
  X <- cbind(1, kronecker(diag(1, P),rep(1, N / P))[, 2:P]) # f.e. design
  Z <- kronecker(diag(U), rep(1, times = N/U)) # r.e. design

  ## known parameters
  precision_mu0 <- settings$precision_mu0 # prior precision for mu0
  precision_beta <- settings$precision_beta # prior precision for beta
  precision_u <- settings$precision_u # prior precision for u

  ## simulate
  mu0 <- mvtnorm::rmvnorm(n = 1,
                          mean = rep(0, P),
                          sigma = diag(1/precision_mu0, P))
  beta <- mvtnorm::rmvnorm(n = G,
                           mean = mu0,
                           sigma = diag(1/precision_beta, P))
  u <- mvtnorm::rmvnorm(n = G,
                        mean = rep(0, U),
                        sigma = diag(1/precision_u, U))

  truepars_REH <- list(mu0 = mu0,
                       beta = beta,
                       u = u,
                       precision_mu0 = precision_mu0,
                       precision_beta = precision_beta,
                       precision_u = precision_u)

  C <- cbind(X, Z)
  phi <- cbind(beta, u)
  eta <- c(apply(phi, 1, function(p) C %*% p))
  y_vector = rpois(N * G, lambda = exp(eta))
  y <- vector2list(y = y_vector,
                   G = G,
                   N = N)

  data_REH <- list(y = y,
                   C = C,
                   G = G,
                   P = P,
                   U = U,
                   y_vector = y_vector)



  init_phi <- list()

  for (i in seq(1, G)) {

    init_phi[[i]] <- list(
      mu = rep(3, (P + U)),
      Sigma = diag((P + U)))

  }

  # conflicting names of r, t, and s is very annoying --
  ## To do: fix variance component naming scheme
  ## To do: write functions that don't depend on particular
  ##   names for params (bigger issue)
  init_theta <- list(
    M = rep(0, (P + U)),
    R = diag(P),
    Tau = diag(c(rep(precision_beta, P), rep(precision_u, U))),
    precision_beta = precision_beta,
    precision_u = precision_u,
    precision_mu0 = precision_mu0
  )

  init_REH <- list(
    phi = init_phi,
    theta = init_theta
  )
  X_stan = kronecker(rep(1, G), X)
  Z_stan = kronecker(rep(1, G), Z)
  group <- rep(1:G, each = N)
  data_stan <- list(N = N,
                    G = G,
                    P = P,
                    U = U,
                    y = y_vector,
                    group = group,
                    X = X_stan,
                    Z = Z_stan,
                    sigB = sqrt(1/precision_beta),
                    sigu = sqrt(1/precision_u),
                    sig = sqrt(1/precision_mu0))

  list(data = data_REH, init = init_REH, data_stan = data_stan)
}
