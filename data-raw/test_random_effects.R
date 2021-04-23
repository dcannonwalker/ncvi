## code to prepare `test_random_effects` dataset goes here
## constructing random effects model
set.seed(112)

P <- 2 # number fixed effects / treatments
N <- 12 # number of samples / observations per group or gene
U <- 6 # number of random effects
G <- 10 # number of groups or genes
X <- cbind(1, kronecker(diag(1, P),rep(1, N / P))[, 2:P]) # f.e. design
Z <- kronecker(diag(U), rep(1, times = N/U)) # r.e. design

## known parameters
precision_mu0 <- 1 # prior precision for mu0
precision_beta <- 1/2 # prior precision for beta
precision_u <- 1/3 # prior precision for u

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

# priors_REH <- list(
#   a = ar,
#   b = br,
#   at = at,
#   bt = bt
# )

usethis::use_data(data_REH, init_REH, truepars_REH, overwrite = TRUE)
