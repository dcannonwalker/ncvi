## code to prepare `test_hierarchical_poisson` dataset goes here

set.seed(558)
at <- 1
bt <- 2
ar <- 1
br <- 2
P <- 2
G <- 10
N <- 10

# t <- rgamma(1, shape = at, rate = bt)
#
# r <- rgamma(1, shape = ar, rate = br)

s <- 1

t <- 1/3

mu0 <- mvtnorm::rmvnorm(1, mean = rep(0, P), sigma = diag(1/s, P))

beta <- mvtnorm::rmvnorm(G, mean = mu0, sigma = diag(1/t, P))

X <- cbind(1, kronecker(diag(1, P),rep(1, N/P))[, 2:P])

H <- diag(P)

eta_mat <- apply(beta, 1, function(b) X%*%b)

eta <- c(apply(beta, 1, function(b) X%*%b))

y <- rpois(N*G, lambda = exp(eta))

eta
y

y_list <- list()

for (i in seq(1, G)) {

  y_list[[i]] <- y[((i - 1) * N + 1):(i * N)]

}

data_legacy <- list(y = y_list,
                etc = list(X = X,
                           t = t,
                           s = s,
                           G = G))

SigmaBinit <- list()
for (i in seq(1,G)) {
  SigmaBinit[[i]] <- diag(P)
}

EBinit <- list()
for (i in seq(1,G)) {
  EBinit[[i]] <- rep(0,P)
}

init_legacy <- list(phi = list(EB = EBinit,
                           SigmaB = SigmaBinit),
                theta = list(EM = rep(0,P)))

differentials_legacy <- list(SigmaB = differential_SigmaB,
                             EB = differential_EB)

data_new <- list(
  y = y_list,
  C = X,
  H = diag(P),
  G = G
)

init_phi <- list()

for (i in seq(1, G)) {

  init_phi[[i]] <- list(
    mu = rep(0, P),
    Sigma = diag(P))

}

# conflicting names of r, t, and s is very annoying --
## To do: fix variance component naming scheme
## To do: write functions that don't depend on particular
##   names for params (bigger issue)
init_theta <- list(
  M = rep(0, P),
  R = diag(P),
  Tau = diag(t, P),
  t = s,
  r = t
)

init_new <- list(
  phi = init_phi,
  theta = init_theta
)

priors_new <- list(
  a = ar,
  b = br,
  at = at,
  bt = bt
)

differentials_new <- list(
  Sigma = function(data, pars) {
    d_mvn_cov(data, pars, d_mvn_cov_poisson)
  },
  mu = function(data, pars) {
    d_mvn_mean(data, pars, d_mvn_mean_poisson)
  }
)


test_hierarchical_poisson <- list(data_legacy = data_legacy,
                                  init_legacy = init_legacy,
                                  differentials_legacy = differentials_legacy,
                                  data_new = data_new,
                                  init_new = init_new,
                                  priors_new = priors_new,
                                  differentials_new = differentials_new)

usethis::use_data(test_hierarchical_poisson, overwrite = TRUE)
