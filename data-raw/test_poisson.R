## code to prepare `DATASET` dataset goes here

N <- 10
P <- 2
C <- cbind(1, kronecker(diag(1, P),rep(1, N / P))[, 2:P])
M <- c(1, 1)
t <- 1 / 3
Tau <- diag(t, P)
mu <- c(0.2, 1.5)
Sigma <- diag(P)
y <- rpois(N, 3)

data_legacy <- list(y = y, etc = list(X = C, t = t))
data_new <- list(y = y, C = C)
pars_legacy <- list(phi = list(EB = mu, SigmaB = Sigma),
                    theta = list(EM = M))
pars_new <- list(phi = list(mu = mu, Sigma = Sigma),
                 theta = list(Tau = Tau, M = M))

differentials_legacy <- list(SigmaB = differential_SigmaB,
                             EB = differential_EB)

differentials_new <- list(
  Sigma = function(data, pars) {
    d_mvn_cov(data, pars, d_mvn_cov_poisson)
  },
  mu = function(data, pars) {
    d_mvn_mean(data, pars, d_mvn_mean_poisson)
  }
)

test_poisson <- list(data_legacy = data_legacy,
                     data_new = data_new,
                     pars_legacy = pars_legacy,
                     pars_new = pars_new,
                     differentials_legacy = differentials_legacy,
                     differentials_new = differentials_new)
usethis::use_data(test_poisson, overwrite = TRUE)


