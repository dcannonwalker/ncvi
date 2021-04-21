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

t <- 1

r <- 1/3

mu0 <- mvtnorm::rmvnorm(1, mean = rep(0, P), sigma = diag(1/t, P))

beta <- mvtnorm::rmvnorm(G, mean = mu0, sigma = diag(1/r, P))

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

data <- list(
  y = y_list,
  C = X,
  H = H,
  G = G
)

init_phi <- list()
lapply(y_list, function(x) log(mean(x)))
for (i in seq(1, G)) {

  init_phi[[i]] <- list(
    mu = rep(0, P),
    Sigma = diag(P))

}
init_phi
init_theta <- list(
  M = rep(0, P),
  R = diag(P),
  Tau = diag(r, P),
  t = t,
  r = r
)

init <- list(
  phi = init_phi,
  theta = init_theta
)

priors <- list(
  a = ar,
  b = br,
  at = at,
  bt = bt
)

update_phi_pkg <- function(data, pars, differentials) {

  phi <- pars$phi
  theta <- pars$theta
  y <- data$y
  C <- data$C
  G <- data$G

  for (i in seq(1, G)) {

    phi[[i]] <- nc_update_mvn(data = list(y = y[[i]], C = C),
                              pars = list(phi = phi[[i]], theta = theta),
                              differentials = differentials)

  }

  # return phi
  phi

}

update_theta_pkg <- function(data, pars, priors) {

  theta <- pars$theta

  pars_mu0 <- conjugate_update_mvn_hierarchical_mean(pars)

  theta$M <- pars_mu0$M

  theta$R <- pars_mu0$R

  # return theta
  theta

}

differentials <- list(
  Sigma = function(data, pars) {
    d_mvn_cov(data, pars, d_mvn_cov_poisson)
  },
  mu = function(data, pars) {
    d_mvn_mean(data, pars, d_mvn_mean_poisson)
  }
)

options = list(max_iter = 1000,
               elbo_delta = 0.1,
               verbose = T)
