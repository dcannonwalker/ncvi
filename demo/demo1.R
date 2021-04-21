library(devtools)
load_all()

data <- list()
init <- list()
elbo <- function(data, pars) 1
update_pars <- function(data, pars, args) {
  print(args)
  pars
}

fit_ncvi(data, init, elbo = elbo,
         update_pars = update_pars, test = "test")

set.seed(111111)
N <- 100
P <-  2
sig2 <- 3

beta <- mvtnorm::rmvnorm(1, sigma = diag(rep(sig2, P)))

X <- cbind(1, kronecker(diag(1, P),rep(1, N/P))[, 2:P])

H <- diag(P)

eta <- X %*% t(beta)

y <- rpois(N, lambda = exp(eta))

beta

data = list(y = y, C = X, H = H)

init = list(phi = list(mu = rep(5, P),
                       Sigma = diag(nrow = P)),
            theta = list(r = 1 / sig2,
                         Tau = H * 1 / sig2,
                         M = rep(0, P)))

differentials <- list(
  Sigma = function(data, pars) {
    d_mvn_cov(data, pars, d_mvn_cov_poisson)
  },
  mu = function(data, pars) {
    d_mvn_mean(data, pars, d_mvn_mean_poisson)
  }
)

update_phi <- nc_update_mvn

update_pars <- function(data, pars, args) {

  update_phi(data, pars, args$differentials)

}

elbo <- function(data, pars){
  y <- data$y
  C <- data$C
  mu <- pars$phi$mu
  Sigma <- pars$phi$Sigma
  A <- c(C %*% mu + 0.5 * diag(C %*% Sigma %*% t(C)))
  sig2 <- 1 / pars$theta$r

  # returns:
  y %*% C %*% mu - sum(exp(A)) -
    (sum(mu^2) + sum(diag(Sigma))) / (2 * sig2) +
    0.5 * log(det(Sigma))
}

update_pars <- function(data, pars, args) {

  phi <- update_phi(data, pars, args$differentials)

  list(phi = phi, theta = pars$theta)

}

options = list(max_iter = 10,
               elbo_delta = 0.0001,
               verbose = T)

fit_ncvi(data, init,
         update_pars = update_pars,
         elbo = elbo,
         options = options,
         differentials = differentials)
beta
