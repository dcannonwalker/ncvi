library(devtools)
load_all()

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

data_vb <- list(y = y_list,
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

init_vb <- list(phi = list(EB = EBinit,
                           SigmaB = SigmaBinit),
                theta = list(EM = rep(0,P)))

differentials <- list(SigmaB = differential_SigmaB,
                      EB = differential_EB)

update_pars <- function(data, pars, args) {
  pars$phi <- update_phi_ex4(data, pars, args$differentials)
  pars$theta <- update_theta_ex4(data, pars)
  pars
}

fit <- fit_ncvi(data = data_vb,
         init = init_vb,
         update_pars = update_pars,
         elbo = elbo_ex4,
         differentials = differentials)

fit
beta
