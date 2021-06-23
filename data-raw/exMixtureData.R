## code to prepare `exMixtureData` dataset goes here

# code to simulate data for model with
# mixture prior on one fixed effect parameter

set.seed(112)

P <- 2 # number fixed effects / treatments
N <- 10 # number of samples / observations per group or gene
U <- 2 # number of random effects
G <- 10 # number of groups or genes
X <- cbind(1, kronecker(diag(1, P),rep(1, N / P))[, 2:P]) # f.e. design
Z <- kronecker(diag(U), rep(1, times = N/U)) # r.e. design

## known parameters
precision_mu0 <- 1 # prior precision for mu0
precision_beta <- 1/2 # prior precision for beta
precision_u <- 5 # prior precision for u
pi0 <- 0.2

# function to simulate beta with mixture prior
sim_beta <- function(G, pi, mean, sigma) {
  beta <- mvtnorm::rmvnorm(G, mean = mean, sigma = sigma)
  D <- rbinom(G, 1, pi)
  beta[,2] <- D * beta[,2]

  beta
}


mu0 <- mvtnorm::rmvnorm(1, sigma = diag(1 / precision_mu0, P))
beta <- sim_beta(G, pi0, mean = mu0, sigma = diag(1 / precision_beta, P))
u <- mvtnorm::rmvnorm(n = G,
                      mean = rep(0, U),
                      sigma = diag(1 / precision_u, U))

C <- cbind(X, Z)
phi <- cbind(beta, u)
eta <- c(apply(phi, 1, function(p) C %*% p))
y_vector = rpois(N * G, lambda = exp(eta))
y <- vector2list(y = y_vector,
                 G = G,
                 N = N)

exMixtureData <- list(y = y,
                      G = G,
                      C = C,
                      P = P,
                      U = U,
                      y_vector = y_vector,
                      precision_beta = precision_beta,
                      precision_mu0 = precision_mu0,
                      precision_u = precision_u,
                      pi0 = pi0)


usethis::use_data(exMixtureData, overwrite = TRUE)
