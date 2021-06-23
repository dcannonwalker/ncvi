## code to prepare `exMixtureDataFixedOnly` dataset goes here

set.seed(5562)

P <- 2 # number fixed effects / treatments
N <- 10 # number of samples / observations per group or gene
G <- 200 # number of groups or genes
X <- cbind(1, kronecker(diag(1, P),rep(1, N / P))[, 2:P]) # f.e. design

## known parameters
precision_mu0 <- 1 # prior precision for mu0
precision_beta <- 2 # prior precision for beta
pi0 <- 0.2

# function to simulate beta with mixture prior
sim_beta <- function(G, pi, mean, sigma) {
  beta <- mvtnorm::rmvnorm(G, mean = mean, sigma = sigma)
  D <- rbinom(G, 1, pi)
  beta[,2] <- D * beta[,2]

  beta
}


mu0 <- c(1, 5)
beta <- sim_beta(G, pi0, mean = mu0, sigma = diag(1 / precision_beta, P))

eta <- c(apply(beta, 1, function(p) X %*% p))
y_vector = rpois(N * G, lambda = exp(eta))
y <- ncvi::vector2list(y = y_vector,
                 G = G,
                 N = N)

exMixtureDataFixedOnly <- list(y = y,
                      G = G,
                      C = X,
                      P = P,
                      U = 0,
                      y_vector = y_vector,
                      precision_beta = precision_beta,
                      precision_mu0 = precision_mu0,
                      pi0 = pi0,
                      beta = beta,
                      mu0 = mu0)

usethis::use_data(exMixtureDataFixedOnly, overwrite = TRUE)
