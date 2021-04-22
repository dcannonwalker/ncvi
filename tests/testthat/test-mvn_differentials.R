test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

N <- 10
P <- 2
C <- cbind(1, kronecker(diag(1, P),rep(1, N / P))[, 2:P])
M <- c(1, 1)
t <- 1 / 3
Tau <- diag(t, P)
mu <- c(0.2, 1.5)
Sigma <- diag(P)
y <- rpois(N, exp(C%*%mu))

data1 <- list(y = y, etc = list(X = C, t = t))
data2 <- list(y = y, C = C)
pars1 <- list(phi = list(EB = mu, SigmaB = Sigma),
              theta = list(EM = M))
pars2 <- list(phi = list(mu = mu, Sigma = Sigma),
              theta = list(Tau = Tau, M = M))
test_that("different mean differentials functions agree", {
  expect_equal(
    differential_EB(data1, pars1),
    d_mvn_mean(data2, pars2, d_mvn_mean_poisson))
})

test_that("different cov differentials functions agree", {
  expect_equal(
    d_mvn_cov(data2, pars2, d_mvn_cov_poisson),
    differential_SigmaB(data1, pars1)
  )
})
