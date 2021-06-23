library(ncvi)
devtools::load_all()
settings <- list(
  P = 2, # number fixed effects / treatments
  N = 12, # number of samples / observations per group or gene
  U = 6, # number of random effects
  G = 100, # number of groups or genes
  X = cbind(1, kronecker(diag(1, 2),rep(1, 12 / 2))[, 2:2]), # f.e. design
  Z = kronecker(diag(6), rep(1, times = 12 / 6)), # r.e. design

  ## known parameters
  precision_mu0 = 1, # prior precision for mu0
  precision_beta = 1/2, # prior precision for beta
  precision_u = 1/3 # prior precision for u
)
sim1 <- sim_data_REH(settings = settings)
str(sim1)
m <- rstan::stan("stan/REH.stan", data = sim1$data_stan,
          chains = 4,
          warmup = 1000,
          iter = 2000,
          cores = 2,
          refresh = 0)
m <- rstan::stan_model("stan/REH.stan")

bench <- rbenchmark::benchmark("fit_ncvi" = {
  sim1 <- sim_data_REH(settings)
  fit_ncvi(data = sim1$data, init = sim1$init, update_pars = update_pars_REH,
           elbo = elbo_REH,
           differentials = differentials_new)
},
"fit_stanvb" = {
  sim1 <- sim_data_REH(settings)
  rstan::vb(m, data = sim1$data_stan)
},
replications = 10)
bench

vb <- rstan::vb(m, data = sim1$data_stan)
summary_vb <- summary(vb)
print(summary_vb, pars = "mu")
print(vb, pars = "mu")
ncvi <- fit_ncvi(data = sim1$data, init = sim1$init, update_pars = update_pars_REH,
                 elbo = elbo_REH,
                 differentials = differentials_new)
ncvi$pars$theta
