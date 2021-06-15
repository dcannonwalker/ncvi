## code to prepare `fit_compare2stan` dataset goes here
settings <- list(
  P = 2, # number fixed effects / treatments
  N = 12, # number of samples / observations per group or gene
  U = 6, # number of random effects
  G = 10, # number of groups or genes
  X = cbind(1, kronecker(diag(1, 2),rep(1, 12 / 2))[, 2:2]), # f.e. design
  Z = kronecker(diag(6), rep(1, times = 12 / 6)), # r.e. design

  ## known parameters
  precision_mu0 = 1, # prior precision for mu0
  precision_beta = 1/2, # prior precision for beta
  precision_u = 1/3 # prior precision for u
)
sim1 <- sim_data_REH(settings)

fit_vb <- fit_ncvi(data = sim1$data, init = sim1$init, update_pars = update_pars_REH,
         elbo = elbo_REH, differentials = differentials_new)

str(sim1)

fit_stan <- rstan::stan("stan/REH.stan",
     data = sim1$data_stan,
     chains = 4,
     warmup = 1000,
     iter = 2000,
     cores = 2,
     refresh = 0)

fit_compare2stan <- list(fit_vb = fit_vb, fit_stan = fit_stan)

print(fit_stan, pars = c("beta", "mu"))
usethis::use_data(fit_compare2stan, overwrite = TRUE)
