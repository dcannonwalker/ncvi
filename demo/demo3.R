library(devtools)
load_all()

data("test_hierarchical_poisson")

update_pars_legacy <- function(data, pars, args) {
  pars$phi <- update_phi_ex4(data, pars, args$differentials)
  pars$theta <- update_theta_ex4(data, pars)
  pars
}

data_legacy

fit_legacy <- fit_ncvi(data = data_legacy,
         init = init_legacy,
         update_pars = update_pars_legacy,
         elbo = elbo_ex4,
         differentials = differentials_legacy)

# update_pars_new <- function(data, pars, args) {
#   pars$phi <- update_phi
# }

