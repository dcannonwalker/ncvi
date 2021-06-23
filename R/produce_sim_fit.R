produce_sim_fit <- function(settings, generate_sim,
                            generate_init,
                            differentials,
                            update_pars,
                            elbo,
                            fit_options) {
  sim <- generate_sim(settings)
  init <- generate_init(settings)
  fit <- fit_ncvi(data = sim, init = init, update_pars = update_pars,
                  elbo = elbo,
                  options = fit_options,
                  differentials = differentials)

  fit
}
