library(devtools)
load_all()

data("test_hierarchical_poisson")


fit_legacy <- fit_ncvi(data = test_hierarchical_poisson$data_legacy,
         init = test_hierarchical_poisson$init_legacy,
         update_pars = update_pars_legacy,
         elbo = elbo_legacy,
         differentials = test_hierarchical_poisson$differentials_legacy,
         options = list(max_iter = 100,
                        elbo_delta = 0.0001,
                        verbose = T))


fit_new <- fit_ncvi(data = test_hierarchical_poisson$data_new,
                    init = test_hierarchical_poisson$init_new,
                    update_pars = update_pars_new,
                    elbo = elbo_hierarchical,
                    differentials = test_hierarchical_poisson$differentials_new,
                    options = list(max_iter = 100,
                                   elbo_delta = 0.0001,
                                   verbose = T))

lapply(fit_new$phi, function(x) x$mu)
fit_legacy$phi$EB

