
devtools::load_all()
data("data_REH")
data("init_REH")
init_REH
conjugate_update_mvn_REH

conjugate_update_mvn_REH(data = data_REH, pars = init_REH)
update_theta_REH
fit <- fit_ncvi(data = data_REH,
         init = init_REH,
         update_pars = update_pars_REH,
         elbo = elbo_REH,
         differentials = differentials_new,
         priors = NULL
         )
data("truepars_REH")
truepars_REH$beta
t(sapply(fit$phi, function(x) x$mu[1:2]))

truepars_REH$u
t(sapply(fit$phi, function(x) x$mu[3:8]))


