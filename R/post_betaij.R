post_betaij <- function(fit, i, j) {
  list(mean = fit$pars$phi[[i]]$mu[j],
       sd = sqrt(fit$pars$phi[[i]]$Sigma[j, j]))
}
