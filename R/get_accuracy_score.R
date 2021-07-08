get_accuracy_score <- function(ncvi_list, jags_samples, n) {
  jdens <- density(jags_samples, n = n)
  ncvi_list$pars$x <- jdens$x
  nvec <- do.call(ncvi_list$density_fxn, ncvi_list$pars)
  100 * (1 - mean(abs(nvec - jdens$y)) / 2)
}
