get_beta3_post <- function(fit) {
  w_post <- get_w_post(fit)
  pi <- fit$pars$pi
  mean <- as.numeric(w_post$mu_list) * (1 - pi)

  list(mean = mean, w_post = w_post, pi = pi)
}
