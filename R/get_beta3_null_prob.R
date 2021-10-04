get_w_post <- function(fit) {
  mu_list <- list()
  sig2_list <- list()
  for (i in seq(1, fit$data$G)) {
    mu_list[[i]] <- fit$pars$phi[[i]]$mu[4]
    sig2_list[[i]] <- fit$pars$phi[[i]]$Sigma[4, 4]
  }

  list(mu_list = mu_list, sig2_list = sig2_list)
}

get_beta3_null_prob <- function(fit, cutoff) {
  w_post <- get_w_post(fit)
  null_prob <- c()
  for (i in seq(1, length(w_post$mu_list))) {
    mu <- w_post$mu_list[[i]]
    sig2 <- w_post$sig2_list[[i]]
    p <- pnorm(cutoff, mean = mu, sd = sqrt(sig2), lower.tail = F) +
      pnorm(-cutoff, mean = mu, sd = sqrt(sig2))
    null_prob[i] <- 1 - p
  }

  (1 - fit$pars$pi) * null_prob + fit$pars$pi
}

