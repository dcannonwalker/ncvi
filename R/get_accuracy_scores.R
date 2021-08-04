get_accuracy_scores <- function(vifit, jagsfit, parnames, n) {
  P <- vifit$data$P
  G <- vifit$data$G
  post <- posterior(jagsfit)
  out <- list()
  if ("beta" %in% parnames) {
    post_beta <- post %>% select(starts_with("beta"))
    d <- list()
    for(j in 1:P) {
      for(i in 1:G) {
        colname <- paste0("beta", i, j)
        colnames(post)[(j - 1) * G + i] <- colname
        d[[(j - 1) * G + i]] <- post_betaij(vifit, i, j)
      }
    }

    vipostbeta <- list2DF(d)
    for(j in 1:P) {
      for(i in 1:G) {
        colnames(vipostbeta)[(j - 1) * G + i] <- paste0("beta", i, (j - 1))
      }
    }


    # get accuracy scores for all beta
    accuracy_scores_beta <- numeric()
    for (i in 1:(G * P)) {
      mean <- unlist(vipostbeta[1, i])
      sd <- unlist(vipostbeta[2, i])
      jags_samples <- unlist(post[, i])
      accuracy_scores_beta[i] <-
        get_accuracy_score(ncvi_list = list(density_fxn = dnorm,
                                            pars = list(mean = mean,
                                                        sd = sd)),
                           jags_samples = jags_samples,
                           n = n)
    }
    accuracy_scores_betaM <- matrix(accuracy_scores_beta, nrow = G, ncol = P)

    out$beta <- accuracy_scores_betaM
  }





  #hist(accuracy_scores, breaks = 100)

  # get accuracy scores for mu
  if ("mu0" %in% parnames) {
    postmu <- post %>% select(starts_with("mu_beta"))
    vipostmu <- matrix(nrow = 2, ncol = P)
    vipostmu[1, ] <- vifit$pars$theta$M[1:P]
    vipostmu[2, ] <- sqrt(diag(vifit$pars$theta$R))
    accuracy_scores_mu <- numeric()
    for (i in 1:P) {
      mean <- vipostmu[1, i]
      sd <- vipostmu[2, i]
      jags_samples <- unlist(postmu[, i])
      accuracy_scores_mu[i] <- get_accuracy_score(ncvi_list =
                                                    list(density_fxn = dnorm,
                                                         pars = list(mean = mean,
                                                                     sd = sd)),
                                                  jags_samples = jags_samples,
                                                  n = n)
    }
    out$mu <- accuracy_scores_mu
  }

  if ('precision_beta' %in% parnames) {
    postprecisionbeta <- post %>% select(starts_with("precision_beta"))
    a <- vifit$pars$theta$list_beta$a
    b <- vifit$pars$theta$list_beta$b
    accuracy_scores_precision_beta <- numeric()
    for (i in 1:P) {
      shape <- a[i]
      rate <- b[i]
      jags_samples <- unlist(postprecisionbeta[, i])
      accuracy_scores_precision_beta[i] <-
        get_accuracy_score(ncvi_list = list(density_fxn = dgamma,
                                            pars = list(shape = shape,
                                                        rate = rate)),
                           jags_samples = jags_samples,
                           n = n)
    }
    out$precision_beta <- accuracy_scores_precision_beta
  }

  out

}
