
#' New update function for phi param in hierarchical mean poisson example
#'
#' @param data List of observed variables, with `y, C, G`
#' @param pars List of pars for update, with `phi, theta`
#' @param differentials List of differentials functions, with `mu, Sigma`
#' @export
update_phi_new <- function(data, pars, differentials) {

  phi <- pars$phi
  theta <- pars$theta
  y <- data$y
  C <- data$C
  G <- data$G

  for (i in seq(1, G)) {

    phi[[i]] <- nc_update_mvn(data = list(y = y[[i]], C = C),
                              pars = list(phi = phi[[i]], theta = theta),
                              differentials = differentials)

  }

  # return phi
  phi

}

#' New update function for theta param in hierarchical example;
#' treats everything but `theta$M` as fixed
#'
#' @param data List of observed variables, with `y, C, G`
#' @param pars List of pars for update, with `phi, theta`
#' @param priors
#' @export
update_theta_new <- function(data, pars, priors) {

  theta <- pars$theta

  pars_mu0 <- conjugate_update_mvn_hierarchical_mean(data, pars)

  theta$M <- pars_mu0$M

  theta$R <- pars_mu0$R

  # return theta
  theta

}

#' New wrapper update wrapper for hierarchical poisson / ex4
#'
#' @param data
#' @param init
#' @param args
#' @export
update_pars_new <- function(data, pars, args) {
  pars$phi <- update_phi_new(data, pars, args$differentials)
  pars$theta <- update_theta_new(data, pars, args$priors)

  # return pars
  pars
}


