## NCVMP algorithm for 2-group hierarchical model
## y_i | beta_i is Poisson(exp(XB_i))
## beta_i | mu, T is MVN(mu, (1/T)I)
## mu (or M) | sig2 is MVN(0, sig2I)
## t (or T, or 1/sig2B) can be known or can be Gamma(a,b), with a and b known
## 1/sig2 (or S, or s) can be known or can be Gamma(c,d) with c and d known

## pars = list(theta, phi)
## pars$phi = list(EB, SigmaB)
## EB = vector(EB1,...,EBG)
## SigmaB = list(SigmaB1,...,SigmaBG)
## pars$theta = list(EM, VM, ET, ES)


## Non-conjugate updates for q(beta_i)

# update the variational covariance matrix
# for beta_i
# data and pars should be for a particular i, e.g.
# data = list(y = yi, etc)
# pars = list(phi = phii, theta)

#' Legacy update for mvn covariance
#'
#' @inheritParams differential_SigmaB
#' @export
update_SigmaB <- function(data, pars, differential_SigmaB){
  return(
    MASS::ginv(
      -2*differential_SigmaB(data, pars)
    )
  )
}

#' Legacy update for mvn mean
#'
#' @inheritParams differential_SigmaB
#' @export
update_EB <- function(data, pars, differential_EB){
  return(
    c(pars$phi$EB + pars$phi$SigmaB%*%differential_EB(data, pars))
  )
}

#' Legacy differential for mvn covariance
#'
#' @param data List with `etc$X, y, etc$t`
#' @param pars List with `phi$EB, phi$SigmaB`
#' @export
differential_SigmaB <- function(data, pars){
  X <- data$etc$X
  A <- c(X%*%pars$phi$EB + 0.5*diag(X%*%pars$phi$SigmaB%*%t(X)))
  return(
    -0.5*
      (t(X)%*%diag(exp(A))%*%X +
         diag(data$etc$t, nrow = ncol(X)))
  )
}

#' Legacy differential for mvn mean
#'
#' @inheritParams differential_SigmaB
#' @export
differential_EB <- function(data, pars){
  X <- data$etc$X
  EB <- pars$phi$EB
  A <- c(X%*%EB + 0.5*diag(X%*%pars$phi$SigmaB%*%t(X)))
  return(
    t(X)%*%(data$y - exp(A)) -
      (EB - pars$theta$EM)*data$etc$t
  )
}

## conjugate updates for other parameters

# update the variational mean for the hierarchical
# mean parameter
# known sig2 and sig2B

#' Legacy update for conjugate mvn hierarchical mean
#'
#' @param data List with `etc$t, etc$G, etc$s`
#' @param pars List with `phi$EB`
#' @export
update_EM <- function(data, pars){
  EB <- data.frame(pars$phi$EB)
  return(
    rowSums(EB)/(data$etc$s/data$etc$t + data$etc$G)
  )
}

# for unknown sig2 and sig2B
# replace s and t with ES and ET
# update_EM <- function(data, pars){
#   return(
#     sum(EB)/(ES/ET+G)
#   )
# }

# update the variance for the hierarchical mean parameter
# this would be fixed for known sig2 and sig2B

update_VM <- function(data, pars){
  return(
    (pars$theta$ES + data$etc$G*pars$theta$ET)^-1
  )
}

# update the variational mean for T

# update the variational mean for 1/sig2

## ELBO
# for known sig2, sig2B
# convenient use of elbo function for non-hierarchical version of the model

#' ELBO contribution for a single group
#'
#' @inheritParams differential_SigmaB
#' @export
elbo_i_legacy <- function(data, pars){
  y <- data$y
  X <- data$X
  EB <- pars$EB
  SigmaB <- pars$SigmaB
  A <- c(X%*%EB + 0.5*diag(X%*%SigmaB%*%t(X)))
  t <- data$t
  return(
    y%*%X%*%EB - sum(exp(A)) -
      t*(sum(EB^2) + sum(diag(SigmaB)))/2 +
      0.5*log(det(SigmaB))
  )
}

#' Extra ELBO component outside group contribution
#'
#' @param data The entire `data` arg
#' @param pars The entire `pars` arg
#' @export
elbo_extra_legacy <- function(data, pars){
  s <- data$etc$s
  t <- data$etc$t
  EM <- pars$theta$EM
  EB <- data.frame(pars$phi$EB)
  G <- data$etc$G
  return(
    -0.5*s*sum(EM^2) + t*rowSums(EB)%*%EM -
      0.5*G*t*sum(EM^2)
  )
}

#' Entire ELBO function for ex4
#'
#' @inheritParams elbo_extra_legacy
#' @export
elbo_legacy <- function(data, pars){
  L <- elbo_extra_legacy(data, pars)
  SigmaB <- pars$phi$SigmaB
  EB <- pars$phi$EB
  theta <- pars$theta
  X <- data$etc$X
  t <- data$etc$t
  y <- data$y
  G <- data$etc$G
  for (i in seq(1,G)) {
    L <- L + elbo_i_legacy(data = list(y = y[[i]], X = X, t = t),
                   pars = list(EB = EB[[i]], SigmaB = SigmaB[[i]]))
  }
  return(L)
}
## wrapper functions

# function to update the non-conjugate distribution
# of the beta_i parameters

#' Update using the legacy format/differentials for ex4
#'
#' @export
update_phi_legacy <- function(data, pars, differentials){
  SigmaB <- pars$phi$SigmaB
  EB <- pars$phi$EB
  theta <- pars$theta
  etc <- data$etc
  y <- data$y
  G <- data$etc$G
  for (i in seq(1,G)) {
    SigmaB[[i]] <- update_SigmaB(data = list(y = y[[i]],
                                             etc = etc),
                                 pars = list(phi = list(EB = EB[[i]],
                                                        SigmaB = SigmaB[[i]]),
                                             theta = theta),
                                 differential_SigmaB = differentials$SigmaB)
    EB[[i]] <- update_EB(data = list(y = y[[i]],
                                     etc = etc),
                         pars = list(phi = list(EB = EB[[i]],
                                                SigmaB = SigmaB[[i]]),
                                     theta = theta),
                         differential_EB = differentials$EB)
  }
  return(list(EB = EB, SigmaB = SigmaB))
}

# function to update the conjugate distributions
# of mu, sig2beta, sig2...

# known sig2, sig2B

#' Legacy format update for `theta` in ex4
#'
#' @export
update_theta_legacy <- function(data, pars){
  theta <- pars$theta
  theta$EM <- update_EM(data, pars)
  return(theta)
}

#' Wrapper for legacy format updates to pass to `fit_ncvi()`
#'
#' @export
update_pars_legacy <- function(data, pars, args) {
  pars$phi <- update_phi_legacy(data, pars, args$differentials)
  pars$theta <- update_theta_legacy(data, pars)
  pars
}
# unknown sig2, sig2B
# update_theta <- function(data, pars){
#   theta <- pars$theta
#   theta$EM <- update_EM(data, pars)
#   theta$VM <- update_VM(data, pars)
#   theta$ET <- update_ET(data, pars)
#   theta$ES <- update_ES(data, pars)
#   return(theta)
# }



