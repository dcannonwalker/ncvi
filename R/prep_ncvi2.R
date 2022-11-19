#' Prepare list of objects required for `fit_ncvi()`. This is an
#' updated version of `prep_ncvi()`, which will replace it once
#' I've done some debugging.
#'
#' @param counts Data frame of counts with first column the
#' gene ID
#' @param S Vector of size factors with `length(S) = ncol(counts) - 1`
#' @param X Fixed effects design matrix (expects preparation in column 2
#' and treatment in column 3)
#' @param Z Random effects design matrix
#'
#' @export

prep_ncvi2 <- function(counts, S, X, Z) {

  N <- ncol(counts) - 1
  P <- ncol(X)
  U <- ncol(Z)
  y <- list()
  for (i in seq(1, nrow(counts))) {
    y[[i]] <- as.numeric(counts[i, 2:ncol(counts)])
  }
  data <- list(y = y, G = length(y),
               C = cbind(X, Z), P = P, U = U)

  preparation <- X[, 2]
  treatment <- X[, 3]
  sample <- factor(rep(1:(N / 2), 2))
  design <- model.matrix(~treatment + preparation + preparation:treatment +
                           sample)
  design <- design[, c(1:(3 + N / 2 - 2), ncol(design))]
  group <- factor(c(1 + treatment[1:(length(treatment) / 2)],
                    3 + treatment[1:(length(treatment) / 2)]))

  # Incorporate or calculate normalization factors
  if(missing(S)) {
    edger_list <- edgeR::DGEList(counts = counts[, 2:(N + 1)],
                 group = group)
    edger_list <- edgeR::calcNormFactors(edger_list)
    S <- edger_list$samples$norm.factors
  }
  else edger_list <- edgeR::DGEList(counts = counts[, 2:(N + 1)],
                    group = group,
                    norm.factors = S)
  data$S <- -log(S)

  # Estimate coefficients
  edger_list <-
    edgeR::estimateDisp(edger_list, design = design)
  edger_fit <-
    edgeR::glmQLFit(edger_list, design = design)
  coefs <- edger_fit$coefficients
  beta <- coefs[, c('(Intercept)', 'treatment',
                    'preparation', 'treatment:preparation')]
  beta[, '(Intercept)'] <- beta[, '(Intercept)'] + mean(edger_fit$offset)
  u <- coefs[, 4:(ncol(coefs) - 1)]

  beta_prec <- 1 / apply(beta, 2, var)
  u_prec <- 1 / var(c(u))
  beta_mean <- apply(beta, 2, mean)

  # * generate inits ----

  # * * init theta ----
  M <- c(beta_mean, rep(0, U))
  R <- diag(1, P)
  Tau <- diag(c(beta_prec, rep(u_prec, U)))
  a_beta <- beta_prec * 10
  b_beta <- rep(10, length(beta_prec))
  list_beta <- list(mean = beta_prec, a = a_beta,
                    b = b_beta)
  a_u <- 10 * u_prec
  b_u <- 10
  list_u <- list(mean = u_prec, a = a_u, b = b_u)
  theta <- list(M = M,
                R = R,
                Tau = Tau,
                precision_beta = beta_prec,
                precision_u = u_prec,
                pi0 = 0.8,
                list_beta = list_beta,
                list_u = list_u,
                precision_mu0 = 1)

  # * * init phi ----
  phi <- list()

  for (i in seq(1, length(data$y))) {
    phi[[i]] <- list(
      mu = c(beta[i, ], 0, u[i, ], 0),
      Sigma = theta$Tau
    )
  }

  # * * init pi ----
  pi <- rep(0.8, length(data$y))

  # * * collect inits ----
  init <- list(
    theta = theta,
    phi = phi,
    pi = pi
  )
  # * priors ----
  priors <- list(a_beta = a_beta, b_beta = b_beta,
                 a_u = a_u, b_u = b_u)


  list(
    data = data,
    init = init,
    priors = priors
  )
}
