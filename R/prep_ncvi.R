#' Prepare list of objects required for `fit_ncvi()`
#'
#' @param counts Data frame of counts with first column the
#' gene ID
#' @param X Fixed effects design matrix
#' @param Z Random effects design matrix
#'
#' @export
prep_ncvi <- function(counts, X, Z) {
    P <- ncol(X)
    U <- ncol(Z)
    y <- list()
    for (i in seq(1, nrow(counts))) {
        y[[i]] <- as.numeric(counts[i, 2:ncol(counts)])
    }
    data <- list(y = y, G = length(y),
                      C = cbind(X, Z), P = P, U = U)

    # * generate inits ----
    beta <- list()
    for (i in seq(1, length(y))) {
        beta[[i]] <- glm(y[[i]] ~ 0 + X ,family = poisson)$coefficients
    }



    beta.df <- list2DF(beta) |> data.table::transpose()
    beta.sd <- beta.df |> apply(2, sd)
    beta.prec <- beta.sd^(-2)
    beta.mean <- beta.df |> apply(2, mean)

    # * * init theta ----
    M <- c(beta.mean, rep(0, U))
    R <- diag(1, P)
    u.prec = 1 # needs improvement
    Tau <- diag(c(beta.prec, rep(u.prec, U)))
    a.beta <- beta.prec * 10
    b.beta <- 10
    list_beta <- list(mean = beta.prec, a = a.beta,
                           b = b.beta)
    list_u <- list(mean = u.prec, a = 2, b = 2) # needs improvement
    theta <- list(M = M,
                       R = R,
                       Tau = Tau,
                       precision_beta = beta.prec,
                       precision_u = u.prec,
                       pi0 = 0.8,
                       list_beta = list_beta,
                       list_u = list_u,
                       precision_mu0 = 1)

    # * * init phi ----
    phi <- list()

    # might want to estimate the random effects?
    for (i in seq(1, length(y))) {
        phi[[i]] <- list(
            mu = c(beta[[i]], rep(0, U)),
            Sigma = theta$Tau
        )
    }

    # * * init pi ----
    pi <- rep(0.8, length(y))

    # * * collect inits ----
    init <- list(
        theta = theta,
        phi = phi,
        pi = pi
    )
    # * priors ----
    priors <- list(a_beta = a.beta, b_beta = b.beta,
                        a_u = 2, b_u = 2)


    list(
        data = data,
        init = init,
        priors = priors
    )
}