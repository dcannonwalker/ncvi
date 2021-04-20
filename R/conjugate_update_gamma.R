## updates for Gamma-distributed precision parameters

conjugate_update_gamma <- function() {
  b <- b + 0.5 * (
    sum(mu^2 + sig2 + M^2 + (mu * M) + r2)
    )
}
