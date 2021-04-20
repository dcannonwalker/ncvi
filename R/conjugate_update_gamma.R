## updates for Gamma-distributed precision parameters

conjugate_update_gamma <- function() {
  a <- a + sum(S)/2
  b <- b + 0.5 * (
    sum(mu^2 + diag(Sigma) + M^2 + (mu * M) + diag(R))
    )
}
