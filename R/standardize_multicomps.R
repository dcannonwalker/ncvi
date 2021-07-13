standardize_multicomps <- function(v1, v2, alpha) {

  intervals <- findInterval(v1, alpha, rightmost.closed = T)

  populated_interval_means <- tapply(v2, as.factor(intervals), mean)
  interval_means <- rep(NA, length(alpha))
  interval_means[unique(intervals)] <- populated_interval_means

  interval_means

}
