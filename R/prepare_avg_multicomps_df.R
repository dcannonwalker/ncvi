prepare_avg_multicomps_df <- function(fits, methods, true_nulls,
                                      alpha, tfdr = T) {
  mat <- matrix(ncol = 2, nrow = length(alpha) * length(fits))
  for (i in seq(1, length(fits))) {
    std_multicomps <- mapply(function(x, y) {
      extract_std_multicomps(x, y, method = methods[i],
                             alpha = alpha, tfdr)
    }, fits[[i]], true_nulls, SIMPLIFY = T)
    avg_multicomps <- apply(std_multicomps, 1, mean, na.rm = T)
    mat[((i - 1) * length(alpha) + 1):(i * length(alpha)), ] <-
      cbind(avg_multicomps, alpha)
  }

  df <- data.frame(mat, rep(methods, each = length(alpha)))
  if (tfdr == T) {
    colnames(df) <- c("tfdr", "fdr", "method")
  }
  else colnames(df) <- c("tpr", "fpr", "method")

  df
}
