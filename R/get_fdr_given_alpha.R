get_fdr_given_alpha <- function(mat, alpha) {
  mat_filter <- mat[mat[, "p"] <= alpha, ]
  if (is.null(nrow(mat_filter))) matrix(mat_filter[c("fdr", "true_fdr")], nrow = 1)
  else if(nrow(mat_filter) == 0) c(fdr = 0, true_fdr = 0)
  else {mat_filter[nrow(mat_filter), c("fdr", "true_fdr")]}
}

get_fdr_given_alpha_vector <- function(mat, alpha_vector) {
  mat <- sapply(alpha_vector, function(a) {
    get_fdr_given_alpha(mat, a)
  })
  cbind(t(mat), alpha_vector)
}

# ggplot2::ggplot(df, ggplot2::aes(fpr, tpr)) + ggplot2::geom_line()
