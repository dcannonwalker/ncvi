get_fdr_vector <- function(mat) {
  p_vector <- mat[, "p"]
  fdr <- sapply(seq(1, length(p_vector)), function(n) {
    sum(p_vector[1:n]) / n
  })
  fdr
}

get_true_fdr <- function(mat, n) {
  sum(mat[1:n, "true_negative"]) / n
}

get_true_fdr_vector <- function(mat) {
  sapply(seq(1, nrow(mat)), function(n) get_true_fdr(mat, n))
}


