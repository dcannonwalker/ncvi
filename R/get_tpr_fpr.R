get_tpr <- function(mat, n) {
    sum(mat[1:n, "true_positive"]) / sum(mat[, "true_positive"])
}

get_fpr <- function(mat, n) {
  sum(mat[1:n, "true_negative"]) / sum(mat[, "true_negative"])
}

get_tpr_vector <- function(mat) {
  sapply(seq(1, nrow(mat)), function(n) get_tpr(mat, n))
}

get_fpr_vector <- function(mat) {
  sapply(seq(1, nrow(mat)), function(n) get_fpr(mat, n))
}



