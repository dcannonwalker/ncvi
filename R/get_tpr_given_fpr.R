get_tpr_given_fpr <- function(mat, alpha) {
  mat_filter <- mat[mat[, "fpr"] <= alpha, ]
  if (is.null(nrow(mat_filter))) matrix(mat_filter[c("fpr", "tpr")],
                                        nrow = 1)
  else if(nrow(mat_filter) == 0) c(fpr = 0, tpr = 0)
  else {mat_filter[nrow(mat_filter), c("fpr", "tpr")]}
}

get_tpr_given_fpr_vector <- function(mat, alpha_vector) {
  mat <- sapply(alpha_vector, function(a) {
    get_tpr_given_fpr(mat, a)
  })
  cbind(t(mat), fpr = alpha_vector)
}

get_tpr_list_given_fpr_vector <- function(fit_list, alpha_vector) {
  lapply(fit_list, function(f) {
    mat <- get_fdr_tpr_fpr_matrix(f)
    get_tpr_given_fpr_vector(mat, alpha = alpha_vector)
  })
}

get_tpr_vector_edgeR <- function(fit, alpha_vector) {
  p <- get_edgeR_p(fit$data)
  beta1 <- fit$data$beta[, 2]
  true_positive <- as.numeric(beta1 != 0)
  true_negative <- as.numeric(beta1 == 0)
  tag <- seq(1, fit$data$G)
  mat <- cbind(p = p, true_positive = true_positive,
               true_negative = true_negative, tag = tag)
  mat <- mat[order(mat[, "p"]), ]

  mat <- cbind(mat,
               tpr = get_tpr_vector(mat),
               fpr = get_fpr_vector(mat),
               fdr = get_fdr_vector(mat),
               true_fdr = get_true_fdr_vector(mat))

  get_tpr_given_fpr_vector(mat, alpha = alpha_vector)
}
get_tpr_list_edgeR <- function(fit_list, alpha_vector) {
  lapply(fit_list, function(f) {
    get_tpr_vector_edgeR(f, alpha_vector)
  })
}
