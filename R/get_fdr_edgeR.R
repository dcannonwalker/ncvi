# get fdr averages for edgeR fit
get_fdr_vector_edgeR <- function(fit, alpha_vector) {
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

  get_fdr_given_alpha_vector(mat, alpha = alpha_vector)
}

get_fdr_list_edgeR <- function(fit_list, alpha_vector) {
  lapply(fit_list, function(f) {
    get_fdr_vector_edgeR(f, alpha_vector)
  })
}
