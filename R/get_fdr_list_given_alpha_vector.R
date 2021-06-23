get_fdr_list_given_alpha_vector <- function(fit_list, alpha_vector) {
  lapply(fit_list, function(f) {
    mat <- get_fdr_tpr_fpr_matrix(f)
    get_fdr_given_alpha_vector(mat, alpha = alpha_vector)
    })
}

get_fdr_vector_given_alpha_vector <- function(fit, alpha_vector) {
  mat <- get_fdr_tpr_fpr_matrix(fit)
  get_fdr_given_alpha_vector(mat, alpha = alpha_vector)
}

