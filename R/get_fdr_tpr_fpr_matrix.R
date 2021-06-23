
get_fdr_tpr_fpr_df <- function(fit) {
  p <- unlist(fit$pars$pi)
  beta1 <- fit$data$beta[, 2]
  true_positive <- beta1 != 0
  true_negative <- beta1 == 0
  df <- data.frame(p = p, true_positive = true_positive,
                   true_negative = true_negative)
  df <- dplyr::arrange(df, p)
  data.frame(p = df$p,
             tpr = get_tpr_vector(df),
             fpr = get_fpr_vector(df),
             fdr = get_fdr(df),
             true_fdr = get_true_fdr_vector(df),
             true_positive = df$true_positive,
             true_negative = df$true_negative)
}

get_fdr_tpr_fpr_matrix <- function(fit) {
  p <- unlist(fit$pars$pi)
  beta1 <- fit$data$beta[, 2]
  true_positive <- as.numeric(beta1 != 0)
  true_negative <- as.numeric(beta1 == 0)
  tag <- seq(1, fit$data$G)
  mat <- cbind(p = p, true_positive = true_positive,
                   true_negative = true_negative, tag = tag)
  mat <- mat[order(mat[, "p"]), ]

  cbind(mat,
        tpr = get_tpr_vector(mat),
        fpr = get_fpr_vector(mat),
        fdr = get_fdr_vector(mat),
        true_fdr = get_true_fdr_vector(mat))

}

