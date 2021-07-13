extract_std_multicomps <- function(fit, true_null, method, alpha, tfdr = T) {
  if (is.null(method)) message("method required")

  else if (method == "ncvi") {
    p <- fit$pars$pi
  }

  else if (method == "mcmc") {
    post <- posterior(fit)
    p  <- post %>% select(starts_with("D", ignore.case = F)) %>%
      colMeans()
  }

  else if (method == "baySeq") {
    p <- fit$p
  }

  else if (method == "edgeR") {
    p <- fit$padj
  }
  if (method != "edgeR") method <- NULL
  mult <- get_all_multiple_comparisons(p, true_null, method = method)
  if (tfdr == T) {
    standardize_multicomps(mult$fdr, mult$tfdr, alpha)
  }
  else if (tfdr == F) {
    standardize_multicomps(mult$fpr, mult$tpr, alpha)
  }

}
