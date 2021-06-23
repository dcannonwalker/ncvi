fdr <- function(p_vector, alpha) {

  # subset of p values declared significant at alpha level
  significant <- p_vector[p_vector < alpha]

  # false discovery rate
  sum(significant) / length(significant)

}

fdr_vector <- function(p_vector, alpha_vector) {

  sapply(alpha_vector, function(s) fdr(p_vector, s))

}

