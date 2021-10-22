bfdr <- function(post_prob, return.names = F) {
  p <- sort(post_prob, decreasing = F)
  n <- 1:length(p)
  bfdr <- sapply(n, function(x) {
    sum(p[1:x]) / x
  })
  if(return.names) {
    names <- 1:length(post_prob)
    names <- names[order(post_prob)]
    cbind(bfdr, names)
  }
  else bfdr
}

tfdr <- function(p, true_null) {
  tn <- true_null[order(p)]
  n <- 1:length(tn)
  sapply(n, function(x) {
    sum(tn[1:x]) / x
  })
}

tpr <- function(p, true_null) {
  tp <- 1 - true_null[order(p)]
  p <- p[order(p)]
  n <- 1:length(p)
  sapply(n, function(x) {
    sum(tp[1:x]) / sum(tp)
  })
}

fpr <- function(p, true_null) {
  tn <- true_null[order(p)]
  p <- p[order(p)]
  n <- 1:length(p)
  sapply(n, function(x) {
    sum(tn[1:x]) / sum(tn)
  })
}

get_all_multiple_comparisons <- function(p, true_null, padj = NULL,
                                         method = NULL, gene = NULL) {
  if (is.null(method)) {
    list(
      p = p[order(p)],
      fdr = bfdr(p),
      tfdr = tfdr(p, true_null),
      tpr = tpr(p, true_null),
      fpr = fpr(p, true_null)
  )
  }
  else if (method == "edgeR") {
    list(
      p = p[order(p)],
      fdr = padj[order(p)],
      tfdr = tfdr(p, true_null),
      tpr = tpr(p, true_null),
      fpr = fpr(p, true_null)
    )
  }

  else if (method == "xtail") {
    true_null <- true_null[gene]
    list(
      p = p[order(p)],
      fdr = padj[order(p)],
      tfdr = tfdr(p, true_null),
      tpr = tpr(p, true_null),
      fpr = fpr(p, true_null)
    )
  }
}
