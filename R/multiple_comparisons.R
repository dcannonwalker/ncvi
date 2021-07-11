bfdr <- function(post_prob) {
  p <- sort(post_prob, decreasing = F)
  n <- 1:length(p)
  sapply(n, function(x) {
    sum(p[1:x]) / x
  })
}

tfdr <- function(post_prob, true_null) {
  tn <- true_null[order(post_prob)]
  n <- 1:length(tn)
  sapply(n, function(x) {
    sum(tn[1:x]) / x
  })
}

tpr <- function(post_prob, true_null) {
  tp <- 1 - true_null[order(post_prob)]
  p <- post_prob[order(post_prob)]
  n <- 1:length(p)
  sapply(n, function(x) {
    sum(tp[1:x]) / sum(tp)
  })
}

fpr <- function(post_prob, true_null) {
  tn <- true_null[order(post_prob)]
  p <- post_prob[order(post_prob)]
  n <- 1:length(p)
  sapply(n, function(x) {
    sum(tn[1:x]) / sum(tn)
  })
}
