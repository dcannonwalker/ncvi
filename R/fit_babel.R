#' Function to prepare data as required by `babel`
#' @param data_edgeR
#' @return A list with `ribo`, `rna`, and `group`
prep_babel <- function(data_edgeR) {
  counts <- data_edgeR$counts
  N <- nrow(data_edgeR$design)
  Data_Type <- rep(c("Ribo", "RNA"), N / 2)
  Conditions <- rep(c("Control", "Treatment"), each = N / 4)
  Replicates <- rep(1:(N / 4), 2)
  Samples  <- paste(Conditions, Replicates, sep = "")
  rna <- counts[, Data_Type == "RNA"]
  ribo <- counts[, Data_Type == "Ribo"]
  colnames(rna) <- colnames(ribo) <- Samples
  list(ribo = ribo, rna = rna, group = Conditions)
}

#' Function to fit the `babel` model to a RiboSeq experiment
#' @param data_babel Output of `prep_babel()`
#' @param nreps Integer, number of mc reps
#' @param mc.cores Integer, to set `options(mc.cores = ...)`
#' @return A list, the output of `babel::babel()`
fit_babel <- function(data_babel, nreps, mc.cores = 1) {
  options(mc.cores = mc.cores)
  res <- babel(rna = data_babel$rna, rp = data_babel$ribo,
               group = data_babel$group,
               nreps = nreps)
  res
}
