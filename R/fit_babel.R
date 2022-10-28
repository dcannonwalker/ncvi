#' Function to prepare data as required by `babel`
#' @param data_edgeR
#' @return A list with `ribo`, `rna`, and `group`
prep_babel <- function(data_edgeR) {
  counts <- data_edgeR$counts
  design <- data_edgeR$design
  N <- nrow(design)
  Data_Type <- factor(design[, 'preparation'],
                      labels = c("Ribo", "RNA"))
  Conditions <- factor(design[, 'treatment'],
                       labels = c("Control", "Treatment"))[1:(N / 2)]
  Replicates <- rep(c(1:(N / 4)), 2)
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
  res <- babel::babel(rna = data_babel$rna, rp = data_babel$ribo,
               group = data_babel$group,
               nreps = nreps)
  res
}
