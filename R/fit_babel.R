#' Function to prepare data as required by `babel`
#' @param other_data
#' @return A list with `ribo`, `rna`, and `group`
prep_babel <- function(other_data, spec_Replicates) {
  counts <- other_data$counts
  design <- other_data$design
  N <- nrow(design)
  Data_Type <- factor(design[, 'preparation'],
                      labels = c("Ribo", "RNA"))
  Conditions <- factor(design[, 'treatment'],
                       labels = c("Control", "Treatment"))[1:(N / 2)]
  if(!missing(spec_Replicates)) Replicates <- spec_Replicates
  else Replicates <- rep(c(1:(N / 4)), 2)
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
