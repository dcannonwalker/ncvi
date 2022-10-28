#' Function to prepare the counts table and annotations
#' required by `DESeq2`
#' @param data_edgeR List, as produced by `sim_data_mixture()`
#' @return A list with `counts` and `experiment`
prep_deseq <- function(data_edgeR) {
  counts <- data_edgeR$counts
  design <- data_edgeR$design
  N <- nrow(design)
  Data_Type <- factor(design[, 'preparation'],
                      labels = c("Ribo", "RNA"))
  Conditions <- factor(design[, 'treatment'],
                       labels = c("Control", "Treatment"))
  Replicates <- rep(c(1, 2), N / 2)
  Samples  <- paste(Data_Type, Conditions, Replicates, sep = "")
  Bio_Replicate2 <- c(0, 1, rep(0, N / 2 - 2), 0, 1, rep(0, N / 2 - 2))
  Bio_Replicate3 <- c(0, 0, 1, rep(0, N / 2 - 3), 0, 0, 1, rep(0, N / 2 - 3))
  experiment <- data.frame(Data_Type, Conditions,
                           Bio_Replicate2,
                           Bio_Replicate3, row.names = Samples)
  colnames(counts) <- Samples
  list(experiment = experiment, counts = counts)
}

#' Function to fit the `DESeq2` model to a RiboSeq experiment
#' and return results for the interaction effect
#' @param `data_deseq` A list as produced by `prep_deseq()`
#' @return A list with `dds` (a `DESeq` data object with
#' attached results) and `res` from `results()` with the interaction effect
#' parameter
fit_deseq <- function(data_deseq) {
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = data_deseq$counts,
                                        colData = data_deseq$experiment,
                                        design = ~ Conditions +
                                          Data_Type + Conditions:Data_Type +
                                          Bio_Replicate2 +
                                          Bio_Replicate3,)
  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds, name = "ConditionsTreatment.Data_TypeRNA")
  list(dds = dds, res = res)
}
