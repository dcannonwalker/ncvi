#' Function to prepare the counts table and annotations
#' required by `DESeq2`
#' @param data_edgeR List, as produced by `sim_data_mixture()`
#' @return A list with `counts` and `experiment`
prep_deseq <- function(data_edgeR) {
  counts <- data_edgeR$counts
  N <- nrow(data_edgeR$design)
  Data_Type <- rep(c("Ribo", "RNA"), N / 2)
  Conditions <- rep(c("Control", "Treatment"), each = N / 2)
  Replicates <- rep(c(1, 1, 2, 2), N / 4)
  Samples  <- paste(Data_Type, Conditions, Replicates, sep = "")
  experiment <- data.frame(Data_Type, Conditions, row.names = Samples)
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
                                          Data_Type + Conditions:Data_Type)
  dds <- DESeq2::DESeq(dds)
  res <- results(dds, name = "ConditionsTreatment.Data_TypeRNA")
  list(dds = dds, res = res)
}
