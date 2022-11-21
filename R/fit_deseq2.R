#' Function to prepare the counts table and annotations
#' required by `DESeq2`
#' @param other_data List, as produced by `sim_data_mixture()`
#' @return A list with `counts` and `experiment`
prep_deseq <- function(other_data) {
  counts <- other_data$counts
  design <- other_data$design
  N <- nrow(design)
  Data_Type <- factor(design[, 'preparation'],
                      labels = c("Ribo", "RNA"))
  Conditions <- factor(design[, 'treatment'],
                       labels = c("Control", "Treatment"))
  Replicates <- rep(1:(N / 4), 4)
  Samples  <- paste(Data_Type, Conditions, Replicates, sep = "")
  Bio_Replicate2 <- factor(c(0, 1, rep(0, N / 2 - 2), 0, 1, rep(0, N / 2 - 2)))
  Bio_Replicates <- data.frame(Bio_Replicate2)
  for (i in seq(3, (N / 2) - 1)) {
    Bio_Replicates[, paste0("Bio_Replicate", i)] <-
      factor(rep(c(rep(0, i - 1), 1, rep(0, N / 2 - i)), 2))
  }
  experiment <- data.frame(Data_Type,
                           Bio_Replicates,
                           Conditions,
                           row.names = Samples)
  colnames(counts) <- Samples
  design_formula <- as.formula(paste("~ ", paste(
    "Data_Type", "Conditions",
    paste(paste0("Bio_Replicate", 2:(N / 2 - 1)), collapse = " + "),
    "Data_Type:Conditions", sep = " + "
  )))
  list(experiment = experiment, counts = counts,
       design_formula = design_formula)
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
                                        design = data_deseq$design_formula
                                          )
  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds)
  list(dds = dds, res = res)
}


