#' Function to take the simulated data set and format it
#' as required for `xtail`
#'
#' @param data_edgeR output element of `sim_data_mixture()`
#' @return list with RPF and mRNA data frames and
#' vector of treatment conditions
prep_xtail <- function(data_edgeR) {


  counts <- data_edgeR$counts
  design <- data_edgeR$design
  N <- nrow(design)
  Data_Type <- factor(design[, 'preparation'],
                      labels = c("Ribo", "RNA"))
  Conditions <- factor(design[, 'treatment'],
                       labels = c("Control", "Treatment"))[1:(N / 2)]
  Replicates <- rep(c(1:(N / 4)), 2)
  Samples  <- paste(Conditions, Replicates, sep = "")
  mrna <- counts[, Data_Type == "RNA"]
  rpf <- counts[, Data_Type == "Ribo"]
  colnames(mrna) <- colnames(rpf) <- Samples


  list(rpf = rpf, mrna = mrna, condition = as.character(Conditions),
       Conditions = Conditions)

}

#' Function to fit `xtail` model to paired RiboSeq and RNA count data
#' @param data_xtail list with `rpf` and `mrna` data frames and
#' character vector `condition` of treatment conditions
#' @param bins `bins` argument for `xtail()`
#' @param threads `threads` argument for `xtail()`
#' @param normalize `normalize` argument for `xtail()`
#' @return `xtail` fit object
fit_xtail <- function(data_xtail, bins, threads = 2,
                      normalize = T,
                      minMeanCount = 1) {

  results <- xtail::xtail(data_xtail$mrna,
                          data_xtail$rpf,
                          data_xtail$condition,
                          bins = bins,
                          threads = threads,
                          normalize = normalize,
                          minMeanCount = minMeanCount)

 list(results = results, p = results$resultsTable$pvalue_final,
      padj = results$resultsTable$pvalue.adjust,
      gene = as.numeric(rownames(results$resultsTable)))

}



