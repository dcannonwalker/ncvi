#' Function to take the simulated data set and format it
#' as required for `xtail`
#'
#' @param sim output of `sim_data_mixture()`
#' @return list with RPF and mRNA data frames and
#' vector of treatment conditions
prep_xtail <- function(sim) {

  y_df <- sim$data$y |> list2DF() |> data.table::transpose()
  rpf <- y_df[, 2 * (0:(nrow(sim$data$C) / 2 - 1)) + 1]
  mrna <- y_df[, 2 * (1:(nrow(sim$data$C) / 2))]

  condition <- c(rep("control", nrow(sim$data$C) / 4),
                 rep("treatment", nrow(sim$data$C) / 4))
  colnames(rpf) <- c(paste0("control", 1:(nrow(sim$data$C) / 4)),
                      paste0("trt", 1:(nrow(sim$data$C) / 4)))
  colnames(mrna) <- c(paste0("control", 1:(nrow(sim$data$C) / 4)),
                     paste0("trt", 1:(nrow(sim$data$C) / 4)))

  list(rpf = rpf, mrna = mrna, condition = condition)

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



