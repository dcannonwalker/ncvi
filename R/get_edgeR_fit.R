get_edgeR_p <- function(data) {
  y <- data$y
  counts <- matrix(unlist(y), nrow = length(y), byrow = T)
  design <- data$C1
  group <- factor(design[, 2])
  data_edgeR <- edgeR::DGEList(counts = counts, group = group)
  data_edgeR <- edgeR::estimateDisp(data_edgeR, design)
  edgeR_fit <- edgeR::glmFit(data_edgeR, design = design)
  lrt <- edgeR::glmLRT(edgeR_fit)
  lrt$table$PValue
}
