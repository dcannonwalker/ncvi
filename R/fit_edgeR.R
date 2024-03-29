# fit_edgeR <- function(data_edgeR) {
#   y <- edgeR::DGEList(counts = data_edgeR$counts,
#                       group = data_edgeR$group)
#   y <- edgeR::calcNormFactors(y)
#
#   y <- edgeR::estimateDisp(y, design = data_edgeR$design)
#
#   fit <- edgeR::glmFit(y, design = data_edgeR$design)
#   lrt <- edgeR::glmLRT(fit)
#   qobj <- qvalue::qvalue(lrt$table$PValue)
#   q <- qobj$qvalues
#   list(fit = fit, lrt = lrt, type_str = "edgeR",
#        true_null = data_edgeR$true_null,
#        p = p.adjust(lrt$table$PValue), q = q)
# }
#
# fit_edgeRV <- function(data_edgeR) {
#   y <- edgeR::DGEList(counts = data_edgeR$counts,
#                       group = data_edgeR$group)
#   y <- edgeR::calcNormFactors(y)
#
#   y <- edgeR::estimateCommonDisp(y)
#   tests <- edgeR::exactTest(y)
#   p <- tests$table$PValue
#   padj <- p.adjust(p, method="BH")
#
#   list(partFit = y, tests = tests, p = p, padj = padj)
# }

fit_edgeRribo <- function(data_edgeR) {
  design = data_edgeR$design
  y <- edgeR::DGEList(counts = data_edgeR$counts,
                      group = data_edgeR$group)
  y <- edgeR::calcNormFactors(y)
  y <- edgeR::estimateDisp(y, design = design)
  fit <- edgeR::glmQLFit(y, design = design)
  test <- edgeR::glmQLFTest(fit, coef = ncol(design))

  list(table = test$table, p = test$table$PValue,
       padj = p.adjust(test$table$PValue, method = "BH"))

}

