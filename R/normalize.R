#' Normalize count table using DESeq2
#'
#' @param counts Data frame of counts with column `Geneid`
#' @param group Vector of group labels for each sample
#'
#' @export
normalize <- function(counts, group) {
    Geneid <- counts$Geneid
    counts <- as.matrix(counts[, 2:ncol(counts)])
    samples <- data.frame(group = group,
                          row.names = colnames(counts))

    ds <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
                                         colData = samples,
                                         design = ~ group)
    ds <- DESeq2::estimateSizeFactors(ds)
    normed <- DESeq2::counts(ds, normalized = T)

    data.frame(Geneid = Geneid, normed)
}
