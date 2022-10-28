#' Expects samples to be organized as Ribo (ctrl, trt), RNA (ctrl, trt)
prep_bayseq <- function(data_edgeR) {
  N <- nrow(data_edgeR$design)
  counts <- as.matrix(data_edgeR$counts)
  dim = c(nrow(counts), N / 2, 2)
  replicates = rep(c(1, 2), N / 4)
  groups = list(NDE = rep(1, N / 2),
                DE = rep(c(1, 2), each = N / 4))
  options <- list(samplesize = 1000)
  bayseq_array <- array(c(counts[, 1:(N / 2)], counts[, (N / 2 + 1):N]), dim = dim)

  data_bayseq <- new("countData", data = bayseq_array,
                     replicates = replicates,
                     groups = groups,
                     densityFunction = bbDensity)

  libsizes(data_bayseq) <- getLibsizes(data_bayseq)
  data_bayseq
}

fit_bayseq <- function(data_bayseq, ncl) {
  G <- nrow(data_bayseq@data)
  cl = parallel::makeCluster(ncl)
  data_bayseq <- baySeq::getPriors(data_bayseq, samplesize = 1000, cl = cl)
  data_bayseq <- baySeq::getLikelihoods(data_bayseq, pET = "BIC", cl = cl)
  parallel::stopCluster(cl)
  # look at results --------------------------------------------------------
  data_bayseq@annotation <- data.frame(name = paste("count",
                                                    1:G,
                                                    sep = "_"))
  list(tC = baySeq::topCounts(data_bayseq, number = G, group = "DE"),
       data_bayseq = data_bayseq)
}
