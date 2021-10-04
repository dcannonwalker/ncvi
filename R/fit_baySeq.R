#' This function may be unnecessary
make_baySeq_data <- function()
fit_baySeq <- function(data_baySeq, options = NULL) {
  # this requires baySeq to be loaded
  CD <- new("countData", data = data_baySeq$data,
            replicates = data_baySeq$replicates,
            groups = data_baySeq$groups)
  baySeq::libsizes(CD) <- baySeq::getLibsizes(CD)
  CD@annotation <- data.frame(name = paste("count",
                                           1:data_baySeq$G, sep = "_"))
  CD <- baySeq::getPriors.NB(CD, samplesize = options$samplesize,
                             estimation = "QL", cl = cl)
  CD <- baySeq::getLikelihoods(CD, cl = cl, bootStraps = 3)
  named_p <- baySeq::topCounts(CD, group = "NDE", number = data_baySeq$G) %>%
    tidyr::separate(name, into = c(NA, "name"), convert = T) %>%
    dplyr::arrange(name) %>%
    dplyr::select(name, likes)
  list(fit = CD, named_p = named_p, p = named_p$likes, type_str = "baySeq",
       true_null = data_baySeq$true_null)
}

fit_baySeqRibo <- function(data_baySeq, options = NULL, ncl = 2) {
  cl = parallel::makeCluster(ncl)

  on.exit(parallel::stopCluster(cl))

  CD <- new("countData", data = array(c(data_baySeq$counts[, 1:8],
                                        data_baySeq$counts[, 9:16]),
                                      dim = c(nrow(data_baySeq$counts),
                                                   8, 2)),
                                      replicates = c(1, 1, 1, 1,
                                                     2, 2, 2, 2),
                                      groups = list(NDE = rep(1, 8),
                                                    DE = rep(c(1, 2),
                                                             each = 4)),
            densityFunction = bbDensity)

  baySeq::libsizes(CD) <- baySeq::getLibsizes(CD)

  CD@annotation <- data.frame(name = paste("count",
                                           1:nrow(data_baySeq$counts),
                                           sep = "_"))

  CD <- baySeq::getPriors(CD, samplesize = options$samplesize,
                          cl = cl)

  CD <- baySeq::getLikelihoods(CD, pET = 'BIC',
                               nullData = TRUE,
                               cl = cl)

  named_p <- baySeq::topCounts(CD, group = "NDE",
                               number = nrow(data_baySeq$counts)) %>%
    tidyr::separate(name, into = c(NA, "name"), convert = T) %>%
    dplyr::arrange(name) %>%
    dplyr::select(name, likes)

  list(fit = CD, named_p = named_p, p = named_p$likes, type_str = "baySeq")

}

# replicates <- as.factor(rep(c("A1", "A2", "B1", "B2"), each = 2))
# NDE <- factor(rep(1, 8))
# trt <- c(rep("A", 4), rep("B", 4))
# prot <- c(rep(c("ribo", "rna"), 4))
# trsl <- rep(c("A", "B", "C", "D"), each = 2)
# groups <- list(NDE = NDE, trt = trt, prot = prot, trsl = trsl)