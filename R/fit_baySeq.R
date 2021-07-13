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

fit_baySeqRibo <- function(data_baySeq, options = NULL) {
  CD <- new("countData", data = data_baySeq$data,
            replicates = data_baySeq$replicates,
            groups = data_baySeq$groups)
  baySeq::libsizes(CD) <- baySeq::getLibsizes(CD)
  CD@annotation <- data.frame(name = paste("count",
                                           1:data_baySeq$G, sep = "_"))
  CD <- baySeq::getPriors.NB(CD, samplesize = options$samplesize,
                             estimation = "QL", cl = cl)
  CD <- baySeq::getLikelihoods(CD, cl = cl, bootStraps = 3)
  named_p <- baySeq::topCounts(CD, group = "trsl", number = data_baySeq$G) %>%
    tidyr::separate(name, into = c(NA, "name"), convert = T) %>%
    dplyr::arrange(name) %>%
    dplyr::select(name, likes)
  list(fit = CD, named_p = named_p, p = named_p$likes, type_str = "baySeq",
       true_null = data_baySeq$true_null)
}

# replicates <- as.factor(rep(c("A1", "A2", "B1", "B2"), each = 2))
# NDE <- factor(rep(1, 8))
# trt <- c(rep("A", 4), rep("B", 4))
# prot <- c(rep(c("ribo", "rna"), 4))
# trsl <- rep(c("A", "B", "C", "D"), each = 2)
# groups <- list(NDE = NDE, trt = trt, prot = prot, trsl = trsl)
