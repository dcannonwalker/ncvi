#' make list of lists with means, lower bound,
#'  and upper bound for 1-alpha confidence intervals
#'  for params in `pars$phi`
#'
#' @param fit A fit produced by `fit_ncvi()`
#' @param alpha Confidence level
#' @export
make_ci_list <- function(fit, alpha = 0.05) {
  lapply(fit$pars$phi, function(x){
    mean <- x$mu
    se <- sqrt(diag(x$Sigma))*qnorm(1-alpha/2)
    lower <- mean - se
    upper <- mean + se
    list(mean = mean, lower = lower, upper = upper)
  })
}

#' Construct data frame of confidence intervals
#'
#' @param fit
#' @export
make_ci <- function(fit, alpha = 0.05, par_names = NULL) {
  ci_list <- make_ci_list(fit, alpha)
  data <- fit$data
  G <- data$G
  P <- data$P
  U <- data$U
  intervals_df <- data.frame(param = c(paste0("beta",
                                              c(0,1:(P-1))),
                                       paste0("u", 1:U)),
                             ci_list)

  colnames(intervals_df) <-
    c("param",
      t(sapply(c("mean.",
                 "lower.",
                 "upper."),
               function(x) paste0(x, 1:G))))

  intervals_df <-
    tidyr::pivot_longer(intervals_df,
                        cols = -param,
                        names_to =
                          c("which_estimate",
                            "symbol"),
                        names_sep = "\\.",
                        values_to = "estimate",
                        names_transform =
                          list(symbol = as.integer)) %>%
    tidyr::pivot_wider(id_cols = c(symbol, param),
                       names_from = which_estimate,
                       values_from = estimate) %>%
    dplyr::arrange(symbol)

  if (!is.null(par_names)) {
    intervals_df %>% dplyr::filter(param %in% par_names)
  }

  else intervals_df
}

#' Perform a hypothesis test based on posterior
#'
#' @param fit Fit from `fit_ncvi()`
#' @param cutoff Threshold value
make_test_list <- function(fit, cutoff) {
  lapply(fit$pars$phi, function(x) {
    mean <- x$mu
    sd <- sqrt(diag(x$Sigma))
    list(central = 1 - (pnorm(-cutoff, mean = mean, sd = sd, lower.tail = T) +
                    pnorm(cutoff, mean = mean, sd = sd, lower.tail = F)),
         lowertail = pnorm(cutoff, mean = mean, sd = sd, lower.tail = T),
         uppertail = pnorm(-cutoff, mean = mean, sd = sd, lower.tail = F))
  })
}

make_test <- function(fit, cutoff, par_names = NULL) {
  test_list <- make_test_list(fit, cutoff)
  data <- fit$data
  G <- data$G
  P <- data$P
  U <- data$U
  test_df <- data.frame(param = c(paste0("beta",
                                              c(0,1:(P-1))),
                                       paste0("u", 1:U)),
                             test_list)

  colnames(test_df) <-
    c("param",
      t(sapply(c("central.",
                 "lowertail.",
                 "uppertail."),
               function(x) paste0(x, 1:G))))

  test_df <-
    tidyr::pivot_longer(test_df,
                        cols = -param,
                        names_to =
                          c("which_test",
                            "symbol"),
                        names_sep = "\\.",
                        values_to = "p",
                        names_transform =
                          list(symbol = as.integer)) %>%
    tidyr::pivot_wider(id_cols = c(symbol, param),
                       names_from = which_test,
                       values_from = p) %>%
    dplyr::arrange(symbol)

  if (!is.null(par_names)) {
    test_df %>% dplyr::filter(param %in% par_names)
  }

  else test_df
}
