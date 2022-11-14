#' A function to calculate tagwise dispersion based on empirical
#' function from real data
#' @param beta0
dispersion_est <- function(beta0) {
  logCPM <- beta0 - 1.36
  which_disp <- which.min(abs(logCPM - dispersion_function$AveLogCPM))
  disp <- dispersion_function$dispersion[which_disp]
  disp
}

