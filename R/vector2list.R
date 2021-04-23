#' Transform a vector of length `N * G`
#' to a list of vectors each of length `N`
#'
#' @param y The original vector
#' @param G Integer for the number of groups
#' @param N Integer for the length of each group
#' @export
vector2list <- function(y, G, N) {
  y_list <- list()

  for (i in seq(1, G)) {

    y_list[[i]] <- y[((i - 1) * N + 1):(i * N)]

  }

  y_list
}
