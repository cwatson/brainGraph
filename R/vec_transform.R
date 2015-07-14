#' Transform a vector to have a different range
#'
#' This function takes a vector and transforms it to have a new range, given
#' the input, or the default values of [0, 1].
#'
#' @param x the vector to transform
#' @param min.val the minimum value of the new range
#' @param max.val the maximum value of the new range
#'
#' @return A vector of the transformed input.

vec.transform <- function(x, min.val=0, max.val=1) {
  ((x - min(x, na.rm=T)) * (max.val - min.val) / diff(range(x, na.rm=T))) + min.val
}
