#' Transform a vector to have a different range.
#'
#' This function takes a vector and transforms it to have a new range, given
#' the input, or the default values of [0, 1].
#'
#' @param x the vector to transform
#' @param min.val the minimum value of the new range
#' @param max.val the maximum value of the new range
#' @export
#'
#' @return A vector of the transformed input.

range.transform <- function(x, min.val=0, max.val=1) {
  if (max.val==1) {
    (x - min(x)) / diff(range(x))
  } else {
    ((x - min(x)) * (max.val-1) / diff(range(x))) + min.val
  }
}
