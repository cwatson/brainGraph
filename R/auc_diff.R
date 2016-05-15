#' Difference in the area-under-the-curve of two vectors
#'
#' This function takes two vectors, calculates the area-under-the-curve (AUC),
#' and calculates the difference between the two.
#'
#' @param x Numeric vector of the x-values
#' @param y A 2-column numeric matrix; each column contains the values for one
#'   group
#'
#' @return A numeric value of the difference between two groups

auc_diff <- function(x, y) {
  return(-diff(apply(y, 2, function(z)
                     sum(diff(x) * (head(z, -1) + tail(z, -1))) / 2)))
}
