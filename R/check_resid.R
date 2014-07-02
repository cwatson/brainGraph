#' Check model residuals for each brain region.
#'
#' This function checks the model residuals for each brain region in the
#' analysis. It simply does a qqplot of the studentized residuals.
#'
#' @param resids data frame of model residuals for all brain regions
#' @export
#'

check.resid <- function(resids) {
  kNumRegions <- dim(resids)[2] - 1

  a <- kNumRegions %/% 9
  b <- kNumRegions %% 9

  for (j in 1:a) {
    dev.new()
    par(mfrow=c(3, 3))
    for (i in 1:9) {
      qqnorm(resids[, 9*(j-1) + i], main=colnames(resids)[9*(j-1) + i])
      qqline(resids[, 9*(j-1) + i])
      title(sub=9*(j-1)+i)
    }
  }
  dev.new()
  par(mfrow=c(3, 3))
  for (k in 1:b) {
    qqnorm(resids[, 9*j + k], main=colnames(resids)[9*j + k])
    qqline(resids[, 9*j + k])
    title(sub=9*j+k)
  }
}
