#' Check model residuals for each brain region.
#'
#' This function checks the model residuals for each brain region in the
#' analysis. It simply does a qqplot of the studentized residuals.
#'
#' @param resids Data table of model residuals for all brain regions
#' @export
#'
#' @seealso \code{\link[stats]{qqnorm}, \link[stats]{qqline}}

check.resid <- function(resids) {
  kNumRegions <- ncol(resids) - 1

  kNumPlots <- 9
  a <- kNumRegions %/% kNumPlots
  b <- kNumRegions %% kNumPlots

  for (j in seq_len(a)) {
    dev.new()
    par(mfrow=c(3, 3))
    for (i in seq_len(kNumPlots)) {
      ind <- kNumPlots * (j - 1) + i
      qqnorm(resids[[ind]], main=colnames(resids)[ind])
      qqline(resids[[ind]])
      title(sub=ind)
    }
  }
  dev.new()
  par(mfrow=c(3, 3))
  for (k in seq_len(b)) {
    ind <- kNumPlots * j + k
    qqnorm(resids[[ind]], main=colnames(resids)[ind])
    qqline(resids[[ind]])
    title(sub=ind)
  }
}
