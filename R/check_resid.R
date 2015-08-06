#' Check model residuals for each brain region
#'
#' This function checks the model residuals for each brain region in the
#' analysis. It simply does a qqplot of the studentized residuals (but uses
#' \code{ggplot2} functions).
#'
#' @param resids Data table of model residuals for all brain regions
#' @export
#'
#' @seealso \code{\link[stats]{qqnorm}}

check.resid <- function(resids) {
  region <- x <- ysort <- NULL
  if (!'Group' %in% names(resids)) resids$Group <- 'Group 1'
  resids.melted <- melt(resids, id.vars='Group', variable.name='region',
                      value.name='resid')
  kNumRegions <- length(resids.melted[, levels(region)])

  a <- ifelse(kNumRegions < 9, 1, kNumRegions %/% 9)
  b <- kNumRegions %% 9

  resids.melted[, x := qnorm(ppoints(resid)), by=region]
  resids.melted[, ysort := sort(resid), by=region]
  setkey(resids.melted, region)

  ggQQ <- function(R) {
    p <- ggplot(R, aes(x=x, y=ysort)) +
      geom_point() +
      geom_line(aes(x=x, y=x)) +
      xlab('Theoretical Quantiles') + ylab('Sample Quantiles') +
      facet_wrap( ~ region, nrow=3, ncol=3, scales='free')
    return(p)
  }
    
  for (j in seq_len(a)) {
    dev.new()
    N1 <- 9 * (j - 1) + 1
    N2 <- min(N1 + 8, kNumRegions)
    print(ggQQ(resids.melted[levels(region)[N1:N2]]))
  }

  if (kNumRegions > 9 & ! b == 0) {
    dev.new()
    print(ggQQ(resids.melted[levels(region)[(N2+1):kNumRegions]]))
  }
}
