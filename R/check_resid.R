#' Check model residuals for each brain region
#'
#' This function checks the model residuals for each brain region in the
#' analysis. It simply does a qqplot of the studentized residuals (but uses
#' \code{ggplot2} functions).
#'
#' @param resids Data table of model residuals for all brain regions
#' @param cols Logical indicating whether to color by group (default: FALSE)
#' @export
#'
#' @seealso \code{\link[stats]{qqnorm}}

check.resid <- function(resids, cols=FALSE) {
  region <- x <- ysort <- ind <- mark <- NULL
  if (!'Group' %in% names(resids)) resids$Group <- 'Group 1'
  resids.m <- melt(resids, id.vars='Group', variable.name='region',
                      value.name='resid')
  kNumRegions <- length(resids.m[, levels(region)])

  a <- ifelse(kNumRegions < 9, 1, kNumRegions %/% 9)
  b <- kNumRegions %% 9

  resids.m[, ind := as.character(.SD[, .I]), by=region]
  setkey(resids.m, region, resid)
  resids.m[, x := qnorm(ppoints(resid)), by=region]
  resids.m[, mark := ifelse(abs(resid) < mean(resid) + 2*sd(resid), 0, 1), by=region]
  resids.m[, mark := as.factor(mark)]
  resids.m[mark == 0, ind := '']

  ggQQ <- function(R, cols) {
    Group <- NULL
    if (isTRUE(cols)) {
      p <- ggplot(R, aes(x=x, y=resid)) +
      geom_text(aes(x=x, y=resid+0.3, label=ind, col=Group), size=3) +
      geom_line(aes(x=x, y=x), col='gray50') +
      geom_point(aes(col=Group, shape=mark, size=mark)) +
      xlab('Theoretical Quantiles') + ylab('Sample Quantiles') +
      facet_wrap( ~ region, nrow=3, ncol=3, scales='free') +
      scale_shape_manual(values=c(20, 8)) +
      scale_size_manual(values=c(1.5, 3)) +
      theme(legend.position='none')
    } else {
      p <- ggplot(R, aes(x=x, y=resid)) +
        geom_point() +
        geom_line(aes(x=x, y=x)) +
        xlab('Theoretical Quantiles') + ylab('Sample Quantiles') +
        facet_wrap( ~ region, nrow=3, ncol=3, scales='free')
    }
    return(p)
  }

  for (j in seq_len(a)) {
    dev.new()
    N1 <- 9 * (j - 1) + 1
    N2 <- min(N1 + 8, kNumRegions)
    print(ggQQ(resids.m[levels(region)[N1:N2]], cols))
  }

  if (kNumRegions > 9 & ! b == 0) {
    dev.new()
    print(ggQQ(resids.m[levels(region)[(N2+1):kNumRegions]], cols))
  }
}
