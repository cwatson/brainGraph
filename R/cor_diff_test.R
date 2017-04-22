#' Calculate the p-value for differences in correlation coefficients
#'
#' Given two sets of correlation coefficients and sample sizes, this function
#' calculates and returns the \emph{z-scores} and \emph{p-values} associated
#' with the difference between correlation coefficients. This function was
#' adapted from \url{http://stackoverflow.com/a/14519007/3357706}.
#'
#' @param r1 Numeric (vector or matrix) of correlation coefficients, group 1
#' @param r2 Numeric (vector or matrix) of correlation coefficients, group 2
#' @param n1 Integer; number of observations, group 1
#' @param n2 Integer; number of observations, group 2
#' @param alternative Character string specifying the alternative hypothesis
#'   test to use; one of: 'two.sided' (default), 'less', 'greater'
#' @export
#'
#' @return A list containing:
#' \item{p}{The p-values}
#' \item{z}{The z-score for the difference in correlation coefficients}
#'
#' @family Matrix functions
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' kNumSubjs <- summary(covars$Group)
#' corr.diffs <- cor.diff.test(corrs[[1]][[1]]$R, corrs[[2]][[1]]$R,
#'                             kNumSubjs[1], kNumSubjs[2], alternative='two.sided')
#' edge.diffs <- t(sapply(which(corr.diffs$p < .05), function(x)
#'                        mapply('[[',
#'                               dimnames(corr.diffs$p),
#'                               arrayInd(x, dim(corr.diffs$p)))
#'                               ))
#' }
cor.diff.test <- function(r1, r2, n1, n2,
                          alternative = c('two.sided', 'less', 'greater')) {

  z1 <- 0.5 * log((1 + r1) / (1 - r1))
  z2 <- 0.5 * log((1 + r2) / (1 - r2))

  SEdiff <- sqrt((1 / (n1 - 3)) + (1 / (n2 - 3)))
  diff.z <- (z1 - z2) / SEdiff

  alt <- match.arg(alternative)
  if (alt == 'less') {
    p <- pnorm(diff.z)
  } else if (alt == 'greater') {
    p <- pnorm(diff.z, lower.tail=F)
  } else if (alt == 'two.sided') {
    p <- 2 * pnorm(abs(diff.z), lower.tail=F)
  }

  return(list(p=p, z=diff.z))
}
