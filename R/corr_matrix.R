#' Calculate correlation matrix and threshold
#'
#' This function does a column-by-column correlation of a given data frame, and
#' will threshold r-values based on a given density value; e.g. 0.1 if you want
#' to keep only the 10\% strongest correlations. It also allows for the exclusion
#' of a set of columns (i.e. regions or nodes), given their indices. Also
#' returns the p-values. Essentially a wrapper for \code{\link[Hmisc]{rcorr}},
#' with some added functionality to work with the type of data more easily.
#'
#' @param dat Matrix or data frame of the columns to correlate
#' @param thresh Absolute correlation value to threshold by
#' @param density Keeps the top X\% of correlations
#' @param exclusions Vector of indices (columns) to exclude (optional)
#' @export
#'
#' @return A list with the following components:
#' \item{R}{Pearson correlation coefficients.}
#' \item{P}{Associated p-values.}
#' \item{r.thresh}{Binary matrix indicating correlations that are above a
#' certain threshold.}
#' \item{threshold}{The threshold used.}
#' @seealso \code{\link[Hmisc]{rcorr}}

corr.matrix <- function(dat, thresh=NULL, density=0.1, exclusions=NULL) {
  if (length(exclusions) == 0) {
    corrs <- rcorr(as.matrix(dat))
  } else {
    corrs <- rcorr(as.matrix(dat[, -exclusions]))
  }
  r <- corrs$r
  p <- corrs$P

  # Calculate a threshold so that 10% of possible connections are present
  # Also called 'sparsity'
  if (hasArg('density')) {
    N <- ncol(r)
    emax <- N  * (N - 1) / 2
  
    thresh <- sort(r[lower.tri(r)])[emax - density * emax]
  }
  r.thresh <- ifelse(abs(r) > thresh, 1, 0)
  out <- list(R=r, P=p, r.thresh=r.thresh, threshold=thresh)
}
