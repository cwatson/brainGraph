#' Calculate correlation matrix and threshold
#'
#' This function does a column-by-column correlation of a given data frame, and
#' will threshold the matrix based on a given density; e.g. 0.1 if you want to
#' keep only the 10\% strongest correlations.
#'
#' If you wish to exclude regions from your analysis, you can give the indices
#' of their columns. Also returns the p-values. Essentially a wrapper for
#' \code{\link[Hmisc]{rcorr}}, with some added functionality to work with this
#' type of data more easily. By default, the Pearson correlation coefficients
#' are calculated, but can return Spearman by passing an additional argument.
#'
#' @param dat Matrix or data frame of the columns to correlate
#' @param thresh Absolute correlation value to threshold by
#' @param density Keeps the top \emph{X}\% of correlations
#' @param exclusions Vector of indices (columns) to exclude (optional)
#' @param ... Other arguments to be passed to \code{\link[Hmisc]{rcorr}}
#' @export
#'
#' @return A list with the following components:
#' \item{R}{Correlation coefficients (default: Pearson).}
#' \item{P}{Associated p-values.}
#' \item{r.thresh}{Binary matrix indicating correlations that are above a
#' certain threshold.}
#' \item{threshold}{The threshold used.}
#' @seealso \code{\link[Hmisc]{rcorr}}

corr.matrix <- function(dat, thresh=NULL, density=0.1, exclusions=NULL, ...) {
  Group <- Study.ID <- NULL
  if ('Group' %in% names(dat)) dat[, Group := NULL]
  if ('Study.ID' %in% names(dat)) dat[, Study.ID := NULL]

  if (is.null(exclusions)) {
    corrs <- Hmisc::rcorr(as.matrix(dat), ...)
  } else {
    corrs <- Hmisc::rcorr(as.matrix(dat[, -exclusions, with=F]), ...)
  }
  r <- corrs$r
  p <- corrs$P

  # Calculate a threshold so that "density"% of possible connections are present
  if (hasArg('density')) {
    N <- ncol(r)
    emax <- N  * (N - 1) / 2

    thresh <- sort(r[lower.tri(r)])[emax - density * emax]
  }
  r.thresh <- ifelse(r > thresh, 1, 0)
  out <- list(R=r, P=p, r.thresh=r.thresh, threshold=thresh)
}
