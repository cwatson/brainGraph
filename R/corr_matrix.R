#' Calculate correlation matrix and threshold
#'
#' \code{corr.matrix} calculates the correlation between all column pairs of a
#' given data frame, and thresholds the resultant correlation matrix based on a
#' given density (e.g., \code{0.1} if you want to keep only the 10\% strongest
#' correlations.
#'
#' If you wish to exclude regions from your analysis, you can give the indices
#' of their columns with the \code{exclusions} argument.
#'
#' By default, the Pearson correlation coefficients are calculated, but can
#' return Spearman by passing an additional argument.
#'
#' @param dat Data table of the data to correlate
#' @param density Numeric indicating the resultant network density; keeps the
#'   top \emph{X}\% of correlations
#' @param thresh Numeric; absolute correlation value to threshold by (default:
#'   \code{NULL})
#' @param exclusions Numeric vector of indices (columns) to exclude (default:
#'   \code{NULL})
#' @param ... Other arguments, passed to \code{\link[Hmisc]{rcorr}}
#' @export
#'
#' @return A list with the following components:
#' \item{R}{Numeric matrix of correlation coefficients.}
#' \item{P}{Numeric matrix of p-values.}
#' \item{r.thresh}{Binary matrix indicating correlations that are above a
#'   certain threshold.}
#' \item{threshold}{Numeric; the threshold value used.}
#'
#' @family Volumetric functions
#' @seealso \code{\link[Hmisc]{rcorr}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' corrs <- lapply(groups, function(x) lapply(densities, function(y)
#'   corr.matrix(resids.all[x], density=y)))
#' }

corr.matrix <- function(dat, density, thresh=NULL, exclusions=NULL, ...) {
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
