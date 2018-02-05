#' Calculate correlation matrix and threshold
#'
#' \code{corr.matrix} calculates the correlation between all column pairs of a
#' given data frame, and thresholds the resultant correlation matrix based on a
#' given density (e.g., \code{0.1} if you want to keep only the 10\% strongest
#' correlations). If you want to threshold by a specific correlation coefficient
#' (via the \code{thresholds} argument), then the \code{densities} argument is
#' ignored.
#'
#' If you wish to exclude regions from your analysis, you can give the indices
#' of their columns with the \code{exclusions} argument.
#'
#' By default, the Pearson correlation coefficients are calculated, but you can
#' return Spearman by passing an additional argument.
#'
#' @param resids Data table of the residuals (from \code{\link{get.resid}})
#' @param densities Numeric vector indicating the resultant network
#'   density(ies); keeps the top \emph{X}\% of correlations
#' @param thresholds Numeric; absolute correlation value to threshold by (default:
#'   \code{NULL})
#' @param exclusions Numeric vector of indices (columns) to exclude (default:
#'   \code{NULL})
#' @param ... Other arguments, passed to \code{\link[Hmisc]{rcorr}}
#' @export
#' @importFrom Hmisc rcorr
#'
#' @return A list with the following components:
#'   \item{R}{Numeric matrix of correlation coefficients.}
#'   \item{P}{Numeric matrix of p-values.}
#'   \item{r.thresh}{A 3-d binary array indicating correlations that are above a
#'     certain threshold. The length of the 3rd dimension equals the number of
#'     thresholds/densities supplied.}
#'   \item{thresholds}{Numeric vector; the thresholds supplied.}
#'   \item{densities}{Numeric vector; the densities supplied.}
#'
#' @family Structural covariance network functions
#' @seealso \code{\link[Hmisc]{rcorr}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' corrs <- lapply(groups, function(x) lapply(densities, function(y)
#'   corr.matrix(resids.all[x], densities=y)))
#' }

corr.matrix <- function(resids, densities, thresholds=NULL, exclusions=NULL, ...) {
  resids[, c('Study.ID', 'Group') := NULL]
  if (!is.null(exclusions)) resids <- resids[, -exclusions, with=F]

  corrs <- rcorr(as.matrix(resids), ...)
  r <- corrs$r
  p <- corrs$P

  # Calculate a threshold so that "density"% of possible connections are present
  if (hasArg('densities')) {
    N <- ncol(r)
    emax <- N  * (N - 1) / 2
    thresholds <- sort(r[lower.tri(r)])[emax - densities * emax]
  }
  r.thresh <- array(0, dim=c(dim(r), length(thresholds)),
                    dimnames=list(rownames(r), colnames(r)))
  for (i in seq_along(thresholds)) r.thresh[, , i] <- ifelse(r > thresholds[i], 1, 0)
  out <- list(R=r, P=p, r.thresh=r.thresh, thresholds=thresholds)
  if (hasArg('densities')) out <- c(out, list(densities=densities))

  return(out)
}
