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
#' of their columns with the \code{exclude.reg} argument.
#'
#' By default, the Pearson correlation coefficients are calculated, but you can
#' return Spearman by changing the \code{type} argument.
#'
#' @param resids An object of class \code{brainGraph_resids} (the output from
#'   \code{\link{get.resid}})
#' @param densities Numeric vector indicating the resultant network
#'   density(ies); keeps the top \emph{X}\% of correlations
#' @param thresholds Numeric; absolute correlation value to threshold by
#'   (default: \code{NULL})
#' @param what Character string indicating whether to correlate the residuals or
#'   the raw structural MRI values (default: \code{'resids'})
#' @param exclude.reg Character vector of regions to exclude (default:
#'   \code{NULL})
#' @param type Character string indicating which type of correlation coefficient
#'   to calculate (default: \code{'pearson'})
#' @param rand Logical indicating whether the function is being called for
#'   permutation testing; not intended for general use (default: \code{FALSE})
#' @export
#' @importFrom Hmisc rcorr
#'
#' @return A nested list containing a list for all subject groups; each of these
#'   has the following components:
#'   \item{R}{Numeric matrix of correlation coefficients.}
#'   \item{P}{Numeric matrix of p-values.}
#'   \item{r.thresh}{A 3-d binary array indicating correlations that are above a
#'     certain threshold. The length of the 3rd dimension equals the number of
#'     thresholds/densities supplied.}
#'   \item{thresholds}{Numeric vector; the thresholds supplied.}
#'   \item{densities}{Numeric vector; the densities supplied.}
#'   \item{what}{Residuals or raw values}
#'   \item{exclude.reg}{Excluded regions (if any)}
#'   \item{type}{Pearson or Spearman}
#'
#' @family Structural covariance network functions
#' @seealso \code{\link[Hmisc]{rcorr}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' myResids <- get.resid(lhrh, covars)
#' corrs <- corr.matrix(myResids, densities=densities)))
#' }

corr.matrix <- function(resids, densities, thresholds=NULL, what=c('resids', 'raw'),
                        exclude.reg=NULL, type=c('pearson', 'spearman'), rand=FALSE) {
  Group <- Study.ID <- NULL
  what <- match.arg(what)
  type <- match.arg(type)
  stopifnot(inherits(resids, 'brainGraph_resids'))

  # Different behavior if called for permutation testing
  if (isTRUE(rand)) {
    res.all <- as.matrix(resids$resids.all[, !c('Study.ID', 'Group')])
    corrs <- rcorr(res.all)
    r <- corrs$r
    N <- ncol(r)
    emax <- N  * (N - 1) / 2
    thresholds <- sort(r[lower.tri(r)])[emax - densities * emax]
    r.thresh <- array(0, dim=c(dim(r), length(thresholds)),
                      dimnames=list(rownames(r), colnames(r)))
    for (i in seq_along(thresholds)) r.thresh[, , i] <- ifelse(r > thresholds[i], 1, 0)
    return(list(list(R=r, r.thresh=r.thresh)))
  }

  groups <- resids$groups
  out <- sapply(groups, function(x) NULL)
  if (what == 'resids') {
    res.all <- resids$resids.all[, !'Study.ID']
  } else if (what == 'raw') {
    res.all <- dcast(resids$all.dat.tidy, 'Study.ID + Group ~ region')
    setkey(res.all, Group, Study.ID)
    res.all <- res.all[, !'Study.ID']
  }
  if (!is.null(exclude.reg)) res.all <- res.all[, -exclude.reg, with=FALSE]

  for (g in groups) {
    corrs <- rcorr(as.matrix(res.all[g, !'Group']), type=type)
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
    out[[g]] <- list(R=r, P=p, r.thresh=r.thresh, thresholds=thresholds,
                     what=what, exclude.reg=exclude.reg, type=type)
    if (hasArg('densities')) out[[g]] <- c(out[[g]], list(densities=densities))
  }
  return(out)
}

#plot.brainGraph_corrs <- function(x, type=c('raw', 'thresholded'), ordered=TRUE, order.by='lobe', g=NULL, ...) {
#  return(p)
#}
