#' Run linear models of thickness data to get model residuals.
#'
#' This function runs linear models on the thickness data for each region, in
#' order to adjust for relevant variables (e.g. age, sex, group, etc.).
#'
#' @param thicknesses data frame of all thickness data, including covariates.
#'  First column must be subject ID (or something similar)
#' @param kNumCovars Integer specifying number of covariates to include in the
#' model
#' @param group1 Character string indicating the first subject group
#' @param group2 Character string indicating the second subject group (optional)
#' @export
#'
#' @return A list with components:
#' \item{models}{The \link{lm} objects for each brain region.}
#' \item{all}{All residuals for all subjects and brain regions.}
#' \item{group1}{Residuals for just group 1.}
#' \item{group2}{Residuals for just group 2.}

get.resid <- function(thicknesses, kNumCovars, group1, group2=NULL) {
  regions <- (kNumCovars+2):dim(thicknesses)[2]
  m <- lapply(names(thicknesses)[regions],
              function(x) lm(as.formula(paste(x, '~',
                              paste(colnames(thicknesses)[2:(kNumCovars+1)],
                                    collapse='+'), sep='')),
                             data=thicknesses))

  # Get the studentized residuals
  all.resid <- data.frame(lapply(m, rstudent))
  colnames(all.resid) <- colnames(thicknesses)[regions]
  all.resid$Group <- thicknesses$Group
  group1.resid <- subset(all.resid, Group==group1)

  if (exists('group2')) {
    group2.resid <- subset(all.resid, Group==group2)
    list(models=m, all=all.resid, group1=group1.resid, group2=group2.resid)
  } else {
    list(models=m, all=all.resid, group1=group1.resid)
  }
}
