#' Run linear models of thickness data to get model residuals.
#'
#' This function runs linear models on the thickness data for each region, in
#' order to adjust for relevant variables (e.g. age, sex, group, etc.).
#'
#' @param thicknesses Data frame of all thickness data. First column must be
#' subject ID (or something similar)
#' @param covars Data frame of covariates for the linear model
#' @export
#'
#' @return A list with components:
#' \item{models}{The \link{lm} objects for each brain region.}
#' \item{resids}{Data table of residuals for all subjects and brain regions.}

get.resid <- function(thicknesses, covars) {
  dat <- merge(covars, thicknesses)
  regions <- 2:ncol(thicknesses)
  m <- lapply(names(thicknesses)[regions],
              function(x) lm(as.formula(paste0(x, '~',
                              paste(names(covars)[2:ncol(covars)],
                                    collapse='+'))),
                             data=dat))

  # Get the studentized residuals
  all.resid <- as.data.table(lapply(m, rstudent))
  setnames(all.resid, colnames(thicknesses)[regions])
  if ('Group' %in% names(covars)) {
    all.resid$Group <- dat$Group
    setkey(all.resid, Group)
  }

  list(models=m, resids=all.resid)
}
