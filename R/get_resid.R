#' Linear model residuals across brain regions
#'
#' This function runs linear models across brain regions listed in a
#' \code{data.table} (e.g. cortical thickness), in order to adjust for relevant
#' variables (e.g. age, sex, group, etc.). It adds the \emph{studentized}
#' residuals as a column and returns the data table.
#'
#' @param tidy.dt A \code{data.table} that has been "tidied", containing all
#' covariates and the brain measure of interest
#' @param covars A \code{data.table} of covariates
#' @param use.mean A logical indicating whether to control for the mean
#' hemispheric brain value (e.g. mean LH/RH cortical thickness) (default: NULL)
#' @param exclude A character vector of columns to exclude (default: NULL)
#' @export
#'
#' @return A list with components:
#' \item{all.dat.tidy}{The tidied \code{data.table} with 'resids' column added}
#' \item{formulas}{Character string of the \code{\link{lm}} formulas used}
#' \item{resids.all}{The "wide" \code{data.table} of residuals}
#' @seealso \code{\link{rstudent}}

get.resid <- function(tidy.dt, covars, use.mean=FALSE, exclude=NULL) {
  region <- resids <- Group <- Study.ID <- NULL
  myDT <- copy(tidy.dt)
  exclude <- c('Study.ID', exclude)

  # Adjust by mean hemispheric values
  if (isTRUE(use.mean)) {
    formula.lh <- paste('value ~',
                        paste(names(covars[, !c(exclude, 'mean.rh'), with=F]),
                                         collapse='+'))
    formula.rh <- paste('value ~',
                        paste(names(covars[, !c(exclude, 'mean.lh'), with=F]),
                                         collapse='+'))
    if (nrow(myDT[grep('^l.*', region)]) == 0) {
      lh.string <- '.*\\.L$'
      rh.string <- '.*\\.R$'
    } else {
      lh.string <- '^l.*'
      rh.string <- '^r.*'
    }
    myDT[grep(lh.string, region),
            resids := rstudent(lm(as.formula(formula.lh), .SD)), by=region]
    myDT[grep(rh.string, region),
            resids := rstudent(lm(as.formula(formula.rh), .SD)), by=region]
    formulas <- c(formula.lh, formula.rh)

  # Don't adjust by mean hemispheric values
  } else {
    # In case 'mean.lh' or 'mean.rh' are in 'covars', exclude
    if ('mean.lh' %in% names(covars)) {
      exclude <- c(exclude, 'mean.lh')
    }
    if ('mean.rh' %in% names(covars)) {
      exclude <- c(exclude, 'mean.rh')
    }
    formula.lhrh <- paste('value ~ ', paste(names(covars[, !exclude, with=F]),
                                            collapse='+'))
    myDT[, resids := rstudent(lm(as.formula(formula.lhrh), .SD)), by=region]
    formulas <- formula.lhrh
  }

  # Return data to "wide" format with just the residuals
  all.dat.wide <- dcast.data.table(myDT, paste(paste(names(covars),
                                                     collapse='+'),
                                               '~ region'),
                                   value.var='resids')
  resids.all <- all.dat.wide[, !setdiff(names(covars), c('Study.ID', 'Group')), with=F]
  setkey(resids.all, Group, Study.ID)

  return(list(all.dat.tidy=myDT, formulas=formulas, resids.all=resids.all))
}
