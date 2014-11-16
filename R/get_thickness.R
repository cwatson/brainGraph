#' Load the thickness data and separate into 2 groups.
#'
#' This function loads the thickness data from csv files (separated into left
#' and right hemispheres. It combines data from each hemisphere, then subsets
#' the data frame into separate DF's for each of 2 subject groups.
#'
#' @param fileLH Character string for the file containing left hemisphere
#' thickness data for both groups
#' @param fileRH Character string for the file containing right hemisphere
#' thickness data for both groups
#' @param group1 Character string identifying first subject group
#' @param group2 Character string identifying second subject group (optional)
#' @param covars Data frame of covariates, one of which indicates group
#' assignment, if necessary
#' @export
#'
#' @return A list with components:
#' \item{lh}{Thicknesses of the left hemisphere.}
#' \item{rh}{Thicknesses of the right hemisphere.}
#' \item{all}{Thicknesses of both hemispheres.}
#' \item{group1}{Thicknesses for only group 1.}
#' \item{group2}{Thicknesses for only group 2.}

get.thickness <- function(fileLH, fileRH, group1, group2=NULL, covars=NULL) {
  lhThick <- read.csv(fileLH)
  rhThick <- read.csv(fileRH)
  all.thick <- merge(lhThick, rhThick)
  if (! 'Group' %in% colnames(all.thick)) {
    all.thick.cov <- merge(all.thick, covars)
  } else {
    all.thick.cov <- all.thick
  }

  group1.thick <- subset(all.thick.cov, Group==group1)
  if (exists('group2')) {
    group2.thick <- subset(all.thick.cov, Group==group2)
    all.thick <- subset(all.thick, select=-Group)
    list(lh=lhThick, rh=rhThick, all=all.thick, group1=group1.thick,
      group2=group2.thick)
  } else {
    all.thick <- subset(all.thick, select=-Group)
    list(lh=lhThick, rh=rhThick, all=all.thick, group1=group1.thick)
  }
}
