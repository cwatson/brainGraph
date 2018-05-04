#' Initialize variables for further use in brainGraph
#'
#' Initializes some variables that are important for further analysis of
#' structural covariance networks (e.g., \emph{cortical thickness}). This mostly
#' involves loading CSV files (of covariates/demographics, cortical
#' thickness/volumes, etc.) and returning them as data tables.
#'
#' You can use any atlas that is already present in the package; to check the
#' available atlases, you can type
#' \code{data(package="brainGraph")$results[, 3]} at the R prompt. If you have a
#' custom atlas, specify \code{atlas="custom"} and supply the R object's name
#' for the argument \code{custom.atlas}.
#'
#' The file containing covariates should be named \code{covars.csv}. However,
#' you may also supply a \code{data.table} using the function argument
#' \code{covars}. This is useful if you have multiple covariates in your file
#' and wish to subset the data on your own. It must have, at minimum, columns
#' named \code{Study.ID} and \code{Group} (even if you have only 1 group).
#'
#' The filenames of the structural MRI data should include hemisphere,
#' atlas, and modality separated by the \emph{underscore} character; e.g.
#' \code{lh_dkt_thickness.csv} contains cortical thickness of left hemisphere
#' regions of the DKT atlas. If you would like to include subcortical gray
#' matter, then you will need files \code{covars.scgm.csv} and \code{scgm.csv}.
#'
#' @param atlas Character string indicating which brain atlas you are using.
#'   This can be an atlas present with the package or a custom atlas; in this
#'   case you must specify \code{custom} here and assign the name to the
#'   argument \code{custom.atlas}.
#' @param densities Numeric vector of the graph densities you would like to
#'   investigate
#' @param datadir Character string; the filesystem location of your input files
#' @param modality Character string indicating the volumetric MRI
#'   modality/measure used to create the graphs (default: \code{thickness})
#' @param covars A \code{data.table} of covariates; specify this if you do not
#'   want to load your full covariates file (default: \code{NULL})
#' @param exclude.subs Character vector of the Study ID's of subjects who are to
#'   be excluded from the analysis (default: \code{NULL})
#' @param custom.atlas Character string of the name of the custom atlas you wish
#'   to use, if applicable (default: \code{NULL})
#' @export
#'
#' @return A list containing:
#'   \item{atlas}{Character string of the brain atlas name}
#'   \item{densities}{Numeric vector of the graph densities}
#'   \item{modality}{Character string of the modality you chose}
#'   \item{kNumDensities}{Integer indicating the number of densities}
#'   \item{covars}{A \code{data.table} of covariates}
#'   \item{groups}{Character vector of subject group names}
#'   \item{kNumGroups}{Integer indicating the number of groups}
#'   \item{kNumVertices}{Integer; the number of vertices in the graphs}
#'   \item{lhrh}{A \code{data.table} of left- and right-hemispheric volumetric
#'     data}
#'
#' @family Structural covariance network functions
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' init.vars <- brainGraph_init(atlas='dkt', densities=seq(0.07, 0.50, 0.01),
#'   datadir='/home/cwatson/Data', modality='thickness', exclude.subs=c('Con07',
#'   'Con23', 'Pat15'))
#' }

brainGraph_init <- function(atlas, densities, datadir, modality='thickness',
                            covars=NULL, exclude.subs=NULL, custom.atlas=NULL) {
  Group <- Study.ID <- region <- NULL

  kNumDensities <- length(densities)
  atlas <- match.arg(atlas, choices=c(data(package='brainGraph')$results[, 3], 'custom'))
  if (atlas == 'custom') {
    stopifnot(!is.null(custom.atlas), exists(custom.atlas))
    atlas <- custom.atlas
  }

  if (is.null(covars)) {
    stopifnot(file.exists(paste0(datadir, '/covars.csv')))
    covars <- fread(paste0(datadir, '/covars.csv'))
  }
  covars[, Group := as.factor(Group)]
  setkey(covars, Study.ID, Group)
  groups <- covars[, levels(Group)]
  kNumGroups <- length(groups)

  stopifnot(file.exists(paste0(datadir, '/lh_', atlas, '_', modality, '.csv')),
            file.exists(paste0(datadir, '/rh_', atlas, '_', modality, '.csv')))
  lh <- fread(paste0(datadir, '/lh_', atlas, '_', modality, '.csv'))
  setkey(lh, Study.ID)
  rh <- fread(paste0(datadir, '/rh_', atlas, '_', modality, '.csv'))
  setkey(rh, Study.ID)
  lhrh <- merge(lh, rh)
  kNumVertices <- ncol(lhrh) - 1

  # Remove subjects that are to be excluded
  if (!is.null(exclude.subs)) {
    covars <- covars[!Study.ID %in% exclude.subs]
    lhrh <- lhrh[!Study.ID %in% exclude.subs]
  }

  res <- list(atlas=atlas, densities=densities, modality=modality,
              kNumDensities=kNumDensities, covars=covars, groups=groups,
              kNumGroups=kNumGroups, kNumVertices=kNumVertices, lhrh=lhrh)

  # Get SCGM and its covariates, if included
  if (isTRUE(grepl('scgm', atlas))) {
    scgm <- fread(paste0(datadir, '/scgm.csv'))
    setkey(scgm, Study.ID)
    if (is.null(covars)) { #FIXME: this branch will never be selected
      covars.scgm <- fread(paste0(datadir, '/covars.scgm.csv'))
    } else {
      covars.scgm <- fread(paste0(datadir, '/covars.csv'))
      covars.scgm <- covars.scgm[scgm == 1, !c('scgm', 'thickness', 'lgi.lh', 'lgi.rh', 'tract'), with=F]
    }
    covars.scgm[, Group := as.factor(Group)]
    setkey(covars.scgm, Study.ID, Group)

    if (!is.null(exclude.subs)) {
      covars.scgm <- covars.scgm[!Study.ID %in% exclude.subs]
      scgm <- scgm[!Study.ID %in% exclude.subs]
    }
    res <- c(res, list(covars.scgm=covars.scgm, scgm=scgm))
  }
  return(res)
}
