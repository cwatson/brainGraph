#' Initialize variables for further use in brainGraph
#'
#' Initializes some variables that are important for further analysis of
#' volumetric (e.g., \emph{cortical thickness}) data. This mostly involves
#' loading CSV files (of covariates/demographics, cortical thickness/volumes,
#' etc.) and returning them as data tables.
#'
#' The file containing covariates should be named \code{covars.csv}. However,
#' you may also supply a \code{data.table} using the function argument
#' \code{covars}. This is useful if you have multiple covariates in your file
#' and wish to subset the data on your own.
#'
#' The filenames of files containing volumetric data should include hemisphere,
#' atlas, and modality separated by the \emph{underscore} character, e.g.
#' \code{lh_dkt_thickness.csv}. If you would like to include subcortical gray
#' matter, then you will need files \code{covars.scgm.csv} and \code{scgm.csv}.
#'
#' @param atlas Character string indicating which brain atlas you are using
#' @param densities Numeric vector of the graph densities you would like to
#'   investigate
#' @param datadir Character string; the filesystem location of your input files
#' @param modality Character string indicating the volumetric MRI
#'   modality/measure used to create the graphs (default: \code{thickness})
#' @param use.mean Logical indicating whether or not you would like to
#'   calculate the mean hemispheric volumetric measure (for later use in linear
#'   models) (default: \code{FALSE})
#' @param covars A \code{data.table} of covariates; specify this if you do not
#'   want to load your full covariates file (default: \code{NULL})
#' @param exclude.subs Character vector of the Study ID's of subjects who are to
#'   be excluded from the analysis (default: \code{NULL})
#' @export
#'
#' @return A list containing:
#' \item{atlas}{Character string of the brain atlas name}
#' \item{densities}{Numeric vector of the graph densities}
#' \item{modality}{Character string of the modality you chose}
#' \item{kNumDensities}{Integer indicating the number of densities}
#' \item{covars}{A \code{data.table} of covariates}
#' \item{groups}{Character vector of subject group names}
#' \item{kNumGroups}{Integer indicating the number of groups}
#' \item{kNumVertices}{Integer; the number of vertices in the graphs}
#' \item{lhrh}{A \code{data.table} of left- and right-hemispheric volumetric
#' data}
#' \item{all.dat}{A merged \code{data.table} of \code{covars} and \code{lhrh}}
#' \item{all.dat.tidy}{A "tidied" version of \code{all.dat}}
#'
#' @family Volumetric functions
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' init.vars <- brainGraph_init(atlas='dkt', densities=seq(0.07, 0.50, 0.01),
#'   datadir='/home/cwatson/Data', modality='thickness', exclude.subs=c('Con07',
#'   'Con23', 'Pat15'), use.mean=FALSE)
#' }

brainGraph_init <- function(atlas, densities, datadir,
                            modality=c('thickness', 'volume', 'lgi', 'area'),
                            use.mean=FALSE, covars=NULL, exclude.subs=NULL) {

  Group <- Study.ID <- hemi <- name <- mean.lh <- mean.rh <- group.mean <-
    value <- region <- NULL
  kNumDensities <- length(densities)
  atlas <- match.arg(atlas, choices=data(package='brainGraph')$results[, 3])
  atlas.dt <- eval(parse(text=atlas))
  kNumVertices <- nrow(atlas.dt)

  if (is.null(covars)) {
    stopifnot(file.exists(paste0(datadir, '/covars.csv')))
    covars <- fread(paste0(datadir, '/covars.csv'))
  }
  covars[, Group := as.factor(Group)]
  setkey(covars, Study.ID, Group)
  groups <- covars[, levels(Group)]
  kNumGroups <- length(groups)

  modality <- match.arg(modality)
  stopifnot(file.exists(paste0(datadir, '/lh_', atlas, '_', modality, '.csv')),
            file.exists(paste0(datadir, '/rh_', atlas, '_', modality, '.csv')))
  lh <- fread(paste0(datadir, '/lh_', atlas, '_', modality, '.csv'))
  setkey(lh, Study.ID)
  rh <- fread(paste0(datadir, '/rh_', atlas, '_', modality, '.csv'))
  setkey(rh, Study.ID)
  lhrh <- merge(lh, rh)

  # Remove subjects that are to be excluded
  if (!is.null(exclude.subs)) {
    covars <- covars[!Study.ID %in% exclude.subs]
    lhrh <- lhrh[!Study.ID %in% exclude.subs]
  }

  # Calculate hemispheric means, if desired
  all.dat <- merge(covars, lhrh)
  if (isTRUE(use.mean)) {
    all.dat[, mean.lh := rowMeans(.SD),
            .SDcols=atlas.dt[hemi == 'L', name],
            by=Study.ID]$V1
    all.dat[, mean.rh := rowMeans(.SD),
            .SDcols=atlas.dt[hemi == 'R', name],
            by=Study.ID]$V1
    covars <- subset(all.dat, select=c(names(covars), 'mean.lh', 'mean.rh'))
  }

  # Get SCGM and its covariates, if included
  if (isTRUE(grepl('scgm', atlas))) {
    scgm <- fread(paste0(datadir, '/scgm.csv'))
    setkey(scgm, Study.ID)
    if (is.null(covars)) {
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

    all.dat.scgm <- merge(covars.scgm, scgm)
    all.dat.scgm.tidy <- melt(all.dat.scgm, id.vars=names(covars.scgm),
                              variable.name='region')
    all.dat.scgm.tidy[, modality := modality]
    all.dat.scgm.tidy[, group.mean := mean(value), by=list(Group, region)]
    setkey(all.dat.scgm.tidy, Group, region)
  }

  all.dat.tidy <- melt(all.dat, id.vars=names(covars), variable.name='region')
  all.dat.tidy[, modality := modality]
  all.dat.tidy[, group.mean := mean(value), by=list(Group, region)]
  setkey(all.dat.tidy, Group, region)

  res <- list(atlas=atlas, densities=densities, modality=modality,
              kNumDensities=kNumDensities, covars=covars, groups=groups,
              kNumGroups=kNumGroups, kNumVertices=kNumVertices, lhrh=lhrh,
              all.dat=all.dat, all.dat.tidy=all.dat.tidy)
  if (isTRUE(grepl('scgm', atlas))) {
    res <- c(res, list(covars.scgm=covars.scgm, scgm=scgm,
                       all.dat.scgm=all.dat.scgm,
                       all.dat.scgm.tidy=all.dat.scgm.tidy))
  }
  return(res)
}
