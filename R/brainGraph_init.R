#' Initialize variables for further use in brainGraph
#'
#' This function initializes some variables that are important for further
#' analysis with the \code{brainGraph} package. This mostly involves loading CSV
#' files (of covariates/demographics, cortical thickness/volumes, etc.) and
#' returning them as data tables.
#'
#' The file containing covariates should be names \code{covars.csv}. The files
#' containing volumetric data should include hemisphere, atlas, and modality,
#' e.g. \code{lh_dkt_thickness.csv}. If you would like to include subcortical
#' gray matter, then you will need files \code{covars.scgm.csv} and
#' \code{scgm.csv}.
#'
#' @param atlas A character string indicating which brain atlas you are using
#' @param densities A numeric vector of the graph densities you would like to
#' investigate
#' @param datadir A character string; the filesystem location of your input
#' files
#' @param modality A character string indicating the volumetric MRI
#' modality/measure you are using to create the graphs ('thickness', 'volume',
#' 'lgi', or 'area')
#' @param use.mean A logical indicating whether or not you would like to
#' calculate the mean hemispheric volumetric measure (for later use in linear
#' models) (default: FALSE)
#' @param exclude.subs (optional) A character vector of the Study ID's of
#' subjects who are to be excluded from the analysis
#' @export
#'
#' @return A list containing:
#' \item{atlas}{A character string of the brain atlas name}
#' \item{densities}{A numeric vector of the graph densities}
#' \item{modality}{A character string of the modality you chose}
#' \item{kNumDensities}{An integer indicating the number of densities}
#' \item{covars}{A \code{data.table} of covariates}
#' \item{groups}{A character vector of subject group names}
#' \item{kNumGroups}{An integer indicating the number of groups}
#' \item{kNumVertices}{An integer; the number of vertices in the graphs}
#' \item{lhrh}{A \code{data.table} of left- and right-hemispheric volumetric
#' data}
#' \item{all.dat}{A merged \code{data.table} of \code{covars} and \code{lhrh}}
#' \item{all.dat.tidy}{A 'tidied' version of \code{all.dat}}
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' init.vars <- brainGraph_init(atlas='dkt', densities=seq(0.07, 0.50, 0.01),
#' datadir='/home/cwatson/Data', modality='thickness', exclude.subs=c('Con07',
#' 'Con23', 'Pat15'), use.mean=FALSE)
#' }

brainGraph_init <- function(atlas=c('aal116', 'aal90', 'brainsuite', 'destrieux',
                                    'dk', 'dk.scgm', 'dkt', 'dkt.scgm', 'hoa112',
                                    'lpba40'), densities, datadir,
                            modality=c('thickness', 'volume', 'lgi', 'area'),
                            use.mean=FALSE, exclude.subs=NULL) {

  Group <- Study.ID <- hemi <- name <- mean.lh <- mean.rh <- group.mean <- NULL
  value <- region <- NULL
  kNumDensities <- length(densities)
  atlas <- match.arg(atlas)
  atlas.dt <- eval(parse(text=atlas))
  kNumVertices <- nrow(atlas.dt)

  if (!file.exists(paste0(datadir, '/covars.csv'))) {
    stop(sprintf('File "covars.csv" does not exist in %s', datadir))
  }
  covars <- fread(paste0(datadir, '/covars.csv'))
  covars[, Group := as.factor(Group)]
  setkey(covars, Study.ID, Group)
  groups <- covars[, levels(Group)]
  kNumGroups <- length(groups)

  modality <- match.arg(modality)
  if (!file.exists(paste0(datadir, '/lh_', atlas, '_', modality, '.csv'))) {
    stop(sprintf('File "%s" does not exist in %s',
                 paste0(datadir, '/lh_', atlas, '_', modality, '.csv'),
                 datadir))
  }
  if (!file.exists(paste0(datadir, '/rh_', atlas, '_', modality, '.csv'))) {
    stop(sprintf('File "%s" does not exist in %s',
                 paste0(datadir, '/rh_', atlas, '_', modality, '.csv'),
                 datadir))
  }
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
  } else {
    covars <- subset(all.dat, select=names(covars))
  }
  lhrh <- subset(all.dat, select=names(lhrh))

  # Get SCGM and its covariates, if included
  if (isTRUE(grepl('scgm', atlas))) {
    scgm <- fread(paste0(datadir, '/scgm.csv'))
    setkey(scgm, Study.ID)
    covars.scgm <- fread(paste0(datadir, '/covars.scgm.csv'))
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
