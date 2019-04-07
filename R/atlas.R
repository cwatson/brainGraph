#' Guess the atlas based on the data
#'
#' If \code{atlas} is not explicitly provided, \code{guess_atlas} will try to
#' determine which atlas is being used based on the data, particularly the
#' number of vertices/regions.
#'
#' @param x An object to test
#' @export
#' @return Character string; either the matched atlas or \code{NA}
#' @name AtlasHelpers
#' @rdname atlas

guess_atlas <- function(x) {
  bgAtlases <- data(package='brainGraph')$results[, 3]
  Nv <- sapply(bgAtlases, function(x) nrow(get(x)))

  # Currently works only for 'dt.vol' from `get.resid`
  if (is.data.table(x)) {
    n <- ncol(x)
    if ('Study.ID' %in% names(x)) n <- n - 1
  } else if (is_igraph(x)) {
    n <- vcount(x)
  } else {
    n <- x
  }

  matched <- which(Nv == n)
  if (length(matched > 0)) {
    atlas <- names(matched)
  } else {
    atlas <- NA_character_
  }
  return(atlas)
}

#' Check if an object is a valid atlas, and convert if it isn't
#'
#' \code{as_atlas} will convert an object to a \code{data.table} that is
#' compatible with \code{brainGraph}.
#'
#' There are several things this function tries to do to make it work without
#' error:
#' \itemize{
#'   \item Coerce the object to \code{data.table}
#'   \item Add a column of integers named \code{index}
#'   \item Change columns named \code{'x'}, \code{'y'}, or \code{'z'} to have
#'     \code{.mni} at the end
#' }
#'
#' @param object An object to test
#' @export
#' @return If possible, it returns a \code{data.table} that conforms to other
#'   atlases in the package. Otherwise, exits with an error.
#' @rdname atlas
#' @examples
#' my_atlas <- data.frame(name=paste('Region', 1:10), x.mni=rnorm(10),
#'   y.mni=rnorm(10), z.mni=rnorm(10),
#'   lobe=rep(c('Frontal', 'Parietal', 'Temporal', 'Occipital', 'Limbic'), 2),
#'   hemi=c(rep('L', 5), rep('R', 5)))
#' my_atlas <- as_atlas(my_atlas)
#' str(my_atlas)

as_atlas <- function(object) {
  # First, coerce to data.table
  if (!is.data.table(object)) {
    warning('Atlas is being coerced to "data.table"')
    object <- as.data.table(object)
  }

  # Check column names; add "index" if not present
  xnames <- names(object)
  if (!'index' %in% xnames) object[, index := seq_len(.N)]
  if ('x' %in% xnames) setnames(object, 'x', 'x.mni')
  if ('y' %in% xnames) setnames(object, 'y', 'y.mni')
  if ('z' %in% xnames) setnames(object, 'z', 'z.mni')
  atlas_cols <- c('name', 'x.mni', 'y.mni', 'z.mni', 'lobe', 'hemi')
  if (!all(atlas_cols %in% names(object))) {
    stop(paste('Atlas is invalid. Required columns are:\n', paste(atlas_cols, collapse=' ')))
  }

  if (!is.factor(object$lobe)) {
    warning('Coercing "lobe" column to factor')
    object[, lobe := as.factor(lobe)]
  }
  if (!is.factor(object$hemi)) {
    warning('Coercing "hemi" column to factor')
    object[, hemi := as.factor(hemi)]
  }
  if (!is.character(object$name)) {
    warning('Coercing "name" column to character')
    object[, name := as.character(name)]
  }
  other_names <- setdiff(atlas_cols, names(object))
  setcolorder(object, c(atlas_cols, other_names))
  return(object)
}

#' Create an atlas compatible with brainGraph
#'
#' \code{create_atlas} will create a \code{data.table} object that is compatible
#' with \code{brainGraph}. Required pieces of information are: region names,
#' spatial coordinates, and lobe and hemisphere membership
#'
#' @param regions Character vector of region names
#' @param coords Numeric matrix of spatial coordinates; must have 3 columns
#' @param lobes Character or factor vector of lobe membership
#' @param hemis Character or factor vector of hemisphere membership. There
#'   should probably not be more than 3 unique elements (for left, right, and
#'   bi-hemispheric regions)
#' @export
#' @return A \code{data.table} that is compatible with \code{brainGraph} atlases
#' @rdname atlas
#' @examples
#' regions <- paste('Region', 1:10)
#' xyz <- matrix(rnorm(30), nrow=10, ncol=3)
#' lobe <- rep(c('Frontal', 'Parietal', 'Temporal', 'Occipital', 'Limbic'), 2)
#' hemi <- c(rep('L', 5), rep('R', 5))
#' my_atlas <- create_atlas(regions, xyz, lobe, hemi)
#' str(my_atlas)

create_atlas <- function(regions, coords, lobes, hemis) {
  if (is.list(coords)) {
    warning("'coords' is a list; using only the first element.")
    coords <- coords[[1]]
  }
  if (!is.numeric(coords)) {
    warning('Coercing coordinates to numeric.')
    try(coords <- as.numeric(coords))
  }
  if (!is.matrix(coords)) {
    warning('Converting coordinates to a 3-column matrix.')
    try(coords <- matrix(coords, ncol=3))
  }
  if (ncol(coords) != 3) stop('Coordinates must be a matrix with 3 columns.')

  kNumRegions <- length(regions)
  kNumCoords <- nrow(coords)
  kNumLobes <- length(lobes)
  kNumHemis <- length(hemis)
  if (length(unique(c(kNumRegions, kNumCoords, kNumLobes, kNumHemis))) != 1) {
    stop('All input data must have the same length.')
  }

  out <- data.table(name=as.character(regions),
                    x.mni=coords[, 1],
                    y.mni=coords[, 2],
                    z.mni=coords[, 3],
                    lobe=as.factor(lobes),
                    hemi=as.factor(hemis),
                    index=seq_len(nrow(coords)))
  return(out)
}
