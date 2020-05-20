#' Atlas helper functions
#'
#' \code{guess_atlas} tries to determine which atlas is being used based on the
#' data; i.e., the number of vertices/regions.
#'
#' @section Guessing the atlas from an object:
#' There are several valid inputs to \code{guess_atlas}:
#' \describe{
#'   \item{data.table}{The atlas will be guessed based on the number of columns
#'     (subtracting by 1 if a \dQuote{Study ID} column is present). This is the
#'     same behavior as for \code{data.frame} objects, as well.}
#'   \item{igraph}{The vertex count}
#'   \item{brainGraph}{If there is a \code{atlas} graph-level attribute, it will
#'     return that. Otherwise, the vertex count.}
#'   \item{matrix,array}{The number of rows, which should equal the number of
#'     columns if the input is a connectivity matrix.}
#' }
#'
#' Note that this will only work properly for atlases that are currently in the
#' package. If you are using a custom atlas and you receive errors, please open
#' an issue on \emph{GitHub}.
#'
#' @param x,object An object to test or convert to an atlas data.table
#' @export
#' @return \code{guess_atlas} - Character string; either the matched atlas or
#'   \code{NA}
#' @name Atlas Helpers
#' @rdname atlas_helpers

guess_atlas <- function(x) {
  bgAtlases <- data(package='brainGraph')$results[, 3L]
  Nv <- vapply(bgAtlases, function(x) dim(get(x))[1L], integer(1L))

  n <- switch(class(x)[1L],
              data.frame=, data.table=dim(x)[2L] - (hasName(x, getOption('bg.subject_id'))),
              igraph=vcount(x),
              brainGraph=if (is.null(x$atlas)) vcount(x) else x$atlas,
              matrix=, array=dim(x)[1L],
              x)

  if (is.character(n)) {
    atlas <- n
  } else {
    matched <- which(Nv == n)
    atlas <- if (length(matched) > 0L) names(matched) else NA_character_
  }
  return(atlas)
}

#' Check if an object is a valid atlas, and convert if it isn't
#'
#' \code{as_atlas} and \code{create_atlas} converts/coerces an object to a
#' a \code{data.table}, or creates one, that is compatible with
#' \code{brainGraph}.
#'
#' @section Coercing to an atlas:
#' There are several things \code{as_atlas} tries to do to make it work without
#' error:
#' \itemize{
#'   \item Coerce the object to \code{data.table}
#'   \item Add a column of integers named \code{index}
#'   \item Change columns named \code{'x'}, \code{'y'}, or \code{'z'} to have
#'     \code{.mni} at the end
#' }
#'
#' @export
#' @return \code{as_atlas} and \code{create_atlas} return a \code{data.table}
#'   that conforms to other atlases in the package, or exits with an error.
#' @rdname atlas_helpers
#' @examples
#' my_atlas <- data.frame(name=paste('Region', 1:10), x.mni=rnorm(10),
#'   y.mni=rnorm(10), z.mni=rnorm(10),
#'   lobe=rep(c('Frontal', 'Parietal', 'Temporal', 'Occipital', 'Limbic'), 2),
#'   hemi=c(rep('L', 5), rep('R', 5)))
#' my_atlas2 <- as_atlas(my_atlas)
#' str(my_atlas)
#' str(my_atlas2)

as_atlas <- function(object) {
  index <- NULL
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
  if (!all(hasName(object, atlas_cols))) {
    stop(paste('Atlas is invalid. Required columns are:\n', paste(atlas_cols, collapse=' ')))
  }

  for (col in c('lobe', 'hemi', 'name')) {
    if (!is.factor(object[[col]])) {
      warning('Coercing "', col, '" column to factor')
      object[, eval(col) := as.factor(get(col))]
    }
  }
  other_names <- setdiff(atlas_cols, names(object))
  setcolorder(object, c(atlas_cols, other_names))
  return(object)
}

#' @param regions Character vector of region names
#' @param coords Numeric matrix of spatial coordinates; must have 3 columns
#' @param lobes Character or factor vector of lobe membership
#' @param hemis Character or factor vector of hemisphere membership. There
#'   should probably not be more than 3 unique elements (for left, right, and
#'   bi-hemispheric regions)
#' @export
#' @rdname atlas_helpers
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
    coords <- coords[[1L]]
  }
  if (!is.numeric(coords)) {
    warning('Coercing coordinates to numeric.')
    try(coords <- as.numeric(coords))
  }
  if (!is.matrix(coords)) {
    warning('Converting coordinates to a 3-column matrix.')
    try(coords <- matrix(coords, ncol=3L))
  }
  dims <- dim(coords)
  if (dims[2L] != 3L) stop('Coordinates must be a matrix with 3 columns.')

  kNumRegions <- length(regions)
  kNumCoords <- dims[1L]
  kNumLobes <- length(lobes)
  kNumHemis <- length(hemis)
  if (length(unique(c(kNumRegions, kNumCoords, kNumLobes, kNumHemis))) != 1L) {
    stop('All input data must have the same length.')
  }

  out <- data.table(name=as.character(regions),
                    x.mni=coords[, 1L],
                    y.mni=coords[, 2L],
                    z.mni=coords[, 3L],
                    lobe=as.factor(lobes),
                    hemi=as.factor(hemis),
                    index=seq_len(kNumCoords))
  return(out)
}
