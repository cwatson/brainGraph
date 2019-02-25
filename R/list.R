#' Create a list of brainGraph graphs
#'
#' \code{make_bgList} creates a \code{brainGraphList} object, a list containing
#' a set of graphs for all subjects in a study at a specific threshold (or
#' density), in addition to some graph-level attributes common to those graphs.
#' The elements are:
#' \itemize{
#'   \item \code{threshold}
#'   \item \code{version} -- the version of \code{brainGraph} used when creating
#'     the graphs
#'   \item \code{atlas} -- the atlas common to all the graphs
#'   \item \code{modality} -- the imaging modality
#'   \item \code{weighting} -- a string indicating what edge weights represent
#'     (if applicable)
#'   \item \code{graphs} -- the \emph{list} of graphs
#' }
#'
#' This object can be considered comparable to a 4-D \emph{NIfTI} file,
#' particularly those returned by FSL's \emph{TBSS} prestats step since this
#' contains the FA volumes for all study subjects..
#'
#' @param A 3-D numeric array of all subjects connectivity matrices (for a
#'   single threshold)
#' @param threshold Numeric; the threshold used to create the connectivity
#'   matrices
#' @param atlas Character string; the brain atlas used
#' @param study.ids Character vector of study ID's (default: \code{NULL})
#' @param groups Character (or factor) vector of group names (default:
#'   \code{NULL})
#' @param ... Other arguments passed to \code{\link{set_brainGraph_attr}} (and
#'   subsequently to \code{\link{make_brainGraph}})
#' @export
#'
#' @return A list with elements
#'
#' @name brainGraphList
#' @rdname brainGraphList
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

make_brainGraph_list <- function(A, threshold, atlas, study.ids=NULL, groups=NULL, ...) {
  num.subjects <- dim(A)[3]
  if (!is.null(study.ids)) stopifnot(length(study.ids) == num.subjects)
  if (!is.null(groups)) stopifnot(length(groups) == num.subjects)

  out <- list(threshold=threshold, version=packageVersion('brainGraph'), atlas=atlas)

  # Get any extra arguments
  fargs <- list(...)
  if (hasArg('rand')) out$rand <- fargs$rand
  if (hasArg('modality')) out$modality <- fargs$modality
  if (hasArg('weighting')) out$weighting <- fargs$weighting

  if (is.null(study.ids)) study.ids <- seq_len(num.subjects)
  if (is.null(groups)) {
    groups <- rep(1, num.subjects)
  } else {
    if (is.factor(groups)) groups <- as.character(groups)
  }

  g <- foreach(i=seq_len(num.subjects)) %dopar% {
    res <- graph_from_adjacency_matrix(A[, , i], mode='undirected', diag=FALSE, weighted=TRUE)
    res <- set_brainGraph_attr(res, atlas, use.parallel=FALSE, A=A[, , i], threshold=threshold,
                                  subject=study.ids[i], group=groups[i], ...)
  }
  out$graphs <- g
  names(out$graphs) <- study.ids
  class(out) <- c('brainGraphList', 'list')
  return(out)
}

#' Index individual graphs of a brainGraphList
#'
#' The \code{[} method will let you select individual graphs from a
#' \code{brainGraphList} object. You can subset the graphs with an integer,
#' character (containing Study ID's), or logical vector.
#'
#' @param x A \code{brainGraphList} object
#' @param i Integer, character, or logical vector for subsetting the graphs. If
#'   character, the supplied value(s) for \code{i} should be one of the
#'   \code{study.ids}
#' @export
#' @method [ brainGraphList
#'
#' @name Extract.brainGraphList
#' @rdname brainGraphList

`[.brainGraphList` <- function(x, i) {
  if (is.logical(i) && length(i) != length(x$graphs)) {
    stop('Logical vector is not the correct length')
  }
  if (length(i) > 1) {
    out <- x$graphs[i]
  } else {
    out <- x$graphs[[i]]
  }

  return(out)
}
