#' Create a list of brainGraph graphs
#'
#' \code{make_bgList} creates a \code{brainGraphList} object, a list containing
#' a set of graphs for all subjects in a study at a specific threshold (or
#' density), in addition to some graph-level attributes common to those graphs.
#'
#' In addition to creating the initial \code{igraph} graphs from the
#' connectivity matrices, \code{\link{set_brainGraph_attr}} will be called on
#' each graph.
#'
#' This object can be considered comparable to a 4-D \emph{NIfTI} file,
#' particularly those returned by FSL's \emph{TBSS} prestats step since this
#' contains the FA volumes for all study subjects.
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
#' @return An object of class \code{brainGraphList} with elements:
#'   \item{threshold}{The specified threshold/density}
#'   \item{version}{The version of \code{brainGraph} used when creating the
#'     graphs}
#'   \item{atlas}{The atlas common to all the graphs}
#'   \item{modality}{The imaging modality (if supplied)}
#'   \item{weighting}{A string indicating what edge weights represent (if
#'     applicable)}
#'   \item{graphs}{The \emph{list} of \code{brainGraph} graphs. The names of
#'     this list will correspond to the Study ID's}
#'
#' @name brainGraphList
#' @rdname brainGraphList
#' @family Graph creation functions
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' # Create a list, one for each threshold
#' g <- vector('list', length(thresholds))
#' for (i in seq_along(thresholds)) {
#'   g[[i]] <- make_brainGraphList(A.norm.sub[[i]], thresholds[i], atlas,
#'       covars.dti$Study.ID, covars.dti$Group, modality='dti', weighting='fa')
#' }
#' }

make_brainGraphList <- function(A, threshold, atlas, study.ids=NULL, groups=NULL, ...) {
  num.subjects <- dim(A)[3]
  if (!is.null(study.ids)) stopifnot(length(study.ids) == num.subjects)
  if (!is.null(groups)) stopifnot(length(groups) == num.subjects)

  out <- list(threshold=threshold, version=packageVersion('brainGraph'), atlas=atlas,
              rand=FALSE, modality=NULL, weighting=NULL)

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
#' @examples
#' \dontrun{
#' # Get the first 5 subjects
#' my.bgL[1:5]
#'
#' # Get subjects by Study ID
#' my.bgL[c('s003', 's008', 's054')]
#'}

`[.brainGraphList` <- function(x, i) {
  if (is.logical(i) && length(i) != length(x$graphs)) {
    stop('Logical indexing vector must be of the same length as the number of graphs')
  }
  if (length(i) > 1) {
    out <- x$graphs[i]
  } else {
    out <- x$graphs[[i]]
  }

  return(out)
}

#' Coerce list of graphs to a brainGraphList object
#'
#' \code{as_brainGraphList} will coerce a list of graphs to a
#' \code{brainGraphList} object. It is assumed that certain "bookkeeping"
#' attributes -- threshold, package version, brain atlas, imaging modality, edge
#' weighting, and whether these are random graphs -- are identical for all
#' graphs in the list. For any that are not, the vector of values will be
#' stored; this may be an issue for downstream analysis.
#'
#' @param g.list List of graph objects
#' @export
#' @return A \code{brainGraphList} object
#'
#' @rdname brainGraphList

as_brainGraphList <- function(g.list) {
  if (!inherits(g.list, 'list')) g.list <- list(g.list)
  stopifnot(all(sapply(g.list, inherits, 'brainGraph')))

  #TODO: check for uniqueness, and save the vector if necessary
  out <- list(threshold=NULL, version=g.list[[1]]$version, atlas=g.list[[1]]$atlas,
              rand=FALSE, modality=g.list[[1]]$modality, weighting=g.list[[1]]$weighting,
              graphs=g.list)
  if ('threshold' %in% graph_attr_names(g.list[[1]])) out$threshold <- g.list[[1]]$threshold

  class(out) <- c('brainGraphList', 'list')
  return(out)
}
