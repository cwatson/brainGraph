#' Return group names
#'
#' \code{groups} returns the \dQuote{Group} graph attribute for each graph in
#' the object.
#'
#' @param x,object An object
#' @export
#' @rdname methods

groups.brainGraphList <- function(x) {
  out <- NULL
  if (x$type == 'observed') {
    if ('Group' %in% graph_attr_names(x[1])) {
      out <- vapply(x[], graph_attr, character(1), 'Group')
    }
  } else if (x$type == 'random') {
    subs <- x$graphs
    if ('Group' %in% graph_attr_names(subs[[1]])) {
      out <- vapply(subs, function(g) graph_attr(g[[1]], 'Group'), character(1))
    }
  }
  return(out)
}

#' @rdname methods
groups.brainGraph_resids <- function(x) {
  Group <- Study.ID <- NULL
  structure(x$resids.all[, as.character(Group)], names=x$resids.all[, Study.ID])
}

#' @rdname methods
groups.corr_mats <- function(x) {
  names(x$r.thresh)
}

#' @rdname methods
nobs.brainGraphList <- function(object, ...) {
  length(object$graphs)
}

#' @rdname methods
nobs.brainGraph_resids <- function(object, ...) {
  dim(object$resids.all)[1L]
}

#' @method case.names brainGraph_resids
#' @rdname methods
case.names.brainGraph_resids <- function(object, ...) {
  Study.ID <- NULL
  object$resids.all[, as.character(Study.ID)]
}

#' Extract region names from brainGraph objects
#'
#' \code{region.names} is a generic method for extracting region names from
#' various \code{brainGraph objects. These are generally convenience functions.
#'
#' For a \code{data.table}, this function assumes that it contains a
#' \emph{factor} column named \code{region}.
#'
#' @export
#' @rdname methods

region.names <- function(object) {
  UseMethod('region.names')
}

#' @method region.names data.table
#' @rdname methods
region.names.data.table <- function(object) {
  object[, levels(region)]
}

#' @export
#' @method region.names brainGraph_resids
#' @rdname residuals
region.names.brainGraph_resids <- function(object) {
  object$all.dat.long[, levels(Region)]
}

#' @export
#' @method region.names corr_mats
#' @rdname correlation_matrices
region.names.corr_mats <- function(object) {
  dimnames(object$R)[[1L]]
}

#' Extract the number of regions in a brainGraph object
#'
#' \code{nregions} is a generic method for extracting the number of regions from
#' various \code{brainGraph} objects.
#'
#' @export
#' @rdname methods

nregions <- function(object) {
  UseMethod('nregions')
}

#' @rdname residuals
nregions.brainGraph_resids <- function(object) {
  object$all.dat.long[, nlevels(Region)]
}

#' @rdname correlation_matrices
nregions.corr_mats <- function(object) {
  dim(object$R)[1L]
}
