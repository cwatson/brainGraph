#' Return group names
#'
#' \code{groups} returns the \dQuote{Group} graph attribute for each graph or
#' observation in the object.
#'
#' @param x,object An object
#' @param ... Unused
#' @export
#' @name Methods
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

#' @export
#' @rdname methods
groups.brainGraph_resids <- function(x) {
  sID <- getOption('bg.subject_id')
  gID <- getOption('bg.group')
  structure(x$resids.all[, as.character(get(gID))], names=x$resids.all[, get(sID)])
}

#' @export
#' @rdname methods
groups.corr_mats <- function(x) {
  names(x$r.thresh)
}

#' @export
#' @rdname methods
nobs.brainGraphList <- function(object, ...) {
  length(object$graphs)
}

#' @export
#' @rdname methods
nobs.brainGraph_resids <- function(object, ...) {
  dim(object$resids.all)[1L]
}

#' @export
#' @method case.names brainGraph_resids
#' @rdname methods
case.names.brainGraph_resids <- function(object, ...) {
  sID <- getOption('bg.subject_id')
  object$resids.all[, as.character(get(sID))]
}

#' Extract region names from brainGraph objects
#'
#' \code{region.names} is a generic method for extracting region names from
#' various \code{brainGraph} objects. These are generally convenience functions.
#'
#' For a \code{data.table}, \code{region.names} assumes that it contains a
#' \emph{factor} column named \code{region}.
#'
#' @export
#' @rdname methods

region.names <- function(object) {
  UseMethod('region.names')
}

#' @export
#' @method region.names data.table
#' @rdname methods
region.names.data.table <- function(object) {
  colname <- if (hasName(object, 'Region')) 'Region' else 'region'
  object[, levels(get(colname))]
}

#' @export
#' @method region.names brainGraph_resids
#' @rdname residuals
region.names.brainGraph_resids <- function(object) {
  Region <- NULL
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

#' @export
#' @rdname residuals
nregions.brainGraph_resids <- function(object) {
  Region <- NULL
  object$all.dat.long[, nlevels(Region)]
}

#' @export
#' @rdname correlation_matrices
nregions.corr_mats <- function(object) {
  dim(object$R)[1L]
}
