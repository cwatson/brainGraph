#' brainGraph generic methods
#'
#' These functions are S3 \emph{generics} for various \code{brainGraph}-defined
#' objects.
#'
#' @name brainGraph-methods
NULL

#' Return group names
#'
#' \code{groups} returns the \dQuote{Group} graph attribute for each graph or
#' observation in the object.
#'
#' @param x,object An object
#' @export
#' @method groups brainGraphList
#' @rdname brainGraph-methods

groups.brainGraphList <- function(x) {
  out <- NULL
  if (x$type == 'observed') {
    if ('Group' %in% graph_attr_names(x[1L])) {
      out <- vapply(x[], graph_attr, character(1L), 'Group')
    }
  } else if (x$type == 'random') {
    subs <- x$graphs
    if ('Group' %in% graph_attr_names(subs[[1L]][[1L]])) {
      out <- vapply(subs, function(g) graph_attr(g[[1L]], 'Group'), character(1L))
    }
  }
  return(out)
}

#' @export
#' @method groups brainGraph_resids
#' @rdname residuals
groups.brainGraph_resids <- function(x) {
  sID <- getOption('bg.subject_id')
  gID <- getOption('bg.group')
  structure(x$resids.all[, as.character(get(gID))], names=x$resids.all[, get(sID)])
}

#' @export
#' @method groups corr_mats
#' @rdname brainGraph-methods
groups.corr_mats <- function(x) names(x$r.thresh)

#' Extract region names from brainGraph objects
#'
#' \code{region.names} is a generic method for extracting region names from
#' various \code{brainGraph} objects. These are generally convenience functions.
#'
#' For a \code{data.table}, \code{region.names} assumes that it contains a
#' \emph{factor} column named \code{region}.
#'
#' @export
#' @rdname brainGraph-methods

region.names <- function(object) {
  UseMethod('region.names')
}

#' @export
#' @method region.names data.table
#' @rdname brainGraph-methods
region.names.data.table <- function(object) {
  colname <- if (hasName(object, 'Region')) 'Region' else 'region'
  object[, levels(get(colname))]
}

#' @export
#' @method region.names brainGraph_resids
#' @rdname residuals
region.names.brainGraph_resids <- function(object) names(object$resids.all)[-c(1L, 2L)]

#' @export
#' @method region.names corr_mats
#' @rdname correlation_matrices
region.names.corr_mats <- function(object) dimnames(object$R)[[1L]]

#' @export
#' @method region.names mtpc
#' @include mtpc.R
#' @rdname mtpc
region.names.mtpc <- region.names.bg_GLM

#' Extract the number of regions in a brainGraph object
#'
#' \code{nregions} is a generic method for extracting the number of regions from
#' various \code{brainGraph} objects.
#'
#' @export
#' @rdname brainGraph-methods

nregions <- function(object) {
  UseMethod('nregions')
}

#' @export
#' @rdname residuals
nregions.brainGraph_resids <- function(object) dim(object$resids.all)[2L] - 2L

#' @export
#' @rdname correlation_matrices
nregions.corr_mats <- function(object) dim(object$R)[1L]

#' @export
#' @include mtpc.R
#' @rdname mtpc
nregions.mtpc <- nregions.bg_GLM

#' @export
#' @rdname NBS
nregions.NBS <- function(object) dim(object$T.mat)[1L]

#' Calculate coefficient of variation
#'
#' \code{coeff_var} is a S3 generic that calculates the \emph{coefficient of
#' variation}, defined as
#' \deqn{CV(x) = \frac{sd(x)}{mean(x)}}
#'
#' If \code{x} is a matrix, it will calculate the CV for each \emph{column}. If
#' \code{x} is a 3D array, it will calculate the coefficient of variation for
#' each \emph{row-column} combination. If the input dimensions are \eqn{n \times
#' n \times r}, a matrix with size \eqn{n \times n} will be returned.
#'
#' @param x Numeric vector, matrix, or array
#' @param na.rm Logical indicating whether \code{NA} values should be stripped
#'   when calculating sums. Default: \code{FALSE}
#' @param ... Unused
#' @return A numeric vector or matrix
#' @export

coeff_var <- function(x, na.rm=FALSE, ...) {
  UseMethod('coeff_var')
}

#' @export
#' @rdname coeff_var
coeff_var.default <- function(x, na.rm=FALSE, ...) {
  N <- sum(!is.na(x))
  mu <- sum(x, na.rm=na.rm) / N
  sqrt(1 / (N - 1L) * (sum((x - mu)^2, na.rm=na.rm))) / mu
}

#' @export
coeff_var.matrix <- function(x, na.rm=FALSE, ...) {
  apply(x, 2L, coeff_var)
}

#' @export
coeff_var.array <- function(x, na.rm=FALSE, ...) {
  dimX <- dim(x)
  N <- if (anyNA(x)) apply(x, 1L:2L, function(z) sum(!is.na(z))) else dimX[3L]
  mu <- rowSums(x, na.rm=na.rm, dims=2L) / N
  sqrt(1 / (N - 1L) * (rowSums((x - array(mu, dim=dimX))^2, na.rm=na.rm, dims=2L))) / mu
}

#' Calculate the inverse of the cross product of a design matrix
#'
#' \code{inv} is a \code{S3} generic that calculates the inverse of the cross
#' product of a design matrix, also referred to as the \dQuote{unscaled
#' covariance matrix}.
#'
#' If \code{x} is a matrix, the Cholesky decomposition of the cross product is
#' calculated (or using \code{\link{tcrossprod}} if \code{transpose=TRUE}), and
#' the inverse is calculated from that result. That is,
#' \deqn{inv(X) = (X^T X)^{-1}}
#' \deqn{inv(X, transpose=TRUE) = (X X^T)^{-1}}
#' \deqn{inv(X, y) = (X^T y)^{-1}}
#'
#' If \code{x} is a 3-dimensional array, then the inverse will be calculated for
#' each matrix along the 3rd dimension, with the same input arguments for each.
#'
#' Finally, there is a method for objects with class \code{qr}, and lists of QR
#' decomposition objects.
#'
#' @note These methods should only be used on \emph{full-rank} matrices, as
#' there is no error checking being performed.
#'
#' @param x A numeric matrix or array, a \code{qr} object, or a list of
#'   \code{qr} objects
#' @param ... Unused
#' @export
#' @return A numeric matrix or array
#' @name Inverse
#' @rdname inverse

inv <- function(x, ...) {
  UseMethod('inv')
}

#' @param y A numeric matrix or vector (for the \code{matrix} and \code{array}
#'   methods). If supplied, this will be multiplied by \code{x} before the
#'   inverse is calculated. Default: \code{NULL}
#' @param transpose Logical. If \code{FALSE} (the default), take the cross
#'   product of the arguments. If \code{TRUE}, use \code{\link{tcrossprod}}
#' @export
#' @rdname inverse

inv.matrix <- function(x, y=NULL, transpose=FALSE, ...) {
  xty <- if (isTRUE(transpose)) tcrossprod(x, y) else crossprod(x, y)
  xinv <- chol2inv(chol.default(xty))
  dimnames(xinv) <- if (isTRUE(transpose)) dimnames(x)[c(1L, 1L)] else dimnames(x)[c(2L, 2L)]
  xinv
}

#' @export
#' @rdname inverse

inv.array <- function(x, y=NULL, transpose=FALSE, ...) {
  namesX <- dimnames(x)
  dimX <- dim(x)
  dim_ind <- if (isTRUE(transpose)) 1L else 2L
  if (is.null(y)) {
    dims <- dimX[c(dim_ind, dim_ind)]
    nms <- namesX[c(dim_ind, dim_ind)]
  } else {
    dims <- c(dimX[dim_ind], dim(y)[dim_ind])
    nms <- c(namesX[dim_ind], dimnames(y)[dim_ind])
  }
  xinv <- array(0, dim=c(dims, dimX[3L]), dimnames=c(nms, namesX[3L]))
  for (k in seq_len(dimX[3L])) {
    xinv[, , k] <- inv(abind::adrop(x[, , k, drop=FALSE], drop=3L), y, transpose)
  }
  xinv
}

#' @export
#' @rdname inverse

inv.qr <- function(x, p=x$rank, ...) {
  xinv <- chol2inv(x$qr, p)
  dimnames(xinv) <- dimnames(x$qr)[c(2L, 2L)]
  xinv
}

#' @param p The rank of the original matrix
#' @param r The number of design matrices; i.e., the length of the input list
#' @param vnames Character vector of the design matrix's variable names
#' @param nms The region names; i.e., the names of the input list
#' @export
#' @rdname inverse

inv.list <- function(x, p=x[[1L]]$rank, r=length(x), vnames=dimnames(x[[1L]]$qr)[[2L]],
                     nms=names(x), ...) {
  array(vapply(x, inv, numeric(p * p), p), dim=c(p, p, r), dimnames=list(vnames, vnames, nms))
}
