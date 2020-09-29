#' Matrix/array utility functions
#'
#' These functions are utility/helper functions when working with matrices or
#' arrays.
#'
#' @param x Numeric matrix or array (the latter, for \code{qr.array} and
#'   \code{symmetrize.array})
#' @param n,p Integer; the number of rows or rank (respectively) of
#'   the input matrix or QR decomposition
#' @param ... Arguments passed to either \code{\link{sort}} (for
#'   \code{get_thresholds}) or \code{\link{qr.default}} (for \code{qr.array}).
#'   For the former, this will typically only be \code{decreasing=TRUE}, if that
#'   is the desired behavior
#' @name Matrix utilities
#' @rdname matrix_utils
NULL

#' @export
#' @rdname matrix_utils

colMax <- function(x, n=dim(x)[1L]) {
  do.call(pmax.int, c(lapply(seq_len(n), function(i) x[i, ]), na.rm=TRUE))
}

#' @export
#' @rdname matrix_utils

colMaxAbs <- function(x, n=dim(x)[1L]) colMax(abs(x), n)

#' @export
#' @rdname matrix_utils

colMin <- function(x, n=dim(x)[1L]) {
  do.call(pmin.int, c(lapply(seq_len(n), function(i) x[i, ]), na.rm=TRUE))
}

#' Return the diagonal of a square matrix
#'
#' \code{diag_sq} is a pared-down version of \code{\link{diag}} for square
#' matrices. It does not return any dimnames, does not check if \code{x} is a
#' square matrix, and it cannot be used to \emph{create} a matrix with a given
#' value along the diagonal. Meant to be used in code that is called repeatedly
#' (thousands of times).
#'
#' @param inds Vector-based indices of the diagonal
#' @export
#' @return \code{diag_sq} returns an unnamed numeric vector with the values
#'   along the diagonal of the input matrix
#' @rdname matrix_utils

diag_sq <- function(x, n=dim(x)[1L], inds=1L + 0L:(n - 1L) * (n + 1L)) {
  x[inds]
}

#' Calculate thresholds
#'
#' \code{get_thresholds} calculates the threshold values that would result in a
#' specific graph density. These depend, necessarily on the values in the matrix
#' themselves.
#'
#' Given a vector of densities, \code{get_thresholds} returns the numeric values
#' that will result in graphs of the given densities after thresholding by those
#' values. In the \emph{Examples} section, the thresholds should result in
#' graphs with densities of \eqn{5, 15, \dots, 55} percent.
#'
#' @param densities Numeric vector of densities
#' @param emax Integer; the maximum number of edges
#' @return \code{get_thresholds} returns a numeric vector of the thresholds
#' @export
#' @rdname matrix_utils
#' @examples
#' x <- matrix(runif(25 * 25), 25, 25)
#' x <- symmetrize(x)
#' diag(x) <- 0
#' densities <- seq(0.05, 0.55, by=0.1)
#' threshes <- get_thresholds(x, densities)
#' ## Verify that the densities are correct
#' graphs <- lapply(threshes, function(th) {
#'   graph_from_adjacency_matrix(x * (x > th), mode='undirected',
#'                               diag=FALSE, weighted=TRUE)
#'   })
#' sapply(graphs, graph.density)

get_thresholds <- function(x, densities, emax=dim(x)[1L] * (dim(x)[1L] - 1L) / 2, ...) {
  stopifnot(all(densities >= 0) && all(densities <= 1))
  sort(x[lower.tri(x)], ...)[emax - densities * emax]
}

#' @return \code{is_binary} returns a logical of length 1
#' @export
#' @rdname matrix_utils

is_binary <- function(x) identical(sum(abs(x)) - sum(x == 1), 0)

#' Calculate Moore-Penrose pseudoinverse of a full rank matrix
#'
#' \code{pinv} calculates \eqn{M^{+} = (M^T M)^{-1} M^T} for full (column) rank
#' matrices. However, it does not verify the matrix's rank.
#'
#' @export
#' @return \code{pinv} returns the input matrix's pseudoinverse
#' @rdname inverse

pinv <- function(x) tcrossprod(inv(x), x)

#' Calculate the QR decomposition for each matrix in an array
#'
#' \code{qr.array} will calculate the QR decomposition for each matrix in a 3D
#' array.
#'
#' @export
#' @return \code{qr.array} returns a \emph{list} in which each element is the QR
#'   decomposition of each matrix along \code{x}'s 3rd dimension
#' @rdname matrix_utils

qr.array <- function(x, ...) apply(x, 3L, qr.default, ...)

#' Slightly faster calculation of the Q and R matrices
#'
#' \code{qr_Q2} and \code{qr_R2} are simplified versions of \code{\link{qr.Q}}
#' and \code{\link{qr.R}}.
#'
#' @param QR A \code{qr} object
#' @param y A numeric matrix with \code{1} along the diagonal, of the same size
#'   as the input matrix (i.e., \code{QR$qr})
#' @export
#' @rdname matrix_utils

qr_Q2 <- function(QR, y=diag(1, n, p), n=dim(QR$qr)[1L], p=QR$rank) {
  qr.qy(QR, y)
}

#' @export
#' @rdname matrix_utils

qr_R2 <- function(QR, p=QR$rank) {
  R <- QR$qr[seq.int(p), , drop=FALSE]
  R[row(R) > col(R)] <- 0
  R
}

#' Symmetrize a matrix with the mean of off-diagonal elements
#'
#' \code{symm_mean} returns a symmetric matrix in which the off-diagonal
#' elements \eqn{A[i, j]} and \eqn{A[j, i]} are set to the mean of the values
#' in the input matrix.
#'
#' @export
#' @rdname matrix_utils

symm_mean <- function(x) {
  0.5 * (x + t(x))
}

#' Symmetrize matrices and arrays
#'
#' \code{symmetrize} will symmetrize a numeric matrix (or each matrix in an
#' array) by assigning to the off-diagonal elements either the \code{max}
#' (default), \code{min}, or \code{average} of \eqn{\{A(i, j), A(j, i)\}}.
#'
#' @param symm.by Character string; how to create symmetric off-diagonal
#'   elements. Default: \code{max}
#' @export
#' @rdname matrix_utils

symmetrize <- function(x, ...) {
  UseMethod('symmetrize')
}

#' @export
#' @rdname matrix_utils
symmetrize.matrix <- function(x, symm.by=c('max', 'min', 'avg'), ...) {
  dims <- dim(x)
  stopifnot(dims[1L] == dims[2L])

  symm.by <- match.arg(symm.by)
  switch(symm.by,
         avg=symm_mean(x),
         max=pmax(x, t(x)),
         min=pmin(x, t(x)))
}

#' @export
#' @rdname matrix_utils
symmetrize.array <- function(x, symm.by=c('max', 'min', 'avg'), ...) {
  dims <- dim(x)
  stopifnot(dims[1L] == dims[2L])
  symm.by <- match.arg(symm.by)
  array(apply(x, 3L, symmetrize, symm.by), dim=dims, dimnames=dimnames(x))
}
