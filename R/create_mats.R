#' Create connection matrices for tractography or fMRI data
#'
#' \code{create_mats} creates arrays from connection matrices (e.g.,
#' \emph{fdt_network_matrix} from FSL or \emph{ROICorrelation.txt} from DPABI).
#' You may choose to normalize these matrices by the \emph{waytotal} or
#' \emph{region size} (tractography), or not at all.
#'
#' @section Connection matrix files:
#' The \code{A.files} argument is mandatory and may be specified in a few ways:
#' \enumerate{
#'   \item A character vector of the filenames (preferably with full path).
#'   \item A single character string specifying the \emph{directory} in which
#'     all connectivity matrices are located. This will load \emph{all} files in
#'     the directory.
#'   \item A \emph{named list} in which the names match the arguments to
#'     \code{\link{list.files}}. This will load \emph{all} files in \code{path}
#'     that match the \code{pattern} argument, if present, and will load
#'     \emph{all} files in child directories if \code{recursive=TRUE}. See
#'     examples below.
#' }
#' The same options apply to \code{div.files} as well.
#'
#' @section Thresholding methods:
#' The argument \code{threshold.by} has 4 options:
#' \enumerate{
#'   \item \code{consensus} Threshold based on the raw (normalized, if selected)
#'     values in the matrices. If this is selected, it uses the
#'     \code{sub.thresh} value to perform \dQuote{consensus} thresholding.
#'   \item \code{density} Threshold the matrices to yield a specific graph
#'     density (given by the \code{mat.thresh} argument).
#'   \item \code{mean} Keep only connections for which the cross-subject mean is
#'     at least 2 standard deviations higher than the threshold (specified by
#'     \code{mat.thresh})
#'   \item \code{consistency} Threshold based on the coefficient of variation to
#'     yield a graph with a specific density (given by \code{mat.thresh}). The
#'     edge weights will still represent those of the input matrices. See
#'     Roberts et al. (2017) for more on \dQuote{consistency-based}
#'     thresholding.
#'   \item \code{raw} Threshold each subject's matrix \emph{individually},
#'     irrespective of group membership. Ignores \code{sub.thresh}.
#' }
#'
#' The argument \code{mat.thresh} allows you to choose a numeric threshold,
#' below which the connections will be replaced with 0; this argument will also
#' accept a numeric vector. The argument \code{sub.thresh} will keep only those
#' connections for which at least \emph{X}\% of subjects have a positive entry
#' (the default is 0.5, or 50\%).
#'
#' @param A.files Character vector of the filenames with connection matrices
#' @param modality Character string indicating data modality (default:
#'   \code{dti})
#' @param divisor Character string indicating how to normalize the connection
#'   matrices; either 'none' (default), 'waytotal', 'size', or 'rowSums'
#'   (ignored if \code{modality} equals \code{fmri})
#' @param div.files Character vector of the filenames with the data to
#'   normalize by (e.g. a list of \emph{waytotal} files) (default: \code{NULL})
#' @param threshold.by Character string indicating how to threshold the data;
#'   choose \code{density}, \code{mean}, or \code{consistency} if you want all
#'   resulting matrices to have the same densities (default: \code{consensus})
#' @param mat.thresh Numeric (vector) for thresholding connection matrices
#'   (default: 0)
#' @param sub.thresh Numeric (between 0 and 1) for thresholding by subject
#'   numbers (default: 0.5)
#' @param inds List (length equal to number of groups) of integers; each list
#'   element should be a vector of length equal to the group sizes
#' @param algo Character string of the tractography algorithm used (default:
#'   \code{'probabilistic'}). Ignored if \emph{modality} is \code{fmri}.
#' @param P Integer; number of samples per seed voxel (default: 5000)
#' @param ... Arguments passed to \code{\link{symmetrize_mats}}
#' @export
#' @importFrom abind abind
#'
#' @return A list containing:
#' \item{A}{A 3-d array of the raw connection matrices}
#' \item{A.norm}{A 3-d array of the normalized connection matrices}
#' \item{A.bin}{A list of 3-d arrays of binarized connection matrices, one array
#'   for each threshold}
#' \item{A.bin.sums}{A list of 3-d arrays of connection matrices, with each
#'   entry signifying the number of subjects with a connection present; the
#'   number of list elements equals the length of \code{mat.thresh}, and the
#'   extent of the arrays equals the number of groups}
#' \item{A.inds}{A list of arrays of binarized connection matrices, containing 1
#'   if that entry is to be included}
#' \item{A.norm.sub}{List of 3-d arrays of the normalized connection matrices
#'   for all given thresholds}
#' \item{A.norm.mean}{List of 3-d arrays of connection matrices averaged for
#'   each group}
#'
#' @family Matrix functions
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Roberts, JA and Perry, A and Roberts, G and Mitchell, PB and
#'   Breakspear, M (2017) Consistency-based thresholding of the human
#'   connectome. \emph{NeuroImage}. \bold{145}, 118--129.
#'   \url{https://dx.doi.org/10.1016/j.neuroimage.2016.09.053}
#' @examples
#' \dontrun{
#' thresholds <- seq(from=0.001, to=0.01, by=0.001)
#' fmri.mats <- create_mats(f.A, modality='fmri', threshold.by='consensus',
#'   mat.thresh=thresholds, sub.thresh=0.5, inds=inds)
#' dti.mats <- create_mats(f.A, divisor='waytotal', div.files=f.way,
#'   mat.thresh=thresholds, sub.thresh=0.5, inds=inds)
#'
#' # Specify a directory and filename pattern
#' conn_files <- list(path='~/data', pattern='.*fdt_network_matrix')
#' dti.mats <- create_mats(conn_files, ...)
#' }

create_mats <- function(A.files, modality=c('dti', 'fmri'),
                        divisor=c('none', 'waytotal', 'size', 'rowSums'),
                        div.files=NULL,
                        threshold.by=c('consensus', 'density', 'mean', 'consistency', 'raw'),
                        mat.thresh=0, sub.thresh=0.5, inds=list(seq_along(A.files)),
                        algo=c('probabilistic', 'deterministic'), P=5e3, ...) {

  # Argument checking
  #-----------------------------------------------------------------------------
  A.bin <- A.bin.sums <- A.inds <- NULL
  kNumSubjs <- lengths(inds)
  nz <- sum(kNumSubjs)
  if (is.list(A.files) || length(A.files) == 1) A.files <- dir2files(A.files)
  stopifnot(all(file.exists(A.files)),
            length(A.files) == nz,
            sub.thresh >= 0 && sub.thresh <= 1)
  if (!is.null(div.files)) {
    if (is.list(div.files) || length(div.files == 1)) div.files <- dir2files(div.files)
    stopifnot(length(div.files) == length(A.files))
  }

  A <- read.array(A.files)
  Nv <- dim(A)[1L]
  emax <- Nv * (Nv - 1) / 2
  A[is.nan(A)] <- 0
  A.norm <- A

  modality <- match.arg(modality)
  algo <- match.arg(algo)
  divisor <- match.arg(divisor)

  # Normalize DTI matrices
  #-------------------------------------
  if (modality == 'dti' && algo == 'probabilistic' && divisor != 'none') {
    A.norm <- normalize_mats(A.norm, divisor, div.files, P)
    A.norm[is.nan(A.norm)] <- 0
  }

  # Matrix thresholding
  #-----------------------------------------------------------------------------
  threshold.by <- match.arg(threshold.by)
  if (threshold.by == 'raw') {
    A.norm.sub <- lapply(mat.thresh, function(x) A.norm * (A.norm > x))

  } else if (threshold.by %in% c('density', 'consistency')) {
    stopifnot(all(mat.thresh >= 0) && all(mat.thresh <= 1))
    A.norm.sub <- vector('list', length(mat.thresh))

    if (threshold.by == 'density') {
      Asym <- symmetrize_array(A.norm, ...)
      thrMat <- vapply(mat.thresh, function(x) apply(Asym, 3, get_thresholds, x, emax), numeric(nz))
      Aperm <- aperm(Asym, c(3, 1, 2))
      for (x in seq_along(mat.thresh)) {
        A.norm.sub[[x]] <- Aperm * (Aperm > thrMat[, x])
        A.norm.sub[[x]] <- aperm(A.norm.sub[[x]], c(2, 3, 1))
      }
    } else if (threshold.by == 'consistency') {
      all.cv <- apply(A.norm, 1:2, coeff_var)
      all.cv <- symmetrize_mats(all.cv, 'min')
      threshes <- get_thresholds(all.cv, mat.thresh, emax, decreasing=TRUE)
      A.inds <- lapply(threshes, function(x) ifelse(all.cv < x, 1L, 0L))
      for (z in seq_along(A.inds)) {
        tmp <- replicate(nz, A.inds[[z]])
        A.norm.sub[[z]] <- A.norm * (tmp == 1L)
      }
    }
  } else {
    if (threshold.by == 'consensus') {
      # Use the given thresholds as-is
      #---------------------------------
      # Binarize the array, then keep entries w/ >= "sub.thresh"% for each group
      A.bin <- lapply(mat.thresh, function(x) {A.norm > x} + 0L)
      A.bin.sums <- lapply(seq_along(mat.thresh), function(x)
                           abind(lapply(inds, function(z)
                                        rowSums(A.bin[[x]][, , z], dims=2)), along=3))

      # For deterministic, threshold by size *after* binarizing
      if (modality == 'dti' && algo == 'deterministic' && divisor == 'size') {
        A.norm <- normalize_mats(A.norm, divisor, div.files, P=1)
      }

      # This is a list (# mat.thresh) of Nv x Nv x Ng arrays
      A.inds <- vector('list', length(mat.thresh))
      if (sub.thresh == 0) {
        A.inds <- lapply(seq_along(mat.thresh), function(x)
                         ifelse(A.bin.sums[[x]] > 0, 1L, 0L))
      } else {
        thresh <- sub.thresh * kNumSubjs
        for (x in seq_along(mat.thresh)) {
          tmp <- aperm(A.bin.sums[[x]], c(3, 1, 2))
          A.inds[[x]] <- array(0L, dim=dim(tmp))
          A.inds[[x]][tmp >= thresh] <- 1L
          A.inds[[x]] <- aperm(A.inds[[x]], c(2, 3, 1))
        }
      }

      # Back to a list of arrays for all subjects
      A.norm.sub <- vector('list', length(mat.thresh))
      for (x in seq_along(mat.thresh)) {
        A.norm.sub[[x]] <- array(0, dim=dim(A.norm))
        keep <- array(FALSE, dim=dim(A.norm))
        for (y in seq_along(inds)) {
          tmp <- A.inds[[x]][, , y] == 1L
          keep[, , inds[[y]]] <- replicate(kNumSubjs[y], tmp)
        }
        A.norm.sub[[x]][keep] <- A.norm[keep]
      }

    } else if (threshold.by == 'mean') {
      # Threshold: mean + 2SD > mat.thresh
      #---------------------------------
      all.mean <- rowMeans(A.norm, dims=2)
      all.means <- replicate(nz, all.mean)
      all.sums <- rowSums((A.norm - all.means)^2, dims=2)
      all.sd <- sqrt(all.sums / (nz - 1))
      all.thresh <- all.mean + (2 * all.sd)
      all.thresh <- replicate(nz, all.thresh)

      A.norm.sub <- lapply(mat.thresh, function(x) A.norm * (all.thresh > x))
    }
  }
  if (threshold.by != 'density') {
    for (x in seq_along(mat.thresh)) A.norm.sub[[x]] <- symmetrize_array(A.norm.sub[[x]], ...)
  }

  A.norm.mean <- lapply(seq_along(mat.thresh), function(x)
                        abind(lapply(inds, function(y)
                                     rowMeans(A.norm.sub[[x]][, , y], dims=2)), along=3))
  if (threshold.by == 'density') {
    A.norm.mean <- lapply(seq_along(mat.thresh), function(x)
                          array(apply(A.norm.mean[[x]], 3, function(y)
                                      y * (y > get_thresholds(y, mat.thresh[x], emax))),
                                dim=dim(A.norm.mean[[x]])))
  }

  return(list(A=A, A.norm=A.norm, A.bin=A.bin, A.bin.sums=A.bin.sums,
              A.inds=A.inds, A.norm.sub=A.norm.sub, A.norm.mean=A.norm.mean))
}

#' Create a symmetric matrix
#'
#' \code{symmetrize_mats} will symmetrize a numeric matrix by assigning the
#' off-diagonal elements values of either the \code{max} (default), \code{min},
#' or \code{average} of \eqn{\{A(i, j), A(j, i)\}}.
#'
#' @param A Numeric matrix
#' @param symm.by Character string; how to create symmetric off-diagonal
#'   elements. Default: \code{max}
#' @export
#' @return Either a single symmetrized matrix, or an (3D) array
#'
#' @family Matrix functions
#' @seealso \code{\link[igraph]{graph_from_adjacency_matrix}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

symmetrize_mats <- function(A, symm.by=c('max', 'min', 'avg')) {
  dims <- dim(A)
  stopifnot(dims[1L] == dims[2L])

  symm.by <- match.arg(symm.by)
  A <- switch(symm.by,
              avg=symm_mean(A),
              max=pmax(A, t(A)),
              min=pmin(A, t(A)))
  return(A)
}

#' Symmetrize each matrix in a 3D array
#'
#' \code{symmetrize_array} is a convenience function which applies
#' \code{\link{symmetrize_mats}} along the 3rd dimension of an array.
#'
#' @param ... Arguments passed to \code{\link{symmetrize_mats}}
#' @inheritParams symmetrize_mats
#' @export
#' @rdname symmetrize_mats

symmetrize_array <- function(A, ...) {
  return(array(apply(A, 3, symmetrize_mats, ...), dim=dim(A)))
}

#' Symmetrize a matrix with the mean of off-diagonal elements
#'
#' \code{symm_mean} returns a symmetric matrix in which the off-diagonal
#' elements \eqn{A[i, j]} and \eqn{A[j, i]} are set to the mean of the values
#' in the input matrix.
#' @export
#' @rdname symmetrize_mats

symm_mean <- function(A) {
  0.5 * (A + t(A))
}

read.array <- function(infiles, ncols=NULL) {
  Nv <- length(readLines(infiles[1]))
  if (is.null(ncols)) ncols <- Nv
  A <- array(sapply(infiles, function(x)
                    matrix(scan(x, what=numeric(), n=Nv*ncols, quiet=TRUE),
                           Nv, ncols, byrow=TRUE)),
             dim=c(Nv, ncols, length(infiles)))
  return(A)
}

normalize_mats <- function(A, divisor, div.files, P) {
  div <- read.array(div.files, ncols=1)

  if (divisor == 'waytotal') {
    # Control for streamline count by waytotal
    Nv <- dim(A)[1L]
    W <- array(apply(div, 3, function(x) x[, rep(1, Nv)]), dim=dim(A))
    A.norm <- A / W

  } else if (divisor == 'size') {
    # Control for the size (# voxels) of both regions 'x' and 'y'
    R <- array(apply(div, 3, function(x) outer(x, x, FUN='+')), dim=dim(A))
    A.norm <- 2 * A / (P * R)

  } else if (divisor == 'rowSums') {
    A.norm <- array(apply(A, 3, function(x) x / rowSums(x)), dim=dim(A))
  }
  return(A.norm)
}

#' Threshold additional set of matrices
#'
#' \code{apply_thresholds} thresholds an additional set of matrices (e.g.,
#' FA-weighted matrices for DTI tractography) based on the matrices that have
#' been returned from \code{\link{create_mats}}. This ensures that the same
#' connections are present in both sets of matrices.
#'
#' The argument \code{W.files} accepts the same formats as \code{A.files}; see
#' \code{\link{create_mats}} for details.
#'
#' @param sub.mats List (length equal to number of thresholds) of numeric arrays
#'   (3-dim) for all subjects
#' @param group.mats List (equal to number of thresholds) of lists (equal to
#'   number of groups) of numeric matrices for group-level data
#' @param W.files Character vector of the filenames of the files with
#'   connectivity matrices
#' @inheritParams create_mats
#' @export
#'
#' @return List containing:
#' \item{W}{A 3-d array of the raw connection matrices}
#' \item{W.norm.sub}{List of 3-d arrays of the normalized connection matrices
#'   for all given thresholds}
#' \item{W.norm.mean}{List of lists of numeric matrices averaged for each group}
#'
#' @family Matrix functions
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#'   W.mats <- apply_thresholds(A.norm.sub, A.norm.mean, f.W, inds)
#' }

apply_thresholds <- function(sub.mats, group.mats, W.files, inds) {
  if (is.list(W.files) || length(W.files) == 1) W.files <- dir2files(W.files)
  dimA <- dim(sub.mats[[1]])
  nz <- dimA[3L]
  stopifnot(length(W.files) == nz)
  W <- read.array(W.files)
  W.norm.sub <- W.norm.mean <- vector('list', length(sub.mats))
  for (x in seq_along(sub.mats)) {
    W.norm.sub[[x]] <- array(0, dim=dimA)
    keep <- sub.mats[[x]] > 0
    W.norm.sub[[x]][keep] <- W[keep]
  }
  ng <- dim(group.mats[[1]])[3L]
  W.norm.mean <- lapply(seq_along(group.mats), function(x)
                        lapply(seq_len(ng), function(g)
                               ifelse(group.mats[[x]][, , g] > 0,
                                      rowMeans(W.norm.sub[[x]][, , inds[[g]]], dims=2),
                                      0)))
  W.norm.mean <- lapply(W.norm.mean, function(x) abind(x, along=3))
  return(list(W=W, W.norm.sub=W.norm.sub, W.norm.mean=W.norm.mean))
}

#' Return a vector of filenames based on a directory name or options list
#'
#' @param x Either a character string specifying a directory or a list of
#'   arguments to be passed to \code{\link{list.files}}
#' @keywords internal

dir2files <- function(x) {
  fargs <- formals(list.files)
  if (is.list(x)) {
    if (is.null(attr(x, 'names'))) names(x) <- names(fargs)[seq_along(x)]
    if (!hasName(x, 'full.names')) x <- c(x, full.names=TRUE)
  } else {
    x <- list(path=x, full.names=TRUE)
  }
  x <- c(x, fargs[setdiff(names(fargs), names(x))])
  stopifnot(file_test('-d', x$path))
  files <- do.call(list.files, x)
  return(files)
}
