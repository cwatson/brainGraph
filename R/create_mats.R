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
#' The argument \code{threshold.by} has 5 options:
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
#' @param ... Arguments passed to \code{\link{symmetrize}}
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
  kNumSubjs <- lengths(inds)
  nz <- sum(kNumSubjs)
  A.files <- dir2files(A.files)
  stopifnot(all(file.exists(A.files)),
            length(A.files) == nz,
            sub.thresh >= 0 && sub.thresh <= 1)
  if (!is.null(div.files)) {
    div.files <- dir2files(div.files)
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
  if (modality == 'dti' && algo == 'probabilistic') {
    A.norm <- normalize_mats(A.norm, divisor, div.files, P)
    A.norm[is.nan(A.norm)] <- 0
  }

  # Calculate the thresholded matrices
  #-------------------------------------
  threshold.by <- match.arg(threshold.by)
  mats <- switch(threshold.by,
                 raw=threshold_raw(A.norm, mat.thresh, ...),
                 mean=threshold_mean_sd(A.norm, mat.thresh, ...),
                 density=threshold_density(A.norm, mat.thresh, ...),
                 consistency=threshold_consistency(A.norm, mat.thresh, ...),
                 consensus=threshold_consensus(A.norm, mat.thresh, sub.thresh, inds,
                                               modality, algo, divisor, div.files, ...))

  A.norm.mean <- lapply(seq_along(mat.thresh), function(x)
                        abind(lapply(inds, function(y)
                                     rowMeans(mats$A.norm.sub[[x]][, , y], dims=2L)), along=3L))
  if (threshold.by == 'density') {
    A.norm.mean <- lapply(seq_along(mat.thresh), function(x)
                          array(apply(A.norm.mean[[x]], 3L, function(y)
                                      y * (y > get_thresholds(y, mat.thresh[x], emax))),
                                dim=dim(A.norm.mean[[x]])))
  }

  c(list(A=A), mats, list(A.norm.mean=A.norm.mean))
}

#-------------------------------------------------------------------------------
# Thresholding helper functions
#-------------------------------------------------------------------------------

#' Raw thresholding of a 3D array
#'
#' \code{threshold_raw} applies the input thresholds directly to the array.
#'
#' @return A list with elements
#'   \item{A.inds}{List of (binary) matrices with entries equal to 1 if the
#'     values in the input array (at that row-column) should be kept}
#'   \item{A.norm.sub}{List of 3D arrays of the input array with entries
#'     \dQuote{zeroed out} if they exceed the given threshold values. The list
#'     length will equal the length of \code{mat.thresh}}
#' @noRd

threshold_raw <- function(A, mat.thresh, ...) {
  A.norm.sub <- vector('list', length(mat.thresh))
  for (x in seq_along(mat.thresh)) {
    A.norm.sub[[x]] <- A * (A > mat.thresh[x])
    A.norm.sub[[x]] <- symmetrize(A.norm.sub[[x]], ...)
  }
  list(A.norm=A, A.inds=NULL, A.norm.sub=A.norm.sub)
}

#' Threshold by the mean + 2SD value (across all subjects)
#'
#' This calculates the mean (across subjects) and adds 2 standard deviations of
#' a 3D array to be used for thresholding/filtering of connectivity matrices.
#'
#' @noRd

threshold_mean_sd <- function(A, mat.thresh, ...) {
  A.norm.sub <- vector('list', length(mat.thresh))
  # Threshold: mean + 2SD > mat.thresh
  dimA <- dim(A)
  all.mean <- rowMeans(A, dims=2L)
  all.means <- array(all.mean, dim=dimA)
  all.sums <- rowSums((A - all.means)^2, dims=2L)
  all.sd <- sqrt(all.sums / (dimA[3L] - 1L))
  all.thresh <- array(all.mean + (2 * all.sd), dim=dimA)
  for (x in seq_along(mat.thresh)) {
    A.norm.sub[[x]] <- A * (all.thresh > mat.thresh[x])
    A.norm.sub[[x]] <- symmetrize(A.norm.sub[[x]], ...)
  }
  list(A.norm=A, A.inds=NULL, A.norm.sub=A.norm.sub)
}

#' Apply thresholds that will result in the given densities for all subjects
#'
#' @noRd

threshold_density <- function(A, mat.thresh, ...) {
  stopifnot(all(mat.thresh >= 0), all(mat.thresh <= 1))
  A.norm.sub <- vector('list', length(mat.thresh))
  nz <- dim(A)[3L]
  Asym <- symmetrize(A, ...)
  thrMat <- vapply(mat.thresh, function(x) apply(Asym, 3L, get_thresholds, x), numeric(nz))
  Aperm <- aperm(Asym, c(3L, 1L, 2L))
  for (x in seq_along(mat.thresh)) {
    A.norm.sub[[x]] <- Aperm * (Aperm > thrMat[, x])
    A.norm.sub[[x]] <- aperm(A.norm.sub[[x]], c(2L, 3L, 1L))
  }
  list(A.norm=A, A.inds=NULL, A.norm.sub=A.norm.sub)
}

#' Helper function to threshold based on consistency
#'
#' Creates thresholded 3D arrays using the \emph{coefficient of variation} and
#' user-defined thresholds.
#'
#' @noRd

threshold_consistency <- function(A, mat.thresh, ...) {
  stopifnot(all(mat.thresh >= 0), all(mat.thresh <= 1))
  A.norm.sub <- vector('list', length(mat.thresh))
  dimA <- dim(A)

  all.cv <- coeff_var(A)
  all.cv <- symmetrize(all.cv, 'min')
  threshes <- get_thresholds(all.cv, mat.thresh, decreasing=TRUE)
  A.inds <- lapply(threshes, function(x) ifelse(all.cv < x, 1L, 0L))
  for (z in seq_along(A.inds)) {
    tmp <- array(A.inds[[z]], dim=dimA)
    A.norm.sub[[z]] <- A * (tmp == 1L)
    A.norm.sub[[z]] <- symmetrize(A.norm.sub[[z]], ...)
  }
  list(A.norm=A, A.inds=A.inds, A.norm.sub=A.norm.sub)
}

#' Threshold based on group consensus
#'
#' Threshold data based on group \dQuote{consensus}; i.e., given a threshold for
#' the number/percent of subjects (\code{sub.thresh}) who exceed the
#' \code{mat.thresh} values, then keep those entries and zero out the rest.
#' @noRd

threshold_consensus <- function(A, mat.thresh, sub.thresh, inds, modality,
                                algo, divisor, div.files, ...) {
  A.bin <- A.bin.sums <- A.inds <- A.norm.sub <- vector('list', length(mat.thresh))
  # Binarize the array, then keep entries w/ >= "sub.thresh"% for each group
  for (x in seq_along(mat.thresh)) {
    A.bin[[x]] <- (A > mat.thresh[x]) + 0L
    A.bin.sums[[x]] <- lapply(inds, function(z)
                              rowSums(A.bin[[x]][, , z], dims=2L))
    A.bin.sums[[x]] <- abind(A.bin.sums[[x]], along=3L)
  }

  # For deterministic, threshold by size *after* binarizing
  if (modality == 'dti' && algo == 'deterministic' && divisor == 'size') {
    A <- normalize_mats(A, divisor, div.files, P=1)
  }

  # This is a list (# mat.thresh) of Nv x Nv x Ng arrays
  kNumSubjs <- lengths(inds)
  if (sub.thresh == 0) {
    A.inds <- lapply(seq_along(mat.thresh), function(x)
                     ifelse(A.bin.sums[[x]] > 0, 1L, 0L))
  } else {
    thresh <- sub.thresh * kNumSubjs
    for (x in seq_along(mat.thresh)) {
      tmp <- aperm(A.bin.sums[[x]], c(3L, 1L, 2L))
      A.inds[[x]] <- array(0L, dim=dim(tmp))
      A.inds[[x]][tmp >= thresh] <- 1L
      A.inds[[x]] <- aperm(A.inds[[x]], c(2L, 3L, 1L))
    }
  }

  # Back to a list of arrays for all subjects
  dimA <- dim(A); Nv <- dimA[1L]
  for (x in seq_along(mat.thresh)) {
    A.norm.sub[[x]] <- array(0, dim=dimA)
    keep <- array(FALSE, dim=dimA)
    for (y in seq_along(inds)) {
      tmp <- A.inds[[x]][, , y] == 1L
      keep[, , inds[[y]]] <- array(tmp, dim=c(Nv, Nv, kNumSubjs[y]))
    }
    A.norm.sub[[x]][keep] <- A[keep]
    A.norm.sub[[x]] <- symmetrize(A.norm.sub[[x]], ...)
  }
  list(A.norm=A, A.bin=A.bin, A.bin.sums=A.bin.sums, A.inds=A.inds, A.norm.sub=A.norm.sub)
}

#' Read in a list of text files as a numeric array
#'
#' \code{read.array} reads in a list of text files containing numeric matrices
#' and returns a numeric array. The order along the 3rd dimension will be the
#' same as the order of \code{infiles}.
#'
#' @param infiles Character vector of filenames
#' @param ncols Integer specifying the number of columns. By default, it is
#'   assumed the files contain \emph{square} matrices, but this behavior can be
#'   overridden through this argument
#' @noRd

read.array <- function(infiles, ncols=NULL) {
  Nv <- length(readLines(infiles[1L]))
  if (is.null(ncols)) ncols <- Nv
  array(sapply(infiles, function(x)
               matrix(scan(x, what=numeric(), n=Nv*ncols, quiet=TRUE),
                      Nv, ncols, byrow=TRUE)),
        dim=c(Nv, ncols, length(infiles)))
}

#' Normalize matrices in an array based on values in other matrices/vectors
#'
#' \code{normalize_mats} will normalize the matrices in \code{A} based on values
#' in separate matrices (found in the files specified by \code{div.files}) or
#' based on the matrix's row sums (if \code{divisor='rowSums'}).
#'
#' @param A Numeric array
#' @param divisor Character string. Either \code{rowSums}, \code{waytotal}, or
#'   \code{size}. If \code{none}, then it returns \code{A}.
#' @param div.files Character vector of filenames; required unless
#'   \code{divisor='rowSums'} or \code{divisor='none'}
#' @param P Integer representing the number of samples (in probabilistic
#'   tractography). Only used if \code{divisor='size'}.
#' @noRd

normalize_mats <- function(A, divisor, div.files, P) {
  dimA <- dim(A)
  if (divisor == 'none') {
    A.norm <- A
  } else if (divisor == 'rowSums') {
    A.norm <- array(apply(A, 3L, function(x) x / rowSums(x)), dim=dimA)

  } else {
    stopifnot(!is.null(div.files))
    div <- read.array(div.files, ncols=1L)

    if (divisor == 'waytotal') {
      # Control for streamline count by waytotal
      W <- array(apply(div, 3L, function(x) x[, rep.int(1L, dimA[1L])]), dim=dimA)
      A.norm <- A / W

    } else if (divisor == 'size') {
      # Control for the size (# voxels) of both regions 'x' and 'y'
      R <- array(apply(div, 3L, function(x) outer(x, x, FUN='+')), dim=dimA)
      A.norm <- 2 * A / (P * R)
    }
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
#' @param group.mats List (length equal to number of thresholds) of numeric
#'   arrays (3-dim) for group-level data
#' @param W.files Character vector of the filenames of the files with
#'   connectivity matrices
#' @inheritParams create_mats
#' @export
#'
#' @return List containing:
#' \item{W}{A 3-d array of the raw connection matrices}
#' \item{W.norm.sub}{List of 3-d arrays of the normalized connection matrices
#'   for all given thresholds}
#' \item{W.norm.mean}{List of 3-d arrays of the normalized connection matrices
#'   averaged for each group}
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#'   W.mats <- apply_thresholds(A.norm.sub, A.norm.mean, f.W, inds)
#' }

apply_thresholds <- function(sub.mats, group.mats, W.files, inds) {
  if (is.list(W.files) || length(W.files) == 1L) W.files <- dir2files(W.files)
  dimA <- dim(sub.mats[[1L]])
  stopifnot(length(W.files) == dimA[3L])
  W <- read.array(W.files)
  W.norm.sub <- W.norm.mean <- vector('list', length(sub.mats))
  for (x in seq_along(sub.mats)) {
    W.norm.sub[[x]] <- array(0, dim=dimA)
    keep <- sub.mats[[x]] > 0
    W.norm.sub[[x]][keep] <- W[keep]
  }
  ng <- dim(group.mats[[1L]])[3L]
  W.norm.mean <- lapply(seq_along(group.mats), function(x)
                        lapply(seq_len(ng), function(g)
                               ifelse(group.mats[[x]][, , g] > 0,
                                      rowMeans(W.norm.sub[[x]][, , inds[[g]]], dims=2L),
                                      0)))
  W.norm.mean <- lapply(W.norm.mean, function(x) abind(x, along=3L))
  return(list(W=W, W.norm.sub=W.norm.sub, W.norm.mean=W.norm.mean))
}

#' Return a vector of filenames based on a directory name or options list
#'
#' @param x Either a character string specifying a directory or a list of
#'   arguments to be passed to \code{\link{list.files}}
#' @keywords internal

dir2files <- function(x) {
  if (length(x) > 1L && !is.list(x)) return(x)  # Vector of filenames
  fargs <- formals(list.files)
  if (is.list(x)) {
    if (is.null(attr(x, 'names'))) names(x) <- names(fargs)[seq_along(x)]
    if (!hasName(x, 'full.names')) x <- c(x, full.names=TRUE)
  } else {
    x <- list(path=x, full.names=TRUE)
  }
  x <- c(x, fargs[setdiff(names(fargs), names(x))])
  stopifnot(file_test('-d', x$path))
  do.call(list.files, x)
}
