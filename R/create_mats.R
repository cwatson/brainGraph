#' Create connection matrices for tractography or fMRI data
#'
#' \code{create_mats} will take a vector of filenames which contain connection
#' matrices (e.g. the \emph{fdt_network_matrix} files from FSL or the
#' \emph{ROICorrelation.txt} files from DPABI) and create arrays of this data.
#' You may choose to normalize these matrices by the \emph{waytotal} or
#' \emph{region size} (tractography), or not at all.
#'
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
#' }

create_mats <- function(A.files, modality=c('dti', 'fmri'),
                        divisor=c('none', 'waytotal', 'size', 'rowSums'),
                        div.files=NULL,
                        threshold.by=c('consensus', 'density', 'mean', 'consistency', 'raw'),
                        mat.thresh=0, sub.thresh=0.5, inds=list(1:length(A.files)),
                        algo=c('probabilistic', 'deterministic'), P=5e3, ...) {

  # Argument checking
  #-----------------------------------------------------------------------------
  kNumSubjs <- lengths(inds)
  stopifnot(all(file.exists(A.files)),
            sum(kNumSubjs) == length(A.files),
            sub.thresh >= 0 && sub.thresh <= 1)
  if (!is.null(div.files) && any(!file.exists(div.files))) {
    stop("'div.files' are not all valid files! Please check your data.")
  }
  A.bin <- A.bin.sums <- A.inds <- NULL

  A <- read.array(A.files)
  Nv <- nrow(A)
  nz <- sum(kNumSubjs)
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
    A.norm.sub <- vector('list', length(mat.thresh))
    for (x in seq_along(mat.thresh)) {
      A.norm.sub[[x]] <- array(apply(A.norm, 3, function(z) z * (z > mat.thresh[x])), dim=dim(A.norm))
    }

  } else if (threshold.by %in% c('density', 'consistency')) {
    stopifnot(all(mat.thresh >= 0) && all(mat.thresh <= 1))
    A.norm.sub <- vector('list', length(mat.thresh))

    if (threshold.by == 'density') {
      Asym <- symmetrize_array(A.norm, ...)
      thrMat <- vapply(mat.thresh, function(x) apply(Asym, 3, get_thresholds, x, emax), numeric(nz))
      for (x in seq_along(mat.thresh)) {
        A.norm.sub[[x]] <- array(0, dim=dim(Asym))
        for (k in seq_len(nz)) {
          A.norm.sub[[x]][, , k] <- Asym[, , k] * (Asym[, , k] > thrMat[k, x])
        }
      }
    } else if (threshold.by == 'consistency') {
      all.cv <- apply(A.norm, 1:2, coeff_var)
      all.cv <- symmetrize_mats(all.cv, 'min')
      threshes <- get_thresholds(all.cv, mat.thresh, emax, decreasing=TRUE)
      A.inds <- lapply(threshes, function(x) ifelse(all.cv < x, 1L, 0L))
      for (z in seq_along(A.inds)) {
        A.norm.sub[[z]] <- array(apply(A.norm, 3, `*`, (A.inds[[z]] == 1L)), dim=dim(A.norm))
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
      if (sub.thresh == 0) {
        A.inds <- lapply(seq_along(mat.thresh), function(x)
                         abind(lapply(seq_along(inds), function(z)
                                      ifelse(A.bin.sums[[x]][, , z] > 0, 1L, 0L)),
                               along=3))
      } else {
        A.inds <- lapply(seq_along(mat.thresh), function(x)
                         abind(lapply(seq_along(inds), function(z)
                                      ifelse(A.bin.sums[[x]][, , z] >= sub.thresh * kNumSubjs[z], 1L, 0L)),
                               along=3))
      }

      # Back to a list of arrays for all subjects
      A.norm.sub <-
        lapply(seq_along(mat.thresh), function(x)
               abind(lapply(seq_along(inds), function(y)
                            array(sapply(inds[[y]], function(z)
                                         ifelse(A.inds[[x]][, , y] == 1, A.norm[, , z], 0)),
                                  dim=dim(A.norm[, , inds[[y]]]))), along=3))

      # Re-order A.norm.sub so that it matches the input files, A, A.norm, etc.
      for (x in seq_along(mat.thresh)) {
        tmp <- array(0, dim=dim(A.norm.sub[[x]]))
        tmp[, , unlist(inds)] <- A.norm.sub[[x]]
        A.norm.sub[[x]] <- tmp
      }

    } else if (threshold.by == 'mean') {
      # Threshold: mean + 2SD > mat.thresh
      #---------------------------------
      all.mean <- rowMeans(A.norm, dims=2)
      all.means <- replicate(nz, all.mean)
      all.sums <- rowSums((A.norm - all.means)^2, dims=2)
      all.sd <- sqrt(all.sums / (nz - 1))
      all.thresh <- all.mean + (2 * all.sd)

      A.norm.sub <-
        lapply(mat.thresh, function(z)
               array(apply(A.norm, 3, function(x) x * (all.thresh > z)), dim=dim(A.norm)))
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
#' off-diagonal elements values of either the \code{max}, \code{min}, or
#' \code{average} of \eqn{\{A(i, j), A(j, i)\}}. The default is \code{max}
#' because that is the default for
#' \code{\link[igraph]{graph_from_adjacency_matrix}}.
#'
#' @param A Numeric matrix
#' @param symm.by Character string; how to create symmetric off-diagonal
#'   elements (default: \code{max})
#' @export
#' @return Either a single symmetrized matrix, or an (3D) array
#'
#' @family Matrix functions
#' @seealso \code{\link[igraph]{graph_from_adjacency_matrix}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

symmetrize_mats <- function(A, symm.by=c('max', 'min', 'avg')) {
  stopifnot(nrow(A) == ncol(A))

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
    Nv <- nrow(A)
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
#' \code{apply_thresholds} will threshold an additional set of matrices (e.g.,
#' FA-weighted matrices for DTI tractography) based on the matrices that have
#' been returned from \code{\link{create_mats}}. This ensures that the same
#' connections are present in both sets of matrices.
#'
#' @param sub.mats List (length equal to number of thresholds) of numeric arrays
#'   (3-dim) for all subjects
#' @param group.mats List (equal to number of thresholds) of lists (equal to
#'   number of groups) of numeric matrices for group-level data
#' @param W.files Character vector of the filenames of the files containing your
#'   connectivity matrices
#' @param inds List (length equal to number of groups) of integers; each list
#'   element should be an integer vector of length equal to the group sizes
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
  W <- read.array(W.files)
  W.norm.sub <- lapply(sub.mats, function(x)
                       array(sapply(seq_len(dim(x)[3]), function(y)
                                    ifelse(x[, , y] > 0, W[, , y], 0)),
                             dim=dim(x)))
  W.norm.mean <- lapply(seq_along(group.mats), function(x)
                        lapply(seq_len(dim(group.mats[[x]])[3]), function(y)
                               ifelse(group.mats[[x]][, , y] > 0,
                                      rowMeans(W.norm.sub[[x]][, , inds[[y]]], dims=2),
                                      0)))
  W.norm.mean <- lapply(W.norm.mean, function(x) abind(x, along=3))
  return(list(W=W, W.norm.sub=W.norm.sub, W.norm.mean=W.norm.mean))
}
