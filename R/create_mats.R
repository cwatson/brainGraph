#' Create connection matrices for tractography or fMRI data
#'
#' This function will take a vector of filenames which contain connection
#' matrices (e.g. the \emph{fdt_network_matrix} files from FSL or the
#' \emph{ROICorrelation.txt} files from DPABI) and create arrays of this data.
#' You may choose to normalize these matrices by the \emph{waytotal} or
#' \emph{region size} (tractography), or not at all (fMRI).
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
#'   choose \code{density} if you want all resulting matrices in a subject group
#'   to have the same density (default: \code{raw})
#' @param mat.thresh Numeric (vector) for thresholding connection matrices
#'   (default: 0)
#' @param sub.thresh Numeric (between 0 and 1) for thresholding by subject
#'   numbers (default: 0.5)
#' @param inds List (length equal to number of groups) of integers; each list
#'   element should be a vector of length equal to the group sizes
#' @param algo Character string of the tractography algorithm used (default:
#'   \code{'probabilistic'})
#' @param P Integer; number of samples generated using FSL (default: 5000)
#' @export
#'
#' @return A list containing:
#' \item{A}{A 3-d array of the raw connection matrices}
#' \item{A.norm}{A 3-d array of the normalized connection matrices}
#' \item{A.bin}{A 3-d array of binarized connection matrices}
#' \item{A.bin.sums}{A list of 2-d arrays of connection matrices, with each
#'   entry signifying the number of subjects with a connection present; the
#'   number of list elements equals the length of \code{mat.thresh}}
#' \item{A.inds}{A list of arrays of binarized connection matrices, containing 1
#'   if that entry is to be included}
#' \item{A.norm.sub}{A list of 3-d arrays of the normalized connection matrices
#'   for all given thresholds}
#' \item{A.norm.mean}{A list of 2-d arrays of the normalized connection matrices
#'   averaged for each group}
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' thresholds <- seq(from=0.001, to=0.01, by=0.001)
#' my.mats <- create_mats(f.A, modality='fmri', threshold.by='raw',
#'   mat.thresh=thresholds, sub.thresh=0.5, inds=inds)
#' dti.mats <- create_mats(f.A, divisor='waytotal', div.files=f.way,
#'   mat.thresh=thresholds, sub.thresh=0.5, inds=inds)
#' }

create_mats <- function(A.files, modality=c('dti', 'fmri'),
                        divisor=c('none', 'waytotal', 'size', 'rowSums'),
                        div.files=NULL, threshold.by=c('raw', 'density'),
                        mat.thresh=0, sub.thresh=0.5, inds,
                        algo=c('probabilistic', 'deterministic'), P=5e3) {

  kNumSubjs <- lengths(inds)
  stopifnot(isTRUE(all(sapply(A.files, file.exists))),
            isTRUE(all(sapply(div.files, file.exists))),
            sum(kNumSubjs) == length(A.files),
            sub.thresh >= 0 && sub.thresh <= 1)
  A.bin <- A.bin.sums <- A.inds <- NULL

  Nv <- length(readLines(A.files[1]))
  A <- array(sapply(A.files, function(x)
                    matrix(scan(x, what=numeric(0), n=Nv*Nv, quiet=T),
                           Nv, Nv, byrow=T)),
             dim=c(Nv, Nv, sum(kNumSubjs)))

  A[is.nan(A)] <- 0
  A.norm <- A

  modality <- match.arg(modality)
  algo <- match.arg(algo)
  divisor <- match.arg(divisor)
  if (modality == 'dti' & algo == 'probabilistic' & divisor != 'none') {
    A.norm <- normalize_mats_prob(A, divisor, div.files, Nv, kNumSubjs, P)
    A.norm[is.nan(A.norm)] <- 0
  }

  threshold.by <- match.arg(threshold.by)
  if (threshold.by == 'raw') {
    # Binarize the array, then keep entries w/ >= "sub.thresh"% for each group
    # These are lists of arrays; each list element is a different threshold
    A.bin <- lapply(mat.thresh, function(x) (A.norm > x) + 0)
    A.bin.sums <- lapply(seq_along(mat.thresh), function(y)
                         lapply(inds, function(x)
                                rowSums(A.bin[[y]][, , x], dims=2)))

    # In case deterministic tractography was used and the user would like to
    # adjust for ROI size
    if (modality == 'dti' & algo == 'deterministic' & divisor == 'size') {
      A.norm <- normalize_mats_det(A.norm, div.files, Nv, kNumSubjs)
    }

    # This is a list (# mat.thresh) of lists (# groups) of the Nv x Nv group matrix
    if (sub.thresh == 0) {
      A.inds <- lapply(seq_along(mat.thresh), function(y)
                       lapply(seq_along(inds), function(x)
                              ifelse(A.bin.sums[[y]][[x]] > 0, 1, 0)))
    } else {
      A.inds <- lapply(seq_along(mat.thresh), function(y)
                       lapply(seq_along(inds), function(x)
                              ifelse(A.bin.sums[[y]][[x]] >= sub.thresh * kNumSubjs[x],
                                     1,
                                     0)))
    }

    # Back to a list of arrays for all subjects
    A.norm.sub <-
      lapply(seq_along(mat.thresh), function(z)
             lapply(seq_along(inds), function(x)
                    array(sapply(inds[[x]], function(y)
                                 ifelse(A.inds[[z]][[x]] == 1,
                                        A.norm[, , y],
                                        0)),
                          dim=dim(A.norm[, , inds[[x]]]))))
    A.norm.sub <- lapply(A.norm.sub, function(x) do.call(abind::abind, x))

  } else if (threshold.by == 'density') {
    stopifnot(all(mat.thresh >= 0) && all(mat.thresh <= 1))
    emax <- Nv * (Nv - 1) / 2

    A.norm.sub <-
      lapply(mat.thresh, function(x)
             array(apply(A, 3, function(y) {
                           thresh <- sort(y[lower.tri(y)])[emax - x * emax]
                           ifelse(y > thresh, y, 0)
                         }), dim=dim(A)))
#    A.norm.sub <- lapply(A.bin, function(x)
#                    array(apply(x, 3, function(y)
#                                ifelse(y == 1, y, 0)),
#                          dim=dim(x)))
  }


  # Group means; first take the means of the unthresholded data
  #A.norm.mean <- lapply(inds, function(x) apply(A.norm[, , x], c(1, 2), mean))
  #A.group <- lapply(seq_len(length(inds)), function(x)
  #                  lapply(mat.thresh, function(y) (A.norm.mean[[x]] > y) + 0))
  A.norm.mean <- lapply(seq_along(mat.thresh), function(x)
                        lapply(inds, function(y)
                               rowMeans(A.norm.sub[[x]][, , y], dims=2)))

  return(list(A=A, A.norm=A.norm, A.bin=A.bin, A.bin.sums=A.bin.sums,
              A.inds=A.inds, A.norm.sub=A.norm.sub, A.norm.mean=A.norm.mean))#,
              #A.group=A.group))
}

normalize_mats_prob <- function(A, divisor, div.files, Nv, kNumSubjs, P) {
  div <- array(sapply(div.files, function(x)
                      matrix(scan(x, what=numeric(0), n=Nv*1, quiet=T),
                             Nv, 1, byrow=T)),
               dim=c(Nv, 1, sum(kNumSubjs)))

  if (divisor == 'waytotal') {
    # Control for streamline count by waytotal
    W <- array(apply(div, 3, function(x) x[, rep(1, Nv)]), dim=dim(A))
    A.norm <- A / W

  } else if (divisor == 'size') {
    # Control for the size (# voxels) of both regions 'x' and 'y'
    R <- array(apply(div, 3, function(x)
                     cbind(sapply(seq_len(Nv), function(y) x + x[y]))),
               dim=dim(A))
    A.norm <- 2 * A / (P * R)

  } else if (divisor == 'rowSums') {
    A.norm <- array(apply(A, 3, function(x) x / rowSums(x)), dim=dim(A))
  }
  return(A.norm)
}

normalize_mats_det <- function(A.norm, div.files, Nv, kNumSubjs) {
  div <- array(sapply(div.files, function(x)
                      matrix(scan(x, what=numeric(0), n=Nv*1, quiet=T),
                             Nv, 1, byrow=T)),
               dim=c(Nv, 1, sum(kNumSubjs)))
  R <- array(apply(div, 3, function(x)
                   cbind(sapply(seq_len(Nv), function(y) x + x[y]))), dim=dim(A.norm))
  A.norm <- 2 * A.norm / R
  return(A.norm)
}
