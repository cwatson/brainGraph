#' Coordinates for data from the AAL-based atlases
#'
#' Datasets containing spatial coordinates for the original AAL atlases and the
#' newer AAL2 atlases, along with indices for the major lobes and hemispheres of
#' the brain.
#'
#' @format A data frame with 90 or 116 (for the original AAL atlases), or 94 or
#' 120 (for the newer AAL2 atlases) observations on the following 7 variables:
#' \describe{
#'   \item{\code{name}}{a character vector of region names}
#'   \item{\code{x.mni}}{a numeric vector of x-coordinates (in MNI space)}
#'   \item{\code{y.mni}}{a numeric vector of y-coordinates (in MNI space)}
#'   \item{\code{z.mni}}{a numeric vector of z-coordinates (in MNI space)}
#'   \item{\code{lobe}}{a factor with levels \code{Frontal} \code{Parietal}
#'     \code{Temporal} \code{Occipital} \code{Insula} \code{Limbic} \code{SCGM}
#'     and \code{Cerebellum} (for \code{aal116} and \code{aal2.120})}
#'   \item{\code{hemi}}{a factor with levels \code{L} \code{R}}
#'   \item{\code{index}}{a numeric vector}
#' }
#'
#' @name AAL
#' @aliases aal116
#' @rdname aal
#' @references Tzourio-Mazoyer, N. and Landeau, B. and Papathanassiou, D. and
#'   Crivello, F. and Etard, O. and Delcroix, N. and Mazoyer, B. and Joliot, M.
#'   (2002) Automated anatomical labeling of activations in SPM using a
#'   macroscopic anatomical parcellation of the MNI MRI single-subject brain.
#'   \emph{NeuroImage}, \bold{15(1)}, 273--289.
#'   \url{https://dx.doi.org/10.1006/nimg.2001.0978}
#' @references Rolls, E.T. and Joliot, M. and Tzourio-Mazoyer, N. (2015)
#'   Implementation of a new parcellation of the orbitofrontal cortex in the
#'   automated anatomical labelling atlas. \emph{NeuroImage}, \bold{122}, 1--5.
#'   \url{https://dx.doi.org/10.1016/j.neuroimage.2015.07.075}
"aal116"

#' @aliases aal90
#' @rdname aal
"aal90"

#' @aliases aal2.120
#' @rdname aal
"aal2.120"

#' @aliases aal2.94
#' @rdname aal
"aal2.94"

#' Coordinates for data from Freesurfer atlases
#'
#' Datasets containing spatial coordinates for the Freesurfer atlases:
#' Destrieux, Desikan-Killiany (DK), and Desikan-Killiany-Tourville (DKT). The
#' datasets also contain indices for the major lobes and hemispheres of the
#' brain, in addition to the \emph{class} variable for Destrieux atlases.
#'
#' @format A data frame with 148 or 162 (for Destrieux), 68 or 82 (for DK), or
#' 62 or 76 (for DKT) observations on the following 8 variables:
#' \describe{
#'   \item{\code{name}}{a character vector of region names}
#'   \item{\code{x.mni}}{a numeric vector of x-coordinates (in MNI space)}
#'   \item{\code{y.mni}}{a numeric vector of y-coordinates (in MNI space)}
#'   \item{\code{z.mni}}{a numeric vector of z-coordinates (in MNI space)}
#'   \item{\code{lobe}}{a factor with levels \code{Frontal}, \code{Parietal},
#'     \code{Temporal}, \code{Occipital}, \code{Insula}, \code{Limbic}, and
#'     \code{SCGM} (for atlases ending in \code{.scgm})}
#'   \item{\code{hemi}}{a factor with levels \code{L} \code{R}}
#'   \item{\code{index}}{a numeric vector}
#'   \item{\code{name.full}}{a character vector of full region names, for the DK
#'     and DKT atlases}
#'   \item{\code{class}}{a factor with levels \code{G} \code{G_and_S} \code{S}}
#' }
#'
#' @name FreesurferAtlases
#' @aliases destrieux destrieux.scgm dk dk.scgm dkt dkt.scgm
#' @rdname freesurfer_atlases
#' @references Destrieux, C. and Fischl, B. and Dale, A. and Halgren E. (2010)
#'   Automatic parcellation of human cortical gyri and sulci using standard
#'   anatomic nomenclature. \emph{NeuroImage}, \bold{53(1)}, 1--15.
#'   \url{https://dx.doi.org/10.1016/j.neuroimage.2010.06.010}
#' @references Desikan, R.S. and Segonne, F. and Fischl, B. et al. (2006) An
#'   automated labeling system for subdividing the human cerebral cortex on MRI
#'   scans into gyral based regions of interest. \emph{NeuroImage}, \bold{31},
#'   968--980. \url{https://dx.doi.org/10.1016/j.neuroimage.2006.01.021}
#' @references Klein, A. and Tourville, J. (2012) 101 labeled brain images
#'   and a consistent human cortical labeling protocol. \emph{Front Neurosci},
#'   6. \url{https://dx.doi.org/10.3389/fnins.2012.00171}
"destrieux"

#' @rdname freesurfer_atlases
"destrieux.scgm"

#' @rdname freesurfer_atlases
"dk"

#' @rdname freesurfer_atlases
"dk.scgm"

#' @rdname freesurfer_atlases
"dkt"

#' @rdname freesurfer_atlases
"dkt.scgm"
