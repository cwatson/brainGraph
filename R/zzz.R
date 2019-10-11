#' Default options for brainGraph
#'
#' brainGraph is a package for performing \emph{graph theory analysis} of brain
#' MRI data.
#'
#' @section Package options:
#'
#' brainGraph uses the following \code{\link{options}} to configure behavior:
#' \itemize{
#'   \item \code{bg.subject_id}: character string specifying the name your
#'     project/study uses as a subject identifier. All imported data (e.g.,
#'     covariates tables) \emph{MUST} have a column matching this. One possible
#'     alternative is \code{'participant_id'}, recommended by BIDS. Default:
#'     \code{'Study.ID'}
#'   \item \code{bg.group}: character string specifying the name your
#'     project/study uses as a group identifier. All imported data (e.g.,
#'     covariates tables) \emph{MUST} have a column matching this. One possible
#'     alternative is \code{'group'}, recommended by BIDS. Default:
#'     \code{'Group'}
#'   \item \code{bg.session}: character string specifying the name your
#'     project/study uses as a \dQuote{time} or session identifier, in the case
#'     of longitudinal studies. All imported data (e.g., covariates tables)
#'     \emph{MUST} have a column matching this. One possible alternative is
#'     \code{'session_id'}, recommended by BIDS. Default: \code{'Time'}
#'   \item \code{bg.progress}: logical indicating whether to show progress bars
#'     for functions that provide the option. Default: \code{TRUE}
#'   \item \code{bg.ncpus}: integer indicating the number of cores to use for
#'     parallel operations. Only used if you have not already registered a
#'     parallel backend (see Chapter 5 of the User Guide or
#'     \url{https://github.com/cwatson/brainGraph/README.md} for examples).
#'     Default: \code{2L}
#' }
#' @docType package
#' @name brainGraph
#' @aliases brainGraph-options
NULL

# Default values for options
bg.options <- list(
  bg.subject_id='Study.ID',
  bg.group='Group',
  bg.session='Time',
  bg.progress=TRUE,
  bg.ncpus=2L
)

.onLoad <- function(libname, pkgname) {
  op <- options()
  toset <- !(names(bg.options) %in% names(op))
  if (any(toset)) options(bg.options[toset])

  invisible()
}
