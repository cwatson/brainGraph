#' Default options for brainGraph
#'
#' @section Package options:
#'
#' brainGraph uses the following \code{\link{options}} to configure behavior:
#' \itemize{
#'   \item \code{bg.subject.id}: character string specifying the name your
#'     project/study uses as a subject identifier. All imported data (e.g.,
#'     covariates tables) \emph{MUST} have a column matching this. One possible
#'     alternative is \code{'participant_id'}, recommended by BIDS. Default:
#'     \code{'Study.ID'}
#'   \item \code{bg.group.id}: character string specifying the name your
#'     project/study uses as a group identifier. All imported data (e.g.,
#'     covariates tables) \emph{MUST} have a column matching this. One possible
#'     alternative is \code{'group'}, recommended by BIDS. Default:
#'     \code{'Group'}
#'   \item \code{bg.session.id}: character string specifying the name your
#'     project/study uses as a "time" identifier, in the case of longitudinal
#'     studies. All imported data (e.g., covariates tables) \emph{MUST} have a
#'     column matching this. One possible alternative is \code{'session_id'},
#'     recommended by BIDS. Default: \code{'Time'}
#'   \item \code{bg.color.palette}: character string specifying the color
#'     palette to use for plots. This is only relevant if you are producing
#'     graph plots with colors for vertices and/or edges. The only choices are
#'     \code{'rgb'} (the default) and \code{'cmyk'} (which might be required by
#'     a journal)
#' }
#' @docType package
#' @name brainGraph_options
NULL

#TODO: think of any more that would help ("datadir"/"studydir"? "atlas"?)
bg_default_options <- list(
  bg.subject.id='Study.ID',
  bg.group.id='Group',  #TODO: would this really be necessary? maybe would be nice for foreign users?
  bg.session.id='Time',
  bg.color.palette='rgb'  #TODO: think more about this; see if it's possible to allow only 'rgb'/'cmyk', and what to do from there
)

.onLoad <- function(libname, pkgname) {
  op <- options()
  toset <- !(names(bg_default_options) %in% names(op))
  if (any(toset)) options(bg_default_options[toset])

  invisible()
}
