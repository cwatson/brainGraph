#' Subset observations of a bg_GLM object
#'
#' The \code{[} method allows you to select observations (i.e., rows of \code{X}
#' and \code{y}) and independent variables (i.e., columns of \code{X}) from a
#' \code{bg_GLM} object.
#'
#' @note The \code{[} method is used when calculating \emph{studentized
#' residuals} and other \dQuote{leave-one-out} diagnostics, and typically should
#' not be called directly by the user.
#'
#' @param i Integer/character vector; the observation number(s) or row names to
#'   select or remove
#' @param j Integer/character vector; the design matrix column number(s) or
#'   names to select or remove
#' @export
#' @return A \code{bg_GLM} object with the specified row(s) selected or removed
#'   from both \code{X} and \code{y}, and column(s) selected/removed from
#'   \code{X}
#' @rdname glm

`[.bg_GLM` <- function(x, i, j) {
  dimX <- dim(x$X)
  if (missing(i)) i <- seq_len(dimX[1L])
  if (missing(j)) j <- seq_len(dimX[2L])
  x$X <- if (length(dimX) == 3L) x$X[i, j, , drop=FALSE] else x$X[i, j, drop=FALSE]
  x$y <- x$y[i, , drop=FALSE]
  return(x)
}

#' Extract basic information from a bg_GLM object
#'
#' These functions return the \code{terms}, \emph{term labels}, \emph{model
#' formula}, \dQuote{case names}, \dQuote{variable names}, \emph{region names},
#' and number of observations for a \code{bg_GLM} object. The term labels are
#' used for ANOVA tables.
#'
#' @note \code{terms} has only been tested for 2-way interactions. If your model
#' has higher-order interactions, it may not work properly. Functions affected
#' would include: \code{\link{residuals.bg_GLM}} (partial),
#' \code{\link{vif.bg_GLM}}, and \code{\link{anova.bg_GLM}}. Furthermore,
#' \code{labels} and \code{formula} may also be incorrect.
#' @note \code{formula} returns a character string, not a \code{formula}
#' object.
#'
#' @param x,object A \code{bg_GLM} object
#' @param ... Unused
#' @export
#' @return \code{terms} returns a named integer list in which the names are the
#'   term labels and the list elements are the column(s) of the design matrix
#'   for each term. \code{nobs} returns an integer. The other functions return
#'   character vectors.
#' @name GLM basic info
#' @rdname glm_info

nobs.bg_GLM <- function(object, ...) {
  dim(object$X)[1L]
}

#' @export
#' @rdname glm_info

terms.bg_GLM <- function(x, ...) {
  cvnames <- names(x$covars)[-1L]
  vnames <- variable.names(x)
  if (x$outcome != x$measure && !(x$measure %in% cvnames)) {
    cvnames <- c(cvnames, x$measure)
  }
  cols <- sapply(cvnames, grep, vnames)

  # Determine if there are any interaction terms; TODO: only works for 2-way interactions
  if (any(grepl(':', vnames))) {
    combos <- combn(seq_along(cols), 2)
    ints <- apply(combos, 2, function(x) intersect(cols[[x[1L]]], cols[[x[2L]]]))
    match <- combos[, which(lengths(ints) > 0)]
    int_terms <- cvnames[match]
    int_name <- paste(int_terms, collapse=':')
    nonzero <- Filter(length, ints)[[1L]]
    cols[[int_name]] <- nonzero
    cols[int_terms] <- lapply(cols[int_terms], setdiff, nonzero)
  }
  if (any(grepl('Intercept', vnames))) {
    cols <- c(list(Intercept=grep('Intercept', vnames)), cols)
  }
  cols <- cols[order(unlist(lapply(cols, `[[`, 1)))]
  return(cols)
}

#' @export
#' @rdname glm_info

formula.bg_GLM <- function(x, ...) {
  tlabels <- labels(x)
  if (any(grepl('Intercept', tlabels))) tlabels <- tlabels[-grep('Intercept', tlabels)]
  splits <- strsplit(tlabels, ':')
  splits.int <- which(lengths(splits) > 1L)
  for (i in splits.int) {
    if (all(splits[[i]] %in% tlabels)) {
      inds.remove <- vapply(splits[[i]], function(x) which(x == tlabels), integer(1))
      inds.remove <- c(inds.remove, i)
      tlabels <- tlabels[-inds.remove]
      tlabels <- c(tlabels, paste(splits[[i]], collapse=' * '))
    }
  }
  form <- paste(x$outcome, '~', paste(tlabels, collapse=' + '))
  return(form)
}

#' @export
#' @rdname glm_info

labels.bg_GLM <- function(object, ...) {
  labels(terms(object))
}

#' @export
#' @method case.names bg_GLM
#' @rdname glm_info

case.names.bg_GLM <- function(object, ...) {
  dimnames(object$X)[[1L]]
}

#' @export
#' @method variable.names bg_GLM
#' @rdname glm_info

variable.names.bg_GLM <- function(object, ...) {
  dimnames(object$X)[[2L]]
}

#' @export
#' @method region.names bg_GLM
#' @rdname glm_info

region.names.bg_GLM <- function(object) {
  region <- NULL
  object$DT[, levels(region)]
}

#' @export
#' @rdname glm_info

nregions.bg_GLM <- function(object) {
  region <- NULL
  object$DT[, nlevels(region)]
}
