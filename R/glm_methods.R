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
#' @include brainGraph_GLM.R
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

nobs.bg_GLM <- function(object, ...) dim(object$X)[1L]

#' @export
#' @rdname glm_info

terms.bg_GLM <- function(x, ...) {
  cvnames <- setdiff(names(x$DT.Xy), c(getOption('bg.subject_id'), 'region', x$outcome))
  vnames <- variable.names(x)
  cols <- cols_orig <- sapply(cvnames, function(y) list(grep(y, vnames)))

  # Determine if there are any interaction terms
  rmatches <- gregexpr(':', vnames)
  numMatches <- vapply(rmatches, function(y) length(y[y > 0]), integer(1))
  if (any(numMatches > 0)) {
    combos <- combn(seq_along(cvnames), 2, simplify=FALSE)
    ints <- lapply(combos, function(y) intersect(cols[[y[1L]]], cols[[y[2L]]]))
    matches <- lapply(combos[which(lengths(ints) > 0L)], as.matrix)
    int_terms <- vapply(matches, function(y) cvnames[y], character(2))
    twoWays <- apply(int_terms, 2L, paste, collapse=':')
    nonzero <- Filter(length, ints)
    for (i in seq_along(twoWays)) {
      cols[[twoWays[i]]] <- nonzero[[i]]
      cols[int_terms[, i]] <- lapply(cols[int_terms[, i]], setdiff, nonzero[[i]])
    }
    cols <- Filter(length, cols)
    if (any(numMatches == 2)) {
      combos <- combn(seq_along(cvnames), 3)
      ints <- apply(combos, 2L, function(y) intersect(cols_orig[[y[1L]]], intersect(cols_orig[[y[2L]]], cols_orig[[y[3L]]])))
      matches <- as.matrix(combos[, which(lengths(ints) > 0L)])
      int_terms <- cvnames[matches]
      dim(int_terms) <- dim(matches)
      int_name <- apply(int_terms, 2L, paste, collapse=':')
      nonzero <- Filter(length, ints)
      for (i in seq_along(int_name)) {
        cols[[int_name[i]]] <- nonzero[[i]]
        cols[int_terms[, i]] <- lapply(cols[int_terms[, i]], setdiff, nonzero[[i]])
        cols[twoWays] <- lapply(cols[twoWays], setdiff, nonzero[[i]])
      }
    }
  }
  if (any(grepl('Intercept', vnames))) {
    cols <- c(list(Intercept=grep('Intercept', vnames)), cols)
  }
  cols <- cols[order(unlist(lapply(cols, `[[`, 1L)))]
  return(cols)
}

#' @export
#' @rdname glm_info

formula.bg_GLM <- function(x, ...) {
  tlabels <- labels(x)
  if (any(grepl('Intercept', tlabels))) tlabels <- tlabels[-grep('Intercept', tlabels)]
  rmatches <- gregexpr(':', tlabels)
  intOrder <- vapply(rmatches, function(y) length(y[y > 0]), integer(1)) + 1L
  splits <- strsplit(tlabels, ':')
  if (any(intOrder == 3L)) {
    for (i in which(intOrder == 3L)) {
      inds.remove <- vapply(splits[[i]], function(y) which(y == tlabels), integer(1))
      inds.remove <- c(inds.remove, i)
      twoWays <- combn(splits[[i]], 2L, function(y) paste(y, collapse=':'))
      inds.remove <- c(inds.remove, vapply(twoWays, function(y) which(y == tlabels), integer(1)))
      tlabels <- tlabels[-inds.remove]
      tlabels <- c(tlabels, paste(splits[[i]], collapse=' * '))
    }
  } else {
    for (i in which(intOrder > 1L)) {
      if (all(splits[[i]] %in% tlabels)) {
        inds.remove <- vapply(splits[[i]], function(y) which(y == tlabels), integer(1))
        inds.remove <- c(inds.remove, i)
        tlabels <- tlabels[-inds.remove]
        tlabels <- c(tlabels, paste(splits[[i]], collapse=' * '))
      }
    }
  }
  form <- paste(x$outcome, '~', paste(tlabels, collapse=' + '))
  return(form)
}

#' @export
#' @rdname glm_info
labels.bg_GLM <- function(object, ...) labels(terms(object))

#' @export
#' @method case.names bg_GLM
#' @rdname glm_info
case.names.bg_GLM <- function(object, ...) dimnames(object$X)[[1L]]

#' @export
#' @method variable.names bg_GLM
#' @rdname glm_info
variable.names.bg_GLM <- function(object, ...) dimnames(object$X)[[2L]]

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
