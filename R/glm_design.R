#' Create a design matrix for linear model analysis
#'
#' \code{brainGraph_GLM_design} takes a \code{data.table} of covariates and
#' returns a \emph{design matrix} to be used in linear model analysis.
#'
#' There are three different ways to code factors: \emph{dummy}, \emph{effects},
#' or \emph{cell-means} (chosen by the argument \code{coding}). \emph{Effects}
#' coding is sometimes referred to as \emph{deviation} coding. \emph{Dummy}
#' coding is the default when calling \code{lm}. To understand the difference
#' between these, see Chapter 8 of the User Guide.
#'
#' @section Character variables:
#' The default behavior is to convert all character columns (excluding the Study
#' ID column and any that you list in the \code{binarize} argument) to factor
#' variables. To change this, set \code{factorize=FALSE}. So, if your covariates
#' include multiple character columns, but you want to convert \emph{Scanner} to
#' binary instead of a factor, you may still specify \code{binarize='Scanner'}
#' and get the expected result. \code{binarize} will convert the given factor
#' variable(s) into numeric variable(s), which is performed \emph{before}
#' centering (if applicable).
#'
#' @section Centering:
#' The argument \code{mean.center} will mean-center (i.e., subtract the mean of
#' from each variable) any non-factor variables (including any dummy/indicator
#' covariates). This is done \emph{after} \dQuote{factorizing} and
#' \dQuote{binarizing}. If \code{center.how='all'}, then the \dQuote{grand mean}
#' will be used; otherwise, the groupwise means will be used. The grouping
#' variable is determined by \code{center.by} and is by default \code{'Group'}.
#'
#' @section Interactions:
#' \code{int} specifies which variables should interact with one another. This
#' argument accepts both numeric/continuous (e.g., \emph{Age}) and factor
#' variables (e.g., \emph{Sex}). All interaction combinations will be generated:
#' if you supply 3 variables, all two-way and the single three-way interaction
#' will be generated. This variable \emph{must} have at least two elements; it
#' is otherwise ignored. It is generally recommended that centering be performed
#' when including interaction terms.
#'
#' @param covars A \code{data.table} of covariates
#' @param coding Character string indicating how factor variables will be coded.
#'   Default: \code{'dummy'}
#' @param factorize Logical indicating whether to convert \emph{character}
#'   columns into \emph{factor}. Default: \code{TRUE}
#' @param binarize Character vector specifying the column name(s) of the
#'   covariate(s) to be converted from type \code{factor} to \code{numeric}.
#'   Default: \code{NULL}
#' @param int Character vector specifying the column name(s) of the
#'   covariate(s) to test for an interaction. Default: \code{NULL}
#' @param mean.center Logical indicating whether to mean center non-factor
#'   variables. Default: \code{FALSE}
#' @param center.how Character string indicating whether to use the grand mean
#'   or groupwise means. Default: \code{'all'}
#' @param center.by Character string indicating which grouping variable to use
#'   for calculating means (if applicable). Default: \code{'Group'}
#' @export
#'
#' @return A numeric matrix
#'
#' @name GLM design
#' @rdname glm_design
#' @family GLM functions
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

brainGraph_GLM_design <- function(covars, coding=c('dummy', 'effects', 'cell.means'),
                                  factorize=TRUE, binarize=NULL, int=NULL,
                                  mean.center=FALSE, center.how=c('all', 'within-groups'),
                                  center.by=getOption('bg.group')) {
  sID <- getOption('bg.subject_id')
  covars <- copy(covars)
  if (!sID %in% names(covars)) covars[, eval(sID) := as.character(seq_len(dim(covars)[1L]))]
  covars[, eval(sID) := as.character(get(sID))]
  X <- matrix(1, nrow=dim(covars)[1L], ncol=1)
  dimnames(X) <- list(covars[, get(sID)], 'Intercept')

  attrs <- list(factorize=factorize, mean.center=mean.center)
  if (isTRUE(factorize)) {
    cols <- names(which(vapply(covars, is.character, logical(1))))
    cols <- cols[!is.element(cols, sID)]
    cols <- setdiff(cols, binarize)
    if (length(cols) > 0) covars[, (cols) := lapply(.SD, as.factor), .SDcols=cols]
  }
  if (!is.null(binarize)) {
    stopifnot(all(binarize %in% names(covars)))
    covars[, (binarize) := lapply(.SD, function(x) as.numeric(x) - 1), .SDcols=binarize]
    attrs <- c(attrs, list(binarize=binarize))
  }

  nums <- which(vapply(covars, is.numeric, logical(1)))
  if (isTRUE(mean.center)) {
    center.how <- match.arg(center.how)
    covars[, (nums) := lapply(.SD, as.numeric), .SDcols=nums]
    if (center.how == 'all') {
      covars[, (nums) := lapply(.SD, function(x) x - mean(x)), .SDcols=nums]
    } else if (center.how == 'within-groups') {
      covars[, (nums) := lapply(.SD, function(x) x - mean(x)), .SDcols=nums, by=center.by]
    }
    attrs <- c(attrs, list(center.how=center.how, center.by=center.by))
  }
  if (length(nums) > 0) X <- cbind(X, as.matrix(covars[, ..nums, with=FALSE]))

  factors <- names(which(vapply(covars, class, character(1)) == 'factor'))
  coding <- match.arg(coding)
  attrs <- c(list(coding=coding), attrs)
  base_val <- switch(coding, effects=-1L, dummy=,cell.means=0L)
  starting <- switch(coding, cell.means=1, dummy=,effects=2)
  colRemove <- switch(coding, cell.means=999L, dummy=,effects=1L)  # Hack to keep all names w/ cell means
  if (coding == 'cell.means' && all(X[, 1L] == 1)) X <- X[, -1L, drop=FALSE]  # Remove intercept term
  for (f in factors) X <- create_dummy_vars(covars, X, f, base_val, starting, colRemove)

  if (!is.null(int) && length(int) > 1) {
    stopifnot(all(int %in% names(covars)))
    intcomb <- combn(int, 2, simplify=FALSE)
    if (length(int) == 3) intcomb <- c(intcomb, combn(int, 3, simplify=FALSE))
    for (x in intcomb) X <- get_int(X, coding, factors, x)
    attrs <- c(attrs, list(int=int))
  }

  attributes(X) <- c(attributes(X), attrs)
  return(X)
}

#' Create dummy variables from a factor
#'
#' Create dummy variables from a factor variable in a data.table.
#'
#' @param fact Character string specifying the column name of the factor
#' @param base_val Integer value for the \dQuote{base} group or factor level.
#'   For dummy and cell-means coding, this is 0. For effects coding, it is -1.
#' @param starting Numeric; the starting value to iterate over factor levels.
#'   For dummy and effects coding, this is 2. For cell-means, it is 1.
#' @param colRemove Numeric specifying which factor level to exclude when
#'   assigning names. For dummy and effects coding, this is 1 (i.e., there is
#'   no column for the base factor level in the matrix). For cell-means, this is
#'   999 (a hack to make it work without knowing the number of levels in
#'   advance)
#' @return A numeric matrix with new columns
#' @keywords internal

create_dummy_vars <- function(covars, X, fact, base_val, starting, colRemove) {
  cov.levels <- covars[, levels(get(fact))]
  cov.vec <- covars[, as.numeric(get(fact))]
  for (i in starting:max(cov.vec)) {
    cov.vec.sub <- ifelse(cov.vec == i, 1L, ifelse(cov.vec == 1, base_val, 0))
    X <- cbind(X, cov.vec.sub)
  }
  p <- dim(X)[2L]
  colnames(X)[(p - max(cov.vec) + starting):p] <- paste0(fact, cov.levels[-colRemove])
  return(X)
}

get_int <- function(X, coding, factors, int) {
  get_colnames <- function(string, X) {
    cnames <- grep(string, colnames(X), value=TRUE)
    if (any(grepl(':', cnames))) cnames <- cnames[-grep(':', cnames)]
    return(cnames)
  }

  p <- dim(X)[2L]
  intnames <- lapply(int, get_colnames, X)

  if (length(int) == 3) {
    X <- cbind(X, X[, intnames[[1]]] * X[, intnames[[2]]] * X[, intnames[[3]]])
    colnames(X)[(p + 1):dim(X)[2L]] <- paste(intnames[[1]], intnames[[2]], intnames[[3]], sep=':')
  } else {
    if (all(int[1:2] %in% factors) && coding == 'cell.means') {
      for (j in seq_along(intnames[[1]])) {
        X <- cbind(X, X[, intnames[[1]][j]] * X[, intnames[[2]]])
        colnames(X)[(p + 1):dim(X)[2L]] <- paste(intnames[[1]][j], intnames[[2]], sep=':')
        p <- dim(X)[2L]
      }
      X <- X[, -which(colnames(X) %in% intnames[[1]])]
      X <- X[, -which(colnames(X) %in% intnames[[2]])]
    } else {  # One of the int. terms is numeric
      X <- cbind(X, X[, intnames[[1]]] * X[, intnames[[2]]])
      colnames(X)[(p + 1):dim(X)[2L]] <- paste(intnames[[1]], intnames[[2]], sep=':')
      if (coding == 'cell.means') {
        if (!int[1] %in% factors) {
          X <- X[, -which(colnames(X) %in% intnames[[1]])]
        } else if (!int[2] %in% names(factors)) {
          X <- X[, -which(colnames(X) %in% intnames[[2]])]
        }
      }
    }
  }
  return(X)
}
