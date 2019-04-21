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
#' @name GLMdesign
#' @aliases brainGraph_GLM_design
#' @rdname glm_design
#' @family GLM functions
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

brainGraph_GLM_design <- function(covars, coding=c('dummy', 'effects', 'cell.means'),
                                  factorize=TRUE, binarize=NULL, int=NULL,
                                  mean.center=FALSE, center.how=c('all', 'within-groups'),
                                  center.by='Group') {
  Study.ID <- NULL
  covars <- copy(covars)
  covars[, Study.ID := as.character(Study.ID)]
  X <- matrix(1, nrow=nrow(covars), ncol=1)
  colnames(X) <- 'Intercept'

  if (isTRUE(factorize)) {
    cols <- names(which(sapply(covars, is.character)))
    cols <- cols[!is.element(cols, 'Study.ID')]
    cols <- setdiff(cols, binarize)
    for (z in cols) covars[, eval(z) := as.factor(get(z))]
  }
  if (!is.null(binarize)) {
    stopifnot(all(binarize %in% names(covars)))
    for (b in binarize) covars[, eval(b) := as.numeric(get(b)) - 1]
  }

  center.how <- match.arg(center.how)
  nums <- which(sapply(covars, is.numeric))
  if (isTRUE(mean.center)) {
    covars[, (nums) := lapply(.SD, as.numeric), .SDcols=nums]
    if (center.how == 'all') {
      covars[, (nums) := lapply(.SD, function(x) x - mean(x)), .SDcols=nums]
    } else if (center.how == 'within-groups') {
      covars[, (nums) := lapply(.SD, function(x) x - mean(x)), .SDcols=nums, by=center.by]
    }
    for (n in nums) covars[[n]] <- covars[[n]] - mean(covars[[n]])
  }
  if (length(nums) > 0) X <- cbind(X, as.matrix(covars[, nums, with=F]))

  factors <- which(sapply(covars, class) == 'factor')
  coding <- match.arg(coding)
  for (f in factors) {
    cov.name <- names(covars)[f]
    cov.levels <- covars[, levels(get(cov.name))]
    cov.vec <- covars[, as.numeric(get(cov.name))]

    if (coding == 'cell.means') {
      for (i in 1:max(cov.vec)) {
        cov.vec.sub <- ifelse(cov.vec == i, 1, 0)
        X <- cbind(X, cov.vec.sub)
      }
      if (all(X[, 1] == 1)) X <- X[, -1]  # Remove intercept term
      p <- ncol(X)
      colnames(X)[(p - max(cov.vec) + 1):p] <- paste0(cov.name, cov.levels)

    } else {
      code <- ifelse(coding == 'dummy', 0, -1)
      for (i in 2:max(cov.vec)) {
        cov.vec.sub <- ifelse(cov.vec == i, 1, 0)
        cov.vec.sub <- ifelse(cov.vec == 1, code, cov.vec.sub)
        X <- cbind(X, cov.vec.sub)
      }
      p <- ncol(X)
      colnames(X)[(p - (max(cov.vec) - 1) + 1):p] <- paste0(cov.name, cov.levels[-1])
    }
  }

  if (!is.null(int) & length(int) > 1) {
    stopifnot(all(int %in% names(covars)))
    intcomb <- combn(int, 2, simplify=FALSE)
    if (length(int) == 3) intcomb <- c(intcomb, combn(int, 3, simplify=FALSE))
    for (x in intcomb) X <- get_int(X, coding, names(factors), x)
  }

  return(X)
}

get_int <- function(X, coding, factors, int) {
  get_colnames <- function(string, X, factors) {
    cnames <- colnames(X)[grep(string, colnames(X))]
    if ((!string %in% factors) && length(cnames) > 1) {
      cnames <- colnames(X)[grep(paste0('^', string, '$'), colnames(X))]
    }
    if (any(grepl(':', cnames))) cnames <- cnames[-grep(':', cnames)]
    return(cnames)
  }

  p <- ncol(X)
  intnames <- lapply(int, get_colnames, X, factors)

  if (length(int) == 3) {
    X <- cbind(X, X[, intnames[[1]]] * X[, intnames[[2]]] * X[, intnames[[3]]])
    colnames(X)[(p + 1):ncol(X)] <- paste(intnames[[1]], intnames[[2]], intnames[[3]], sep=':')
  } else {
    if (int[1] %in% factors && int[2] %in% factors && coding == 'cell.means') {
      for (j in seq_along(intnames[[1]])) {
        X <- cbind(X, X[, intnames[[1]][j]] * X[, intnames[[2]]])
        p2 <- ncol(X)
        colnames(X)[(p + 1):p2] <- paste(intnames[[1]][j], intnames[[2]], sep=':')
        p <- ncol(X)
      }
      X <- X[, -which(colnames(X) %in% intnames[[1]])]
      X <- X[, -which(colnames(X) %in% intnames[[2]])]
    } else {  # One of the int. terms is numeric
      X <- cbind(X, X[, intnames[[1]]] * X[, intnames[[2]]])
      p2 <- ncol(X)
      colnames(X)[(p + 1):p2] <- paste(intnames[[1]], intnames[[2]], sep=':')
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
