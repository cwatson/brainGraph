#' Calculate the AUC across densities of given attributes
#'
#' Given a list of \code{brainGraphList} objects, this function will calculate
#' the area under the curve (AUC) across all thresholds/densities for each
#' subject or group.
#'
#' @param g.list A list of \code{brainGraphList} objects
#' @param g.attr A character vector of graph attribute name(s). Default:
#'   \code{NULL}
#' @param v.attr A character vector of vertex attribute name(s). Default:
#'   \code{NULL}
#' @param norm Logical indicating whether to normalize threshold values to
#'   between 0 and 1. Default: \code{FALSE}
#' @export
#' @return A \code{brainGraphList} object in which the graphs are for each
#'   subject
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach getDoParRegistered
#'
#' @examples
#' \dontrun{
#' g.auc <- make_auc_brainGraph(g.fa, g.attr='E.global.wt')
#' }

make_auc_brainGraph <- function(g.list, g.attr=NULL, v.attr=NULL, norm=FALSE) {
  threshold <- i <- NULL
  # Check if components are 'brainGraphList' objects
  matches <- vapply(g.list, inherits, logical(1), 'brainGraphList')
  if (any(!matches)) stop("Input must be a list of 'brainGraphList' objects.")

  # Get the meta variables first
  attrs <- c('atlas', 'type', 'level', 'modality', 'weighting')
  out <- setNames(vector('list', length(attrs)), attrs)
  for (a in attrs) out[[a]] <- g.list[[1L]][[a]]

  if (!is.null(g.list[[1L]]$threshold)) {
    thresholds <- vapply(g.list, with, numeric(1), threshold)
  } else {
    thresholds <- seq(from=0, to=1, length.out=length(g.list))
  }
  if (isTRUE(norm)) thresholds <- vec.transform(thresholds)
  subjects <- names(g.list[[1L]]$graphs)
  grps <- groups(g.list[[1L]])

  if (!getDoParRegistered()) {
    cl <- makeCluster(getOption('bg.ncpus'))
    registerDoParallel(cl)
  }
  g.auc <- foreach(i=seq_along(subjects)) %dopar% {
    g.subj <- lapply(g.list, `[`, i)
    g.tmp <- with(out,
        make_empty_brainGraph(atlas, type=type, level=level, modality=modality,
                              weighting=weighting, name=subjects[i], Group=grps[i]))
    if (!is.null(g.attr)) {
      for (k in g.attr) {
        y <- sapply(g.subj, graph_attr, k)
        g.tmp <- set_graph_attr(g.tmp, k, abs(auc_diff(thresholds, y)))
      }
    }
    if (!is.null(v.attr)) {
      for (k in v.attr) {
        y <- t(sapply(g.subj, vertex_attr, k))
        g.tmp <- set_vertex_attr(g.tmp, k, value=abs(auc_diff(thresholds, y)))
      }
    }
    g.tmp
  }
  out <- get_metadata(out)
  out$graphs <- g.auc
  names(out$graphs) <- subjects
  class(out) <- c('brainGraphList', 'list')
  return(out)
}

#' Difference in the area-under-the-curve of two vectors
#'
#' This function takes two vectors, calculates the area-under-the-curve (AUC),
#' and calculates the difference between the two (if applicable).
#'
#' There are 4 different behaviors for this function:
#' \enumerate{
#'   \item If \code{x} is a single numeric value, then \code{y} should be a
#'     vector of 2 values and the difference is returned. This generally should
#'     not occur and may be removed in the future.
#'   \item If \code{y} has 1 column (or is a vector), then the AUC of \code{y}
#'     is returned.
#'   \item If \code{y} has exactly 2 columns, then each column should contain
#'     the values of interest for each subject group, and the difference in
#'     AUC's for each group is returned.
#'   \item If \code{y} has multiple columns (e.g., equal to the number of
#'     vertices of a graph), it will calculate the AUC for each column.
#' }
#'
#' @param x Numeric vector of the x-values
#' @param y A numeric vector or matrix
#'
#' @keywords internal
#' @return A numeric value of the difference between two groups, or a numeric
#'   vector of the AUC across vertices

auc_diff <- function(x, y) {
  if (length(x) > 1L) {
    if (NCOL(y) == 1L) {  # A single vector, for MTPC
      return(sum(-diff(x) * (head(y, -1L) + tail(y, -1L))) / 2)
    } else if (ncol(y) > 2L) {
      return(apply(y, 2L, function(z) auc_diff(x, z)))
    } else {
      return(auc_diff_perm(x, y))
    }
  } else {
    return(y[1L] - y[2L])
  }
}

#' @param x Numeric vector
#' @param y Numeric matrix with 2 columns and number of rows equal to the length
#'   of \code{x}
#' @return A single numeric value of the between-group AUC difference
#' @noRd

auc_diff_perm <- function(x, y) {
  -diff(apply(y, 2L, function(z) sum(diff(x) * (head(z, -1L) + tail(z, -1L))) / 2))
}
