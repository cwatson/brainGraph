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
  for (a in attrs) out[[a]] <- g.list[[1]][[a]]

  kNumThresh <- length(g.list)
  if (!is.null(g.list[[1]]$threshold)) {
    thresholds <- vapply(g.list, with, numeric(1), threshold)
  } else {
    thresholds <- seq(from=0, to=1, length.out=kNumThresh)
  }
  if (isTRUE(norm)) thresholds <- vec.transform(thresholds)
  subjects <- names(g.list[[1]]$graphs)
  kNumSubjs <- length(subjects)
  grps <- groups(g.list[[1]])
  g.auc <- foreach(i=seq_along(subjects)) %dopar% {
    g.subj <- lapply(g.list, `[`, i)
    g.tmp <- with(out,
        make_empty_brainGraph(atlas, type=type, level=level, modality=modality,
                              weighting=weighting, name=subjects[i], Group=grps[i]))
    if (!is.null(g.attr)) {
      for (k in seq_along(g.attr)) {
        y <- sapply(g.subj, graph_attr, g.attr[k])
        g.tmp <- set_graph_attr(g.tmp, g.attr[k], abs(auc_diff(thresholds, y)))
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
#' and calculates the difference between the two.
#'
#' If \code{y} has 2 columns, then each column should be the values for each
#' subject group. If \code{y} has multiple columns (e.g., equal to the number of
#' vertices of a graph), it will calculate the AUC for each column.
#'
#' @param x Numeric vector of the x-values
#' @param y A numeric matrix
#'
#' @keywords internal
#' @return A numeric value of the difference between two groups, or a numeric
#'   vector of the AUC across vertices

auc_diff <- function(x, y) {
  if (length(x) > 1) {
    if (NCOL(y) == 1) {  # A single vector, for MTPC
      return(sum(-diff(x) * (head(y, -1) + tail(y, -1))) / 2)
    } else if (ncol(y) > 2) {
      return(apply(y, 2, function(z) auc_diff(x, z)))
    } else {
      return(-diff(apply(y, 2, function(z)
                         sum(diff(x) * (head(z, -1) + tail(z, -1))) / 2)))
    }
  } else {
    return(y[1] - y[2])
  }
}
