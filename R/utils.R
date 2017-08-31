#' Difference in the area-under-the-curve of two vectors
#'
#' This function takes two vectors, calculates the area-under-the-curve (AUC),
#' and calculates the difference between the two.
#'
#' @param x Numeric vector of the x-values
#' @param y A 2-column numeric matrix; each column contains the values for one
#'   group
#'
#' @return A numeric value of the difference between two groups

auc_diff <- function(x, y) {
  if (length(x) > 1) {
    if (is.null(dim(y))) {  # A single vector, for MTPC
      return(sum(-diff(x) * (head(y, -1) + tail(y, -1))) / 2)
    } else {
      return(-diff(apply(y, 2, function(z)
                         sum(diff(x) * (head(z, -1) + tail(z, -1))) / 2)))
    }
  } else {
    return(y[1] - y[2])
  }
}

#' Calculate coefficient of variation
#'
#' Calculates the \emph{coefficient of variation}, defined as
#' \deqn{CV(x) = \frac{sd(x)}{mean(x)}}
#'
#' @param x Numeric vector
#'
#' @return A numeric value

coeff_var <- function(x) {
  N <- length(x)
  mu <- sum(x) / N
  return(sqrt(1 / (N - 1) * (sum((x - mu)^2))) / mu)
}

#' Delete all attributes of a graph
#'
#' Deletes all graph-, vertex-, and edge-level attributes of an \code{igraph}
#' graph object.
#'
#' @param g An \code{igraph} graph object
#' @param keep.names Logical indicating whether to keep the \code{name} vertex
#'   attribute (default: \code{FALSE})
#'
#' @return An \code{igraph} graph object
#' @seealso \code{\link[igraph]{delete_graph_attr},
#'   \link[igraph]{delete_vertex_attr}, \link[igraph]{delete_edge_attr}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

delete_all_attr <- function(g, keep.names=FALSE) {
  for (att in graph_attr_names(g)) g <- delete_graph_attr(g, att)
  for (att in edge_attr_names(g)) g <- delete_edge_attr(g, att)
  if (isTRUE(keep.names)) {
    vattrs <- setdiff(vertex_attr_names(g), 'name')
  } else {
    vattrs <- vertex_attr_names(g)
  }
  for (att in vattrs) g <- delete_vertex_attr(g, att)

  return(g)
}

#' Apply a rotation matrix to a set of points
#'
#' This function takes a set of points and applies a rotation matrix (e.g. will
#' rotate points 90 deg. if given "pi/2" as input)
#'
#' @param x A matrix with 2 columns of the points to rotate
#' @param theta The angle to apply
#'
#' @return A matrix with 2 columns of the points' new locations
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

rotation <- function(x, theta) {
  R <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)),
              nrow=2, ncol=2, byrow=F)
  x.rot <- x %*% R
  return(x.rot)
}

#' Transform a vector to have a different range
#'
#' This function takes a vector and transforms it to have a new range, given
#' the input, or the default values of [0, 1].
#'
#' @param x the vector to transform
#' @param min.val the minimum value of the new range
#' @param max.val the maximum value of the new range
#'
#' @return A vector of the transformed input.

vec.transform <- function(x, min.val=0, max.val=1) {
  if (diff(range(x, na.rm=TRUE)) == 0) {
    return(rep(max.val, length=length(x)))
  } else {
    return(((x - min(x, na.rm=TRUE)) * (max.val - min.val) / diff(range(x, na.rm=TRUE))) + min.val)
  }
}

#' Transform edge weights
#'
#' For distance-based measures, it is important to transform the edge weights so
#' that the \emph{strongest} connections are re-mapped to having the
#' \emph{lowest} weights. Then you may calculate e.g., the \emph{shortest path
#' length} which will include the strongest connections.
#'
#' To transform the weights back to original values, specify \code{invert=TRUE}.
#'
#' @param g An \code{igraph} graph object
#' @param xfm.type Character string specifying how to transform the weights
#'   (default: \code{1/w})
#' @param invert Logical indicating whether or not to invert the transformation
#'   (default: \code{FALSE})
#' @export
#'
#' @return An \code{igraph} graph object with transformed edge weights
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

xfm.weights <- function(g, xfm.type=c('1/w', '-log(w)', '1-w'), invert=FALSE) {
  stopifnot(is_igraph(g), is_weighted(g))
  if (xfm.type == '1/w') {
    E(g)$weight <- 1 / E(g)$weight
  } else if (xfm.type == '-log(w)') {
    if (isTRUE(invert)) {
      E(g)$weight <- exp(-E(g)$weight)
    } else {
      E(g)$weight <- -log(E(g)$weight)
    }
  } else if (xfm.type == '1-w') {
    E(g)$weight <- 1 - E(g)$weight
  }
#  } else {
#    if (xfm.type == '1/w') {
#      E(g)$weight <- 1 / E(g)$weight
#    } else if (xfm.type == '-log(w)') {
#      E(g)$weight <- -log(E(g)$weight)
#    } else if (xfm.type == '1-w') {
#      E(g)$weight <- 1 - E(g)$weight
#    }
#  }
  return(g)
}
