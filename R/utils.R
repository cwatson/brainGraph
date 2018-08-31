#' Convert arguments into a single list
#'
#' \code{args_as_list} converts its arguments into a single list. If the
#' argument is already a single list, it returns the same object. If multiple
#' objects are given, it wraps them in a list.
#'
#' @param ... Graphs or a list of graphs
#' @keywords internal
#' @return A list object

args_as_list <- function(...) {
  graphs <- unlist(recursive=FALSE,
    lapply(list(...), function(l) {
             if (is_igraph(l)) {
               list(l)
             } else {
               l
             }}
    )
  )
  return(graphs)
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
    if (is.null(dim(y))) {  # A single vector, for MTPC
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

#' Check for edge weights
#'
#' \code{check_weights} is a helper function for dealing with edge weights that
#' get passed to different \code{igraph} functions.
#'
#' If \code{weights=NULL} and the graph has a \code{'weight'} attribute, that
#' will be returned. If \code{weights=NA}, then that is returned.
#' @keywords internal

check_weights <- function(g, weights) {
  if (is.null(weights) && 'weight' %in% edge_attr_names(g)) {
    weights <- NULL
  } else {
    if (!is.null(weights) && any(!is.na(weights))) {
      weights <- as.numeric(weights)
    } else {
      weights <- NA
    }
  }
  return(weights)
}

#' Check for presence of a degree attribute
#'
#' Helper function to check if \code{degree} is a vertex attribute of the input
#' graph. Returns a numeric vector of the degree values.
#' @keywords internal

check_degree <- function(g) {
  if ('degree' %in% vertex_attr_names(g)) {
    x <- V(g)$degree
  } else {
    x <- degree(g)
  }
  return(x)
}

#' Check for presence of a strength attribute
#'
#' Helper function to check if \code{strength} is a vertex attribute of the
#' input graph. Returns a numeric vector of the strength values.
#' @keywords internal

check_strength <- function(g) {
  if ('strength' %in% vertex_attr_names(g)) {
    x <- V(g)$strength
  } else {
    x <- strength(g)
  }
  return(x)
}

#' Calculate coefficient of variation
#'
#' Calculates the \emph{coefficient of variation}, defined as
#' \deqn{CV(x) = \frac{sd(x)}{mean(x)}}
#'
#' @param x Numeric vector
#' @export
#'
#' @return A numeric value

coeff_var <- function(x) {
  N <- length(x)
  mu <- sum(x) / N
  return(sqrt(1 / (N - 1) * (sum((x - mu)^2))) / mu)
}

#' Calculate the p-value for differences in correlation coefficients
#'
#' Given two sets of correlation coefficients and sample sizes, this function
#' calculates and returns the \emph{z-scores} and \emph{p-values} associated
#' with the difference between correlation coefficients. This function was
#' adapted from \url{http://stackoverflow.com/a/14519007/3357706}.
#'
#' @param r1 Numeric (vector or matrix) of correlation coefficients, group 1
#' @param r2 Numeric (vector or matrix) of correlation coefficients, group 2
#' @param n1 Integer; number of observations, group 1
#' @param n2 Integer; number of observations, group 2
#' @param alternative Character string specifying the alternative hypothesis
#'   test to use; one of: 'two.sided' (default), 'less', 'greater'
#' @export
#'
#' @return A list containing:
#' \item{p}{The p-values}
#' \item{z}{The z-score for the difference in correlation coefficients}
#'
#' @family Matrix functions
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' kNumSubjs <- summary(covars$Group)
#' corr.diffs <- cor.diff.test(corrs[[1]][[1]]$R, corrs[[2]][[1]]$R,
#'                             kNumSubjs[1], kNumSubjs[2], alternative='two.sided')
#' edge.diffs <- t(sapply(which(corr.diffs$p < .05), function(x)
#'                        mapply('[[',
#'                               dimnames(corr.diffs$p),
#'                               arrayInd(x, dim(corr.diffs$p)))
#'                               ))
#' }
cor.diff.test <- function(r1, r2, n1, n2,
                          alternative = c('two.sided', 'less', 'greater')) {

  z1 <- 0.5 * log((1 + r1) / (1 - r1))
  z2 <- 0.5 * log((1 + r2) / (1 - r2))

  SEdiff <- sqrt((1 / (n1 - 3)) + (1 / (n2 - 3)))
  diff.z <- (z1 - z2) / SEdiff

  alt <- match.arg(alternative)
  if (alt == 'less') {
    p <- pnorm(diff.z)
  } else if (alt == 'greater') {
    p <- pnorm(diff.z, lower.tail=F)
  } else if (alt == 'two.sided') {
    p <- 2 * pnorm(abs(diff.z), lower.tail=F)
  }

  return(list(p=p, z=diff.z))
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
#' @keywords internal
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

#' Symmetrize a matrix with the mean of off-diagonal elements
#'
#' \code{symm_mean} returns a symmetric matrix in which the off-diagonal
#' elements \eqn{A[i, j]} and \eqn{A[j, i]} are equal to the mean of the values
#' in the input matrix.
#' @param A Numeric matrix
#' @keywords internal
#' @return Numeric matrix

symm_mean <- function(A) {
  0.5 * (A + t(A))
}

#' Apply a rotation matrix to a set of points
#'
#' This function takes a set of points and applies a rotation matrix (e.g. will
#' rotate points 90 deg. if given "pi/2" as input)
#'
#' @param x A matrix with 2 columns of the points to rotate
#' @param theta The angle to apply
#'
#' @keywords internal
#' @return A matrix with 2 columns of the points' new locations
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

rotation <- function(x, theta) {
  R <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)),
              nrow=2, ncol=2, byrow=F)
  x.rot <- x %*% R
  return(x.rot)
}

#' Color graph vertices and edges
#'
#' \code{set_vertex_color} takes an integer vector representing membership of
#' some grouping (e.g., a community or connected component) and creates a
#' character vector of colors for each grouping. Isolated vertices will be
#' colored \emph{gray}.
#'
#' @param g An \code{igraph} graph object
#' @param name Character string of the name of the vertex attribute to add
#' @param memb An integer vector representing membership of e.g. a community
#'
#' @return The same graph with additional vertex or edge attribute
#'
#' @keywords internal
#' @name GraphColors
#' @aliases set_vertex_color
#' @rdname color_vertices_edges

set_vertex_color <- function(g, name, memb) {
  big.groups <- which(as.integer(table(memb)) > 1)

  group.cols.memb <- rep('gray', length=max(memb))
  group.cols.memb[big.groups] <- group.cols[big.groups]

  g <- set_vertex_attr(g, name, value=group.cols.memb[memb])
  return(g)
}

#' Color graph edges
#'
#' \code{set_edge_color} assigns a color to each edge (the same as the vertex
#' membership colors). Edges that connect vertices of two different groups are
#' colored gray.
#'
#' @keywords internal
#' @aliases set_edge_color
#' @rdname color_vertices_edges

set_edge_color <- function(g, name, memb) {
  stopifnot(length(memb) == vcount(g))
  big.groups <- as.integer(names(which(table(memb) > 1)))

  newcols <- rep('gray50', length=ecount(g))
  tmp <- vector('list', length=max(big.groups))
  for (i in big.groups) {
    x <- which(memb == i)
    tmp[[i]] <- as.vector(E(g)[x %--% x])
    if (!is.null(tmp[[i]])) newcols[tmp[[i]]] <- group.cols[i]
  }

  g <- set_edge_attr(g, name, value=newcols)
  return(g)
}

#' Subset graphs based on a given logical condition
#'
#' \code{subset_graph} will subset a given graph based on the given logical
#' condition(s). This can be a "simple" logical equation, or can include
#' combinations of \emph{AND} and \emph{OR} (i.e., \code{&} and \code{|}).
#'
#' @param g A graph object
#' @param subgraph Character string specifying an equation for which vertices to
#'   keep
#' @keywords internal
#' @return A graph object

subset_graph <- function(g, subgraph) {
  stopifnot(nchar(subgraph) > 0)
  orig.class <- class(g)
  # Function for creating the condition string to later subset the graph
  get_cond_string <- function(orig) {
    substrings <- strsplit(orig, split='\\s\\&\\s|\\s\\|\\s')[[1]]
    if (length(substrings) > 1) {  # Multiple conditions
      if (!isTRUE(grepl('\\s\\&\\s|\\s\\|\\s', orig))) {
        stop('Logical operators must be surrounded by spaces!')
      }
      nchars <- cumsum(sapply(substrings, nchar))
      splits <- sapply(seq_along(substrings), function(x)
                       substr(orig, start=nchars[x]+(3*x-1), stop=nchars[x]+(3*x-1)))
      substrings <- gsub('^\\s+|\\s+$', '', substrings) # Remove unnecessary whitespace

      cond.string <- paste(sapply(seq_along(substrings), function(x)
                                  paste0('V(g)$', substrings[x], splits[x])),
                           collapse='')
    } else {
      cond.string <- paste0('V(g)$', substrings)
    }
    return(cond.string)
  }

  # Handle when logical expressions are separated by parentheses
  if (isTRUE(grepl('\\(.*\\&.*\\)', subgraph)) || isTRUE(grepl('\\(.*\\|.*\\)', subgraph))) {
    subs <- strsplit(subgraph, split='\\)\\s\\&\\s\\(')[[1]]
    subs <- as.list(gsub('^\\(|\\)$|^\\s+|\\s+$', '', subs))
    cond.strings <- sapply(subs, get_cond_string)
    cond.string <- paste0('(', cond.strings[1], ') & (', cond.strings[2], ')')
  } else {
    cond.string <- get_cond_string(subgraph)
  }

  cond <- eval(parse(text=cond.string))
  if (sum(cond, na.rm=TRUE) == 0) {
    warning('No vertices meet criteria! No graph created')
    return(list(g=NULL, inds=NULL))
  } else {
    inds <- which(cond)
    cond <- setdiff(seq_len(vcount(g)), which(cond))
    g <- delete.vertices(g, cond)
    class(g) <- orig.class
  }
  return(list(g=g, inds=inds))
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
#' @keywords internal
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
#' There are 3 options for the type of transform to apply:
#' \enumerate{
#'   \item \code{1/w}: calculate the inverse
#'   \item \code{-log(w)}: calculate the negative (natural) logarithm
#'   \item \code{1-w}: subtract each weight from 1
#' }
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
#' @return An \code{igraph} graph object with transformed edge weights and a
#'   graph attribute, \code{xfm.type}, of the type of transform
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

xfm.weights <- function(g, xfm.type=c('1/w', '-log(w)', '1-w'), invert=FALSE) {
  stopifnot(is_igraph(g), is_weighted(g))
  xfm.type <- match.arg(xfm.type)
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

  g$xfm.type <- xfm.type
  return(g)
}
