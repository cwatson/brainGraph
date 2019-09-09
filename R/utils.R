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
  x <- if ('degree' %in% vertex_attr_names(g)) V(g)$degree else degree(g)
  return(x)
}

#' Check for presence of a strength attribute
#'
#' Helper function to check if \code{strength} is a vertex attribute of the
#' input graph. Returns a numeric vector of the strength values.
#' @keywords internal

check_strength <- function(g) {
  x <- if ('strength' %in% vertex_attr_names(g)) V(g)$strength else strength(g)
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
  p <- switch(alt,
              less=pnorm(diff.z),
              greater=pnorm(diff.z, lower.tail=FALSE),
              two.sided=2 * pnorm(abs(diff.z), lower.tail=FALSE))

  return(list(p=p, z=diff.z))
}

#' Add metadata to a list-like object
#'
#' Adds metadata to a list-like object.
#'
#' If the object is a graph, graph-level attributes will be added. The
#' elements added are:
#' \itemize{
#'   \item{version}{A list with R, brainGraph, and igraph versions}
#'   \item{sys}{Character vector of system information}
#'   \item{date}{The date and time of creation}
#' }
#'
#' @param object A list-like object
#' @keywords internal
#' @return The same object with version, system, and date information added

get_metadata <- function(object) {
  object$version <- list(r=R.version.string,
                         bg=packageVersion('brainGraph'),
                         ig=packageVersion('igraph'))
  object$sys <- Sys.info()[c(1:3, 5)]
  object$date <- format(as.POSIXct(Sys.time()), '%Y-%m-%dT%H:%M:%OS0')
  return(object)
}

#' Calculate the threshold values to result in a specific density
#'
#' Given a vector of densities, this function will return the numeric values
#' that will result in graphs of the given densities. It is assumed that the
#' density values are all between 0 and 1.
#'
#' @param mat Numeric matrix that will be thresholded
#' @param densities Numeric vector of densities
#' @param emax Integer; the maximum number of edges
#' @param ... Arguments passed to \code{\link{sort}}
#' @return Numeric vector of thresholds
#' @keywords internal

get_thresholds <- function(mat, densities, emax=dim(mat)[1L] * (dim(mat)[1L] - 1) / 2, ...) {
  sort(mat[lower.tri(mat)], ...)[emax - densities * emax]
}

#' Create a data.table with a single graph metric
#'
#' Used in \code{brainGraph_GLM} and \code{brainGraph_mediate} to create a
#' \code{data.table} with the \code{Study.ID} and column(s) for the graph- or
#' vertex-level metric of interest.
#'
#' @inheritParams GLM
#' @return A \code{data.table} with one column containing the \code{Study.ID}
#' and 1 or more columns with the graph- or vertex-level measure of interest.
#' @keywords internal

glm_data_table <- function(g.list, level, measure) {
  if (level == 'vertex') {
    y <- t(vapply(g.list, vertex_attr, numeric(vcount(g.list[[1]])), measure))
    colnames(y) <- V(g.list[[1]])$name
    DT.y <- as.data.table(y, keep.rownames='Study.ID')
  } else if (level == 'graph') {
    DT.y <- as.data.table(vapply(g.list, graph_attr, numeric(1), measure), keep.rownames=TRUE)
    setnames(DT.y, names(DT.y), c('Study.ID', 'graph'))
  }
  return(DT.y)
}

#' Check if a matrix is binary
#'
#' @return Logical of length 1
#' @keywords internal

is_binary <- function(mat) {
  x <- identical(sum(abs(mat)) - sum(mat == 1), 0)
  return(x)
}

#' Convert a matrix to a list of rows
#'
#' Makes working with different contrast types (i.e., t or F) a little simpler.
#' @param mat Numeric matrix in which each row is a single contrast vector
#' @return A list with length equal to the number of rows of \code{C}
#' @keywords internal

matrix2list <- function(mat) {
  l <- lapply(seq_len(dim(mat)[1L]), function(i) mat[i, , drop=FALSE])
  nam <- dimnames(mat)[[1]]
  if (!is.null(nam)) names(l) <- nam
  return(l)
}

#' Helper function to calculate a max or min
#'
#' @keywords internal
maxfun <- function(alternative) {
  fun <- switch(alternative,
                two.sided=function(x) max(abs(x), na.rm=TRUE),
                less=function(x) min(x, na.rm=TRUE),
                greater=function(x) max(x, na.rm=TRUE))
  return(fun)
}

#' Helper function to sort values
#'
#' @keywords internal
sortfun <- function(alternative) {
  fun <- switch(alternative,
                two.sided=function(x) sort(abs(x)),
                less=function(x) sort(x, decreasing=TRUE),
                greater=sort)
  return(fun)
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
#' rotate points 90 deg. if given \dQuote{pi/2} as input)
#'
#' @param x A matrix with 2 columns of the points to rotate
#' @param theta The angle to apply
#'
#' @keywords internal
#' @return A matrix with 2 columns of the points' new locations
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

rotation <- function(x, theta) {
  R <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)),
              nrow=2, ncol=2, byrow=FALSE)
  x.rot <- x %*% R
  return(x.rot)
}

#' Subset graphs based on a given logical condition
#'
#' \code{subset_graph} will subset a given graph based on the given logical
#' condition(s). This can be a \dQuote{simple} logical equation, or can include
#' combinations of \emph{AND} and \emph{OR} (i.e., \code{&} and \code{|}).
#'
#' @param g A graph object
#' @param subgraph Character string specifying an equation for which vertices to
#'   keep
#' @keywords internal
#' @return A graph object

subset_graph <- function(g, subgraph) {
  stopifnot(nzchar(subgraph))
  orig.class <- class(g)
  # Function for creating the condition string to later subset the graph
  get_cond_string <- function(orig) {
    spec <- '\\s\\&\\s|\\s\\|\\s'
    substrings <- strsplit(orig, split=spec)[[1]]
    if (length(substrings) > 1) {  # Multiple conditions
      if (!isTRUE(grepl(spec, orig))) {
        stop('Logical operators must be surrounded by spaces!')
      }
      nchars <- cumsum(nchar(substrings))
      endpts <- nchars + vapply(seq_along(nchars), function(x) (3L * x) - 1L, integer(1))
      splits <- vapply(endpts, function(x) substr(orig, start=x, stop=x), character(1))
      substrings <- trimws(substrings) # Remove unnecessary whitespace

      cond.string <- paste(vapply(seq_along(substrings), function(x)
                                  paste0('V(g)$', substrings[x], splits[x]), character(1)),
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
    cond.strings <- vapply(subs, get_cond_string, character(1))
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
