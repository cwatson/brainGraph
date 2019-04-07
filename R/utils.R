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

#' Add metadata to a list-like object
#'
#' \code{get_metadata} adds metadata to a list-like object.
#'
#' The information added are:
#' \itemize{
#'   \item{version}{R, brainGraph, and igraph versions}
#'   \item{sys}{System information}
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
              nrow=2, ncol=2, byrow=F)
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
    substrings <- strsplit(orig, split='\\s\\&\\s|\\s\\|\\s')[[1]]
    if (length(substrings) > 1) {  # Multiple conditions
      if (!isTRUE(grepl('\\s\\&\\s|\\s\\|\\s', orig))) {
        stop('Logical operators must be surrounded by spaces!')
      }
      nchars <- cumsum(sapply(substrings, nchar))
      splits <- sapply(seq_along(substrings), function(x)
                       substr(orig, start=nchars[x]+(3*x-1), stop=nchars[x]+(3*x-1)))
      substrings <- trimws(substrings) # Remove unnecessary whitespace

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
