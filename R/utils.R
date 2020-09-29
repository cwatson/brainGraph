#' Check if all elements in a list are identical
#'
#' @param l A list
#' @noRd

all.identical <- function(l) all(mapply(identical, head(l, 1L), tail(l, -1L)))

#' Convert arguments into a single list
#'
#' \code{args_as_list} converts its arguments into a single list. If the
#' argument is already a single list, it returns the same object. If multiple
#' objects are given, it wraps them in a list.
#'
#' @param ... Graphs or a list of graphs
#' @noRd

args_as_list <- function(...) {
  ll <- lapply(list(...), function(l) if (is_igraph(l)) list(l) else l)
  unlist(ll, recursive=FALSE)
}

#' Check for vertex or edge attributes
#'
#' \code{check_weights} is a helper function for dealing with edge weights that
#' get passed to different \code{igraph} functions.
#'
#' @return \code{check_weights} - If \code{weights=NULL} and the graph has a
#'   \code{'weight'} attribute, then \code{NULL} will be returned. If
#'   \code{weights=NA}, then \code{NA} is returned.
#' @keywords internal
#' @rdname check_attributes

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
  weights
}

#' Check for presence of a degree attribute
#'
#' \code{check_degree} is a helper function to check if \code{degree} is a
#' vertex attribute of the input graph. Returns a numeric vector of the degree
#' values.
#' @keywords internal
#' @rdname check_attributes

check_degree <- function(g) {
  if ('degree' %in% vertex_attr_names(g)) V(g)$degree else degree(g)
}

#' Check for presence of a strength attribute
#'
#' \code{check_strength} is a helper function to check if \code{strength} is a
#' vertex attribute of the input graph. Returns a numeric vector of the strength
#' values.
#' @keywords internal
#' @rdname check_attributes

check_strength <- function(g) {
  if ('strength' %in% vertex_attr_names(g)) V(g)$strength else strength(g)
}

#' Test if an object is a character vector of numbers
#'
#' \code{check_sID} is a convenience function to test if a vector (typically the
#' \emph{subject ID} column in a \code{data.table}) is a character vector of
#' numbers, a factor vector of numbers, or a numeric vector. If so, it will
#' zero-pad the variable to have equal width.
#'
#' This function is meant to avoid issues that arise when sorting a vector of
#' numbers that have been converted to \code{character}. For example,
#' \code{\link{import_scn}} automatically reads in the first column (with
#' \emph{FreeSurfer} outputs this is the column of subject IDs) as a
#' \code{character} variable. If the subject IDs had been all numbers/integers,
#' then sorting (i.e., setting the \code{key} in a \code{data.table}) would be
#' incorrect: e.g., it might be \code{'1', '10', '2', ...}.
#'
#' @return \code{check_sID} returns either the input vector or a character
#'   vector padded with \code{0}
#' @export
#' @rdname pad_zeros

check_sID <- function(x) {
  cls <- class(x)
  if (cls == 'factor') {
    test <- suppressWarnings(as.numeric(as.character(x)) == x)
  } else {
    test <- suppressWarnings(as.character(as.numeric(x)) == x)
  }
  if (isTRUE(all(test))) x <- pad_zeros(x)
  return(x)
}

#' Calculate the p-value for differences in correlation coefficients
#'
#' Given two sets of correlation coefficients and sample sizes, this function
#' calculates and returns the \emph{z-scores} and \emph{p-values} associated
#' with the difference between correlation coefficients. This function was
#' adapted from \url{https://stackoverflow.com/a/14519007/3357706}.
#'
#' @param r1,r2 Numeric (vector or matrix) of correlation coefficients for both
#'   groups
#' @param n Integer vector; number of observations for both groups
#' @inheritParams GLM
#' @export
#'
#' @return A list with elements \code{p} and \code{z}, the p-values
#'   and z-scores for the difference in correlations.
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' kNumSubjs <- summary(covars$Group)
#' corr.diffs <- cor.diff.test(corrs$R[, , 1], corrs$R[, , 2], kNumSubjs)
#' edge.diffs <- t(sapply(which(corr.diffs$p < .05), function(x)
#'                        mapply('[[',
#'                               dimnames(corr.diffs$p),
#'                               arrayInd(x, dim(corr.diffs$p)))
#'                               ))
#' }
cor.diff.test <- function(r1, r2, n, alternative=c('two.sided', 'less', 'greater')) {
  stopifnot(length(n) == 2L)

  z1 <- 0.5 * log((1 + r1) / (1 - r1))
  z2 <- 0.5 * log((1 + r2) / (1 - r2))

  n <- n - 3L
  SEdiff <- sqrt((1 / n[1L]) + (1 / n[2L]))
  diff.z <- (z1 - z2) / SEdiff

  alt <- match.arg(alternative)
  p <- switch(alt,
              less=pnorm(diff.z),
              greater=pnorm(diff.z, lower.tail=FALSE),
              two.sided=2 * pnorm(abs(diff.z), lower.tail=FALSE))

  return(list(p=p, z=diff.z))
}

#' Utility functions
#'
#' \code{get_metadata} adds metadata to a list-like object.
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
#' @return \code{get_metadata} - the same object with version, system, and date
#'   information added
#' @name Utility functions
#' @rdname utils

get_metadata <- function(object) {
  object$version <- list(r=R.version.string,
                         bg=packageVersion('brainGraph'),
                         ig=packageVersion('igraph'))
  object$sys <- Sys.info()[c(1L:3L, 5L)]
  object$date <- format(Sys.time(), '%Y-%m-%dT%H:%M:%OS0')
  return(object)
}

#' Create a data.table with a single graph metric
#'
#' \code{glm_data_table} is used in \code{brainGraph_GLM} and
#' \code{brainGraph_mediate} to create a \code{data.table} with the
#' \emph{subject IDs} and column(s) for the graph- or vertex-level metric of
#' interest.
#'
#' @inheritParams GLM
#' @return \code{glm_data_table} - A \code{data.table} with one column
#'   containing the subject ID's and 1 or more columns with the graph- or
#'   vertex-level measure of interest.
#' @keywords internal
#' @rdname glm_helpers

glm_data_table <- function(g.list, level, measure) {
  sID <- getOption('bg.subject_id')
  if (level == 'vertex') {
    y <- t(vapply(g.list, vertex_attr, numeric(vcount(g.list[[1L]])), measure))
    colnames(y) <- V(g.list[[1L]])$name
    DT.y <- as.data.table(y, keep.rownames=sID)
  } else if (level == 'graph') {
    DT.y <- as.data.table(vapply(g.list, graph_attr, numeric(1L), measure), keep.rownames=TRUE)
    setnames(DT.y, names(DT.y), c(sID, 'graph'))
  }
  return(DT.y)
}

#' Convert a matrix to a list of rows
#'
#' \code{matrix2list} makes working with different contrast types (i.e., t or F)
#' a little simpler.
#'
#' @param mat Numeric matrix in which each row is a single contrast vector
#' @return \code{matrix2list} -- A list with length equal to the number of rows
#'   of \code{C}
#' @keywords internal
#' @rdname glm_helpers

matrix2list <- function(mat) {
  l <- lapply(seq_len(dim(mat)[1L]), function(i) mat[i, , drop=FALSE])
  nam <- dimnames(mat)[[1L]]
  if (!is.null(nam)) names(l) <- nam
  return(l)
}

#' Helper function to calculate a max or min
#'
#' @keywords internal
#' @rdname glm_helpers
maxfun <- function(alternative) {
  switch(alternative,
         two.sided=function(x) max(abs(x), na.rm=TRUE),
         less=function(x) min(x, na.rm=TRUE),
         greater=function(x) max(x, na.rm=TRUE))
}

#' Helper function to sort values
#'
#' @keywords internal
#' @rdname glm_helpers
sortfun <- function(alternative) {
  switch(alternative,
         two.sided=function(x, ind) sort(abs(x), partial=ind)[ind],
         less=function(x, ind) sort(x, decreasing=TRUE)[ind],
         greater=function(x, ind) sort(x, partial=ind)[ind])
}

#' Faster version of outer for 2 vectors
#'
#' \code{outer_vec} simply performs the cross-product, specifically \code{x %*%
#' t(y)}, and assigns dimnames to the resulting matrix.
#' @keywords internal
#' @rdname utils

outer_vec <- function(x, y) {
  robj <- tcrossprod(x, y)
  dimnames(robj) <- list(names(x), names(y))
  robj
}

#' Pad the front of a numeric or character vector with zeros
#'
#' \code{pad_zeros} pads a vector with zeros to avoid issues with ordering a
#' column of integers or integers converted to \code{character}.
#'
#' If \dQuote{x} is a numeric vector, then the resultant string width will be
#' determined by \code{max(x)} or \code{x} itself if the input is a single
#' integer. For example, if \code{x=10}, it will return \code{'01', '02', ...,
#' '10'}. If \dQuote{x} is a character vector, then the output's string width
#' will be \code{max(nchar(x))}. For example, if \code{x} includes both
#' \code{'1'} and \code{'1000'}, it will return \code{'0001'}, etc.
#'
#' @param x \code{pad_zeros} accepts either a vector (numeric or character) or a
#'   single integer. \code{check_sID} accepts a character, numeric, or factor
#'   vector
#' @return A character vector with zero-padded values
#' @export
#' @examples
#' pad_zeros(10)  # '01' '02' ... '10'
#' x <- c(1, 10, 100)
#' pad_zeros(x)   # '001' '010' '100'
#' x <- as.character(x)
#' pad_zeros(x)   # '001' '010' '100'

pad_zeros <- function(x) {
  if (is.numeric(x)) {
    if (length(x) == 1L) x <- seq_len(x)
    n <- max(x)
    x <- formatC(x, width=floor(log10(n) + 1L), flag='0')
  } else if (is.character(x)) {
    nc <- nchar(x)
    if (length(unique(nc)) == 1L) return(x)
    spec <- paste0('%0', max(nc), 's')
    x <- gsub(' ', '0', sprintf(spec, x))
  }
  return(x)
}

#' Apply a rotation matrix to a set of points
#'
#' This function takes a set of points and applies a rotation matrix (e.g. will
#' rotate points 90 deg. if given \dQuote{pi/2} as input)
#'
#' @param x A matrix with 2 columns of the points to rotate
#' @param theta The angle to apply
#' @return A matrix with 2 columns of the points' new locations
#' @noRd

rotation <- function(x, theta) {
  R <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow=2L)
  x %*% R
}

#' Capitalize the first letter of a character string
#'
#' @noRd
simpleCap <- function(x) paste0(toupper(substring(x, 1L, 1L)), substring(x, 2L))

#' Add newlines to a character string for printing
#'
#' \code{split_string} inserts separator characters in a character string to
#' truncate strings for printing.
#'
#' The \code{delim} argument determines \emph{where} to insert the separator.
#' For example, given the default, it will attempt to insert newline characters
#' only following a plus sign.
#'
#' @param x A character string
#' @param max.len Integer; the max length of one line. Default: 80
#' @param delim Character specifying where to end a line if it is longer than
#'   \code{max.len}. Default: \code{'+'}
#' @param sep Character specifying what to split by. Default: \code{'\n'} (a
#'   newline character)
#' @keywords internal
#' @rdname utils

split_string <- function(x, max_len=80L, delim='\\+', sep='\n') {
  str_len <- nchar(x)
  if (str_len > max_len) {
    nlines <- (str_len %/% max_len) + (str_len %% max_len > 0L)
    delims <- gregexpr(delim, x)[[1L]]
    endpts <- rep.int(str_len, nlines)
    for (i in seq_len(nlines - 1L)) endpts[i] <- max(delims[delims < max_len*i])
    startpts <- c(1L, endpts[-nlines] + 1L)
    lines <- paste(Map(function(a, b) substr(x, a, b), startpts, endpts), '\n')
    x <- trimws(paste(lines, collapse=''), 'right')
  }
  return(x)
}

#' Subset graphs based on a given logical condition
#'
#' \code{subset_graph} will subset a given graph based on the given logical
#' condition(s). This can be a \dQuote{simple} logical equation, or can include
#' combinations of \emph{AND} and \emph{OR} (i.e., \code{&} and \code{|}).
#'
#' @param g A graph object
#' @param condition Character string specifying an equation for which vertices to
#'   keep
#' @keywords internal
#' @return A graph object

subset_graph <- function(g, condition) {
  stopifnot(nzchar(condition))

  # Function for creating the condition string to later subset the graph
  get_cond_string <- function(orig) {
    spec <- '\\s\\&\\s|\\s\\|\\s'  # Splits are either " & " or " | "
    conditions <- strsplit(orig, split=spec)[[1L]]
    if (length(conditions) > 1L) {  # Multiple conditions
      if (isFALSE(grepl(spec, orig))) {
        stop('Logical operators must be surrounded by spaces!')
      }
      nchars <- cumsum(nchar(conditions))
      endpts <- nchars + seq.int(from=2L, by=3L, length.out=length(nchars))
      splits <- vapply(endpts, function(x) substr(orig, start=x, stop=x), character(1L))
      conditions <- trimws(conditions) # Remove unnecessary whitespace

      cond.string <- paste(vapply(seq_along(conditions), function(x)
                                  paste0('V(g)$', conditions[x], splits[x]), character(1L)),
                           collapse='')
    } else {
      cond.string <- paste0('V(g)$', conditions)
    }
    return(cond.string)
  }

  # Handle when logical expressions are separated by parentheses
  if (isTRUE(grepl('\\(.*\\&.*\\)', condition)) || isTRUE(grepl('\\(.*\\|.*\\)', condition))) {
    subs <- strsplit(condition, split='\\) & \\(')[[1L]]
    subs <- as.list(trimws(subs, whitespace='[\\(\\) ]'))
    cond.strings <- vapply(subs, get_cond_string, character(1L))
    cond.string <- paste0('(', paste(cond.strings, collapse=') & ('), ')')
  } else {
    cond.string <- get_cond_string(condition)
  }

  cond <- eval(parse(text=cond.string))
  if (sum(cond, na.rm=TRUE) == 0L) {
    warning('No vertices meet criteria! No graph created')
    g <- inds <- NULL
  } else {
    inds <- which(cond)
    cond <- setdiff(seq_len(vcount(g)), inds)
    orig.class <- class(g)
    g <- delete.vertices(g, cond)
    class(g) <- orig.class
  }
  list(g=g, inds=inds)
}

#' Transform a vector to have a different range
#'
#' \code{vec.transform} takes a vector and transforms it to have a new range,
#' given the input, or the default values of [0, 1].
#'
#' @param min.val the minimum value of the new range
#' @param max.val the maximum value of the new range
#'
#' @keywords internal
#' @return A vector of the transformed input.
#' @rdname utils

vec.transform <- function(x, min.val=0, max.val=1) {
  diffrange <- diff(range(x, na.rm=TRUE))
  if (diffrange == 0) {
    out <- rep_len(max.val, length(x))
  } else {
    out <- ((x - min(x, na.rm=TRUE)) * (max.val - min.val) / diffrange) + min.val
  }
  return(out)
}
