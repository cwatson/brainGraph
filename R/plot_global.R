#' Plot global graph measures across densities
#'
#' Create a faceted line plot of global graph measures across a range of graph
#' densities, calculated from a list of \code{brainGraphList} objects. This
#' requires that the variables of interest are graph-level attributes of the
#' input graphs.
#'
#' You can choose to insert a dashed vertical line at a specific
#' density/threshold of interest, rename the variable levels (which become the
#' facet titles), exclude variables, and include a \code{brainGraph_permute}
#' object of permutation data to add asterisks indicating significant group
#' differences.
#'
#' @inheritParams random_graphs
#' @param xvar A character string indicating whether the variable of
#' interest is \dQuote{density} or \dQuote{threshold} (e.g. with DTI data)
#' @param vline Numeric of length 1 specifying the x-intercept if you would like
#'   to plot a vertical dashed line (e.g., if there is a particular density of
#'   interest). Default: \code{NULL}
#' @param level.names Character vector of facet label names, if you wish to
#'   change them. Default: \code{NULL}
#' @param exclude Character vector of variables to exclude. Default: \code{NULL}
#' @param perms A \code{\link{data.table}} of permutation group differences
#' @param alt Character vector of alternative hypotheses; required if
#'   \emph{perms} is provided, but defaults to \dQuote{two.sided} for all
#'   variables
#' @export
#'
#' @return A \code{\link[ggplot2]{ggplot}} object
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

plot_global <- function(g.list, xvar=c('density', 'threshold'), vline=NULL,
                        level.names=NULL, exclude=NULL, perms=NULL, alt=NULL) {
  sig <- trend <- yloc <- value <- variable <- Group <- threshold <- NULL

  # Check if components are 'brainGraphList' objects
  matches <- vapply(g.list, inherits, logical(1), 'brainGraphList')
  if (any(!matches)) stop("Input must be a list of 'brainGraphList' objects.")

  DT <- rbindlist(lapply(g.list, graph_attr_dt))
  idvars <- c('atlas', 'modality', 'weighting', 'Study.ID', 'Group', 'threshold', 'density')
  idvars <- idvars[which(idvars %in% names(DT))]
  DT.m <- melt(DT, id.vars=idvars)
  DT.m <- droplevels(DT.m[!variable %in% exclude])

  # Add asterisks if a data.table of permutation values is provided
  if (!is.null(perms)) {
    DT.m[, c('sig', 'trend') := '']
    DT.m[, yloc := extendrange(value)[1L], by=variable]
    vars <- setdiff(names(perms$DT), 'densities')
    nvars <- length(vars)
    if (is.null(alt)) alt <- rep('two.sided', nvars)
    if (length(alt) < nvars) alt <- rep(alt, length.out=nvars)
    for (i in seq_along(vars)) {
      dt <- plot(perms, measure=vars[i], alternative=alt[i])[[1]]$data[variable == 'obs.diff']
      DT.m[variable == vars[i], sig := dt$sig]
      DT.m[variable == vars[i], trend := dt$trend]
      DT.m[variable == vars[i], yloc := dt$yloc]
    }
  }

  if (!is.null(level.names)) levels(DT.m$variable) <- level.names

  xvar <- match.arg(xvar)
  p <- switch(xvar,
              density=ggplot(DT.m, aes(x=density, y=value, col=Group)),
              threshold=ggplot(DT.m, aes(x=threshold, y=value, col=Group)) + scale_x_reverse())

  if ('Study.ID' %in% names(DT.m)) {
    p <- p + stat_smooth(method='gam', formula=y~s(x))
  } else {
    p <- p + geom_line()
  }
  p <- p +
    facet_wrap(~ variable, scales='free_y') +
    theme(legend.position='bottom', axis.text.x=element_text(angle=45, hjust=1))

  if (!is.null(vline)) p <- p + geom_vline(xintercept=vline, lty=2, col='grey60')
  if (!is.null(perms)) {
    p <- p +
      geom_text(aes(y=yloc, label=sig), col='red', size=3) +
      geom_text(aes(y=yloc, label=trend), col='blue', size=3)
  }
  return(p)
}
