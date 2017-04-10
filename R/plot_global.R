#' Plot global graph measures across densities
#'
#' Create a faceted line plot of global graph measures across a range of graph
#' densities. Given a "tidied" \code{data.table}, you can choose to insert a
#' dashed vertical line at a density of interest, rename the variable levels
#' (which become the facet titles), exclude certain variables, and include a
#' \code{data.table} of permutation data to add asterisks indicating signficant
#' group differences.
#'
#' @param tidy.dt A \code{\link{data.table}} that has been "tidied", containing
#'   global graph measures for all densities and subject groups
#' @param xvar A character string indicating whether the variable of
#' interest is "density" or "threshold" (e.g. with DTI data)
#' @param vline Numeric; required to plot a dashed vertical line (default: NULL)
#' @param level.names Character vector of facet label names, if you wish to
#'   change them (default: NULL)
#' @param exclude Character vector of variables to exclude (default: NULL)
#' @param perms A \code{\link{data.table}} of permutation group differences
#'   (default: NULL)
#' @param g A list of lists of \code{igraph} graph objects; required if
#'   \emph{perms} is provided (default: NULL)
#' @param alt Character vector of alternative hypotheses; required if
#'   \emph{perms} is provided, but defaults to "two.sided" for all variables
#' @export
#'
#' @return A \code{\link[ggplot2]{ggplot}} object
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

plot_global <- function(tidy.dt, xvar=c('density', 'threshold'), vline=NULL,
                        level.names=NULL, exclude=NULL, perms=NULL,
                        g=NULL, alt=NULL) {
  sig <- trend <- yloc <- value <- variable <- Group <- threshold <- NULL

  subDT <- copy(tidy.dt)
  # Add asterisks if a data.table of permutation values is provided
  if (!is.null(perms)) {
    stopifnot(!is.null(g))

    # Add necessary columns for plotting annotations
    subDT[, c('sig', 'trend') := '']
    subDT[, yloc := round(min(value) - 0.05 * diff(range(value)), 3), by=variable]
    vars <- setdiff(names(perms), 'density')
    if (is.null(alt)) alt <- rep('two.sided', length(vars))
    for (i in seq_along(vars)) {
      dt <- plot_perm_diffs(g[[1]], g[[2]], perms, vars[i],
                            groups=subDT[, unique(Group)], alternative=alt[i])$dt
      subDT[variable == vars[i], sig := dt$sig]
      subDT[variable == vars[i], trend := dt$trend]
      subDT[variable == vars[i], yloc := dt$yloc]
    }
  }

  # Change the factor level names if desired
  subDT <- subDT[!variable %in% exclude]
  if (!is.null(level.names)) levels(subDT$variable) <- level.names

  xvar <- match.arg(xvar)
  if (xvar == 'density') {
    p <- ggplot(subDT, aes(x=density, y=value, col=Group))
  } else if (xvar == 'threshold') {
    p <- ggplot(subDT, aes(x=threshold, y=value, col=Group)) + scale_x_reverse()
  }

  if ('Study.ID' %in% names(subDT)) {
    p <- p + stat_smooth(method='gam', formula=y~s(x))
  } else {
    p <- p + geom_line()
  }
  p <- p +
    facet_wrap(~ variable, scales='free_y') +
    theme(legend.position='bottom',
          axis.text.x=element_text(angle=45, hjust=1))

  # Add a dashed vertical line at the density of interest
  if (!is.null(vline)) p <- p + geom_vline(xintercept=vline, lty=2, col='grey60')

  # Add in the asterisks if permutation values were given
  if (!is.null(perms)) {
    p <- p +
      geom_text(aes(y=yloc, label=sig), col='red', size=3) +
      geom_text(aes(y=yloc, label=trend), col='blue', size=3)
  }
  return(p)
}
