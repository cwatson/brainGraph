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
#' @param legend.pos A numeric vector indicating the legend position (default:
#'   c(1, 0))
#' @param vline Numeric; required to plot a dashed vertical line (default: NULL)
#' @param level.names Character vector of facet label names, if you wish to
#'   change them (default: NULL)
#' @param exclude Character vector of variables to exclude (default: NULL)
#' @param perms A \code{\link{data.table}} of permutation group differences
#'   (default: NULL)
#' @param g A list of lists of \code{igraph} graph objects; required if
#'   \emph{perms} is provided (default: NULL)
#' @param alt Character vector of alternative hypotheses; required if
#'   \emph{perms} is provided (default: NULL)
#' @export
#'
#' @return A \code{\link[ggplot2]{ggplot}} object
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

plot_global <- function(tidy.dt, legend.pos=c(1, 0), vline=NULL,
                        level.names=NULL, exclude=NULL, perms=NULL,
                        g=NULL, alt=NULL) {
  sig <- trend <- yloc <- value <- variable <- Group <- NULL

  subDT <- copy(tidy.dt)
  # Add asterisks if a data.table of permutation values is provided
  if (!is.null(perms)) {
    if (is.null(g)) stop(paste0('Must provide a list of (lists of) graphs'))
    if (is.null(alt)) stop(paste0('Must provide a vector of alternative hypotheses'))

    # Add necessary columns for plotting annotations
    subDT[, sig := '']
    subDT[, trend := '']
    subDT[, yloc := round(min(value) - 0.05 * diff(range(value)), 3), by=variable]
    vars <- setdiff(names(perms), 'density')
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

  p <- ggplot(subDT, aes(x=density, y=value, col=Group)) +
    geom_line() +
    facet_wrap(~ variable, scales='free_y') +
    theme(legend.position=legend.pos, legend.justification=legend.pos,
          legend.background=element_rect(fill='gray90', size=0.5))
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
