#' Plot vertex-level graph measures at a single density or threshold
#'
#' Creates boxplots of a single vertex-level graph measure at a single density
#' or threshold, grouped by the variable specified by \code{group.by} and
#' optionally faceted by another variable (e.g., \emph{lobe} or \emph{network}).
#'
#' @param measure A character string of the graph measure to plot
#' @param facet.by Character string indicating the variable to facet by (if
#'   any). Default: \code{NULL}
#' @param group.by Character string indicating which variable to group the data
#'   by. Default: \code{'Group'}
#' @param type Character string indicating the plot type. Default:
#'   \code{'violin'}
#' @param show.points Logical indicating whether or not to show individual data
#'   points (default: FALSE)
#' @param ylabel A character string for the y-axis label
#' @param ... Other arguments passed to \code{\link[ggplot2]{geom_boxplot}} or
#'   \code{\link[ggplot2]{geom_violin}}
#' @inheritParams plot_brainGraph_multi
#' @export
#'
#' @return A \code{ggplot} object
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' p.deg <- plot_vertex_measures(g[[1]], facet.by='network', measure='degree')
#' }

plot_vertex_measures <- function(g.list, measure, facet.by=NULL, group.by='Group',
                                 type=c('violin', 'boxplot'),
                                 show.points=FALSE, ylabel=measure, ...) {
  variable <- value <- NULL

  if (!inherits(g.list, 'brainGraphList')) try(g.list <- as_brainGraphList(g.list))
  DT <- vertex_attr_dt(g.list)
  stopifnot(measure %in% names(DT), group.by %in% names(DT))
  idvars <- c('atlas', 'modality', 'weighting', 'Study.ID', 'Group', 'threshold',
              'density', 'region', 'lobe', 'hemi', 'class', 'network')
  idvars <- idvars[which(idvars %in% names(DT))]
  DT.m <- melt(DT, id.vars=idvars)
  setnames(DT.m, group.by, 'group.by')

  type <- match.arg(type)
  p <- ggplot(DT.m[variable == measure], aes(x=group.by, y=value, fill=group.by))
  p <- if (type == 'violin') p + geom_violin(...) else p + geom_boxplot(...)

  if (!is.null(facet.by)) {
    stopifnot(facet.by %in% names(DT))
    p <- p + facet_wrap(as.formula(paste('~', facet.by)), scales='free_y')
  }
  if (isTRUE(show.points)) {
    p <- p + geom_jitter(position=position_jitter(width=0.1, height=0))
  }
  p <- p + labs(x=group.by, y=ylabel) + theme(legend.position='none')

  return(p)
}
