#' Plot vertex-level graph measures at a single density or threshold
#'
#' This function creates boxplots of a single vertex-level graph measure at a
#' single density or threshold, grouped by the variable specified by
#' \code{facet.by} (e.g., \emph{lobe} or \emph{network}).
#'
#' @param tidy.dt A ``tidied'' \code{data.table} of vertex-level graph measures
#' @param facet.by Character string indicating whether the data should be
#'   plotted separately by a certain variable (default: 'lobe')
#' @param measure A character string of the graph measure to plot (default:
#'   'btwn.cent')
#' @param show.points Logical indicating whether or not to show individual data
#'   points (default: FALSE)
#' @param ylabel A character string for the y-axis label
#' @export
#'
#' @return A \code{ggplot} object
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' ggp.btwn <- plot_vertex_measures(dt.V.tidy, facet.by='network',
#'   measure='E.nodal')
#' }

plot_vertex_measures <- function(tidy.dt,
                                 facet.by='lobe', measure='btwn.cent',
                                 show.points=FALSE, ylabel=NULL) {
  variable <- Group <- value <- NULL

  if (is.null(ylabel)) ylabel <- measure
  my.plot <- ggplot(tidy.dt[variable == measure],
                    aes(x=Group, y=value, col=Group)) +
    geom_boxplot(outlier.shape=NA) +
    facet_wrap(as.formula(paste('~', facet.by)), scales='free_y') +
    ylab(ylabel) +
    theme(legend.position='none')

  if (isTRUE(show.points)) {
    my.plot <- my.plot + geom_point(position=position_jitter(width=0.2, height=0))
  }
  return(my.plot)
}
