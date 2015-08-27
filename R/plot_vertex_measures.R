#' Plot vertex-level graph measures at a single density
#'
#' This function creates boxplots of a single vertex-level graph measure at a
#' single density, grouped by \emph{lobe}; each lobe has a separate \emph{facet}
#' in a \code{ggplot} object.
#'
#' @param dat A ``tidied'' data table of vertex-level graph measures
#' @param cur.density A numeric indicating the graph density
#' @param measure A character string of the graph measure to plot (default:
#' 'btwn.cent')
#' @param ylabel A character string for the y-axis label
#' @export
#'
#' @return A \code{ggplot} object
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' ggp.btwn <- plot_vertex_measures(net.meas.tidy, densities[N], 'btwn.cent')
#' }

plot_vertex_measures <- function(dat, cur.density=0.1, measure='btwn.cent',
                                 ylabel=NULL) {
  variable <- Group <- value <- NULL

  if (is.null(ylabel)) {
    ylabel <- measure
  }
  my.plot <- ggplot(dat[abs(density - cur.density) <= 0.001 & variable == measure],
                    aes(x=Group, y=value, col=Group)) +
    geom_boxplot(outlier.shape=NA) +
    geom_point(position=position_jitter(width=0.2, height=0)) +
    facet_wrap(~ lobe, scales='free_y') +
    ylab(ylabel) +
    theme(legend.position='none')

  return(my.plot)
}
