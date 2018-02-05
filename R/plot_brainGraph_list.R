#' Write PNG files for a list of graphs
#'
#' This function takes a list of \code{igraph} graph objects and plots them over
#' an axial slice of the brain. A \code{png} file is written for each element of
#' the list, which can be joined as a \code{gif} or converted to video using a
#' tool outside of \code{R}.
#'
#' You can choose to highlight edge differences between subsequent list
#' elements; in this case, new/different edges are colored pink.
#'
#' This function may be particularly useful if the graph list contains graphs of
#' a single subject group at incremental densities, or if the graph list
#' contains graphs of each subject in a group.
#'
#' @param g.list A list of \code{igraph} graph objects
#' @param fname.base A character string specifying the base of the filename for
#'   \emph{png} output
#' @param diffs A logical, indicating whether or not to highlight edge
#'   differences (default: \code{FALSE})
#' @param ... Other parameters (passed to \code{\link{plot_brainGraph}})
#' @export
#'
#' @family Plotting functions
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

plot_brainGraph_list <- function(g.list, fname.base, diffs=FALSE, ...) {

  for (i in seq_along(g.list)) {
    png(filename=sprintf('%s_%03d.png', fname.base, i))

    plot(g.list[[i]], main=g.list[[i]]$Group, ...)
    if (isTRUE(diffs)) {
      if (i > 1) {
        g.diff <- graph.difference(g.list[[i]], g.list[[i-1]])
        class(g.diff) <- c('brainGraph', class(g.diff))
        if (hasArg('edge.color')) {
          ecols <- list(...)$edge.color
        } else {
          ecols <- rep('deeppink', ecount(g.diff))
        }
        ewidth <- 5
        plot(g.diff, add=TRUE, vertex.label=NA, vertex.shape='none',
                        edge.width=ewidth, edge.color=ecols, mni=FALSE)
      }
    }
    dev.off()
  }
}
