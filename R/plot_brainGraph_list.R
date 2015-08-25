#' Write PNG files for a list of graphs
#'
#' This function takes a list of \code{igraph} graph objects and plots them over
#' an axial slice of the brain. A \emph{png} file is written for each element of
#' the list, which can be joined as a \emph{gif} or converted to video.
#'
#' You can choose to highlight edge differences between subsequent list
#' elements, and whether to color vertices by \emph{lobe}, \emph{community}
#' membership, or \emph{lightblue} (the default). By default, the vertex sizes
#' are equal to vertex degree, and max out at 20.
#'
#' This function may be particularly useful if the graph list contains graphs of
#' a single subject group at incremental densities, or if the graph list
#' contains graphs of each subject in a group.
#'
#' @param g.list A list of \code{igraph} graph objects
#' @param fname.base A character string specifying the base of the filename for
#' \emph{png} output
#' @param diffs A logical, indicating whether or not to highlight edge
#' differences (default: FALSE)
#' @param cols A character string indicating how to color the vertices (default:
#' 'none')
#' @export
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

plot_brainGraph_list <- function(g.list, fname.base, diffs=FALSE,
                                 cols=c('none', 'lobe', 'comm')) {

  for (i in seq_along(g.list)) {
    png(filename=sprintf('%s_%03d.png', fname.base, i))
    cols <- match.arg(cols)
    if (cols == 'none') {
      vcols <- 'lightblue'
      ecols <- 'red'
    } else if (cols == 'lobe') {
      vcols <- V(g.list[[i]])$color.lobe
      ecols <- E(g.list[[i]])$color.lobe
    } else if (cols == 'comm') {
      vcols <- V(g.list[[i]])$color.comm
      ecols <- E(g.list[[i]])$color.comm
    }
    plot_brainGraph_mni('axial')
    plot_brainGraph(g.list[[i]], vertex.label=NA,
                    vertex.size=pmin(V(g.list[[i]])$degree, 20),
                    vertex.color=vcols,
                    edge.width=1, edge.color=ecols)

    if (isTRUE(diffs)) {
      if (i > 1) {
        g.diff <- graph.difference(g.list[[i]], g.list[[i-1]])
        if (cols == 'lobe') {
          ecols <- E(g.diff)$color.lobe
        } else if (cols == 'comm') {
          ecols <- E(g.diff)$color.comm
        } else {
          ecols <- 'deeppink'
        }
        ewidth <- 5
        plot_brainGraph(g.diff, add=T, vertex.label=NA, vertex.shape='none',
                        edge.width=ewidth, edge.color=ecols)
      }
    }
    dev.off()
  }
}
