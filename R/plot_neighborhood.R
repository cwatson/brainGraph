#' Plot the neighborhood of a given vertex.
#'
#' This function plots only the neighborhood of a single given vertex, e.g. if
#' you want to see the neighborhood of the left precuneus and nothing else.
#'
#' @param g The graph to be plotted
#' @param n Index (integer) of the vertex whose neighborhood is to be plotted
#' @param ... Other parameters (passed to \code{\link{plot.adj}})
#' @export

plot.neighborhood <- function(g, n, ...) {
  g.sub <- graph.neighborhood(g, nodes=n, order=1)[[1]]

  inds <- sort(c(n, neighbors(g, n)))

  fargs <- list(...)
  if (hasArg(vertex.label)) {
    vertex.label <- fargs$vertex.label
    if (!is.na(vertex.label)) {
      vertex.label <- vertex.label[inds]
      vertex.label.cex <- 0.75
    } else {
      vertex.label.cex <- NA
    }
  }

  if (hasArg(vertex.color)) {
    vertex.color <- fargs$vertex.color
    if (length(vertex.color) > 1) {
      vertex.color <- vertex.color[inds]
    }
  }

  if (hasArg(edge.color)) {
    edge.color <- fargs$edge.color
    if (length(edge.color) > 1) {
      es <- get.edge.ids(g, combn(inds, 2))
      es <- es[which(es != 0)]
      edge.color <- edge.color[es]
    }
  }

  if (hasArg(vertex.size)) {
    vertex.size <- fargs$vertex.size
    if (length(vertex.size) > 1) {
      vertex.size <- vertex.size[inds]
    }
  }

  Nv <- vcount(g.sub)
  Ne <- ecount(g.sub)
  g.density <- round((2 * Ne) / (Nv * (Nv - 1)), digits=3)
  plot.adj(g.sub,
           vertex.size=vertex.size, vertex.color=vertex.color,
           edge.color=edge.color, vertex.label=vertex.label, 
           vertex.label.cex=vertex.label.cex, ...)
  par(new=T, mar=c(5, 0, 3, 0)+0.1)
  title(paste('Neighborhood of', V(g)[n]$name), col.main='white',
        sub=paste('# vertices: ', Nv, '# edges: ', Ne, '\n',
                   'Density: ', g.density),
        col.sub='white')
}
