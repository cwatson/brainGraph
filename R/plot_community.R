#' Plot a given community
#'
#' This function plots a single community, if you want to see e.g. the first
#' community and nothing else.
#'
#' @param g The graph
#' @param n The community to be plotted (integer)
#' @param ... Other parameters (passed to \code{\link{plot.adj}})
#' @export

plot.community <- function(g, n, ...) {
  M <- max(V(g)$comm)
  if (n > M) {
    stop(sprintf('Community number out of bounds; max = %i.', M))
  }

  v <- as.numeric(names(rev(sort(table(V(g)$comm)))[n]))
  v <- V(g)[V(g)$comm==v]
  g.sub <- induced.subgraph(g, v)

  args <- list(...)
  if (hasArg(vertex.label)) {
    vertex.label <- args[['vertex.label']]
    if(!is.na(vertex.label)) {
      vertex.label <- V(g.sub)$name
      vertex.label.cex <- 0.75
    } else {
      vertex.label.cex <- NA
    }
  } else {
    vertex.label <- NA
    vertex.label.cex <- NA
  }

  plot.adj(g.sub,
           vertex.color=V(g.sub)$color, edge.color=E(g.sub)$color,
           vertex.label=vertex.label, vertex.label.cex=vertex.label.cex, ...)
  par(new=T, mar=c(5, 0, 3, 0)+0.1)
}
