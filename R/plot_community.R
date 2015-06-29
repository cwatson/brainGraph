#' Plot a given community
#'
#' This function plots a single (or multiple) community, if you want to see e.g.
#' the first community and nothing else, or the first two, etc.
#'
#' @param g The graph
#' @param n The community to be plotted (integer or vector)
#' @param ... Other parameters (passed to \code{\link{plot.adj}})
#' @export

plot_community <- function(g, n, ...) {
  M <- max(V(g)$comm)
  if (any(n > M)) {
    stop(sprintf('Community number out of bounds; max = %i.', M))
  }

  comm <- as.numeric(names(rev(sort(table(V(g)$comm)))[n]))
  members <- which(V(g)$comm %in% comm)
  g.sub <- induced.subgraph(g, members)

  fargs <- list(...)
  if (hasArg(vertex.label)) {
    vertex.label <- fargs$vertex.label
    if(!is.na(vertex.label[1])) {
      vertex.label <- V(g.sub)$name
      vertex.label.cex <- 0.75
    } else {
      vertex.label.cex <- NA
    }
  }

  if (hasArg(vertex.label.dist)) {
    vertex.label.dist <- fargs$vertex.label.dist
    vertex.label.dist <- vertex.label.dist[members]
  }

  if (hasArg(vertex.color)) {
    vertex.color <- fargs$vertex.color
    if (length(vertex.color) > 1) {
      vertex.color <- vertex.color[members]
    }
  }

  if (hasArg(edge.color)) {
    edge.color <- fargs$edge.color
    if (length(edge.color) > 1) {
      es <- get.edge.ids(g, combn(members, 2))
      es <- es[which(es != 0)]
      edge.color <- edge.color[es]
    }
  }

  if (hasArg(vertex.size)) {
    vertex.size <- fargs$vertex.size
    if (length(vertex.size) > 1) {
      vertex.size <- vertex.size[members]
    }
  }

  g.density <- round(graph.density(g), digits=3)
  plot.adj(g.sub,
           vertex.size=vertex.size, vertex.color=vertex.color,
           edge.color=edge.color, vertex.label=vertex.label,
           vertex.label.cex=vertex.label.cex,
           vertex.label.dist=vertex.label.dist, ...)
}
