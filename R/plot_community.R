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

  comm <- as.numeric(names(rev(sort(table(V(g)$comm)))[n]))
  members <- which(V(g)$comm==comm)
  g.sub <- induced.subgraph(g, members)

  fargs <- list(...)
  if (hasArg(vertex.label)) {
    vertex.label <- fargs$vertex.label
    if(!is.na(vertex.label)) {
      vertex.label <- V(g.sub)$name
      vertex.label.cex <- 0.75
    } else {
      vertex.label.cex <- NA
    }
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

  Nv <- vcount(g.sub)
  Ne <- ecount(g.sub)
  g.density <- round((2 * Ne) / (Nv * (Nv - 1)), digits=3)
  plot.adj(g.sub,
           vertex.size=vertex.size, vertex.color=vertex.color,
           edge.color=edge.color, vertex.label=vertex.label,
           vertex.label.cex=vertex.label.cex, ...)
  par(new=T, mar=c(5, 0, 3, 0)+0.1)
  title(sub=paste('# vertices: ', Nv, '# edges: ', Ne, '\n',
                  'Density: ', g.density),
        col.sub='white')
}
