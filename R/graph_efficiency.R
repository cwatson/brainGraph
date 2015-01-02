#' Calculate the global, local, or nodal efficiency of a graph.
#'
#' This function calculates the global efficiency of a graph or the local or
#' nodal efficiency of each vertex of a graph. The global efficiency is equal
#' to the mean of all nodal efficiencies.
#'
#' @param g The adjacency graph on which to calculate efficiency
#' @param type A character string; either 'local', 'nodal', or 'global'
#' @export
#'
#' @return A vector of the local efficiencies for each node of the graph (if 
#' type='local|nodal') or a number (if type='global').
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Latora V., Marchiori M. (2001) \emph{Efficient behavior of
#' small-world networks}. Phys Rev Lett, 87.19:198701.

graph.efficiency <- function(g, type='local') {
  if ('degree' %in% vertex_attr_names(g)) {
    degs <- V(g)$degree
  } else {
    degs <- degree(g)
  }

  if (type == 'local') {
    eff <- numeric(length(degs))
    nodes <- which(degs > 1)
  
    eff[nodes] <- simplify2array(mclapply(nodes, function(x) {
      neighbs <- neighbors(g, v=x)
      g.sub <- induced.subgraph(g, neighbs)
      Nv <- vcount(g.sub)
  
      paths <- shortest.paths(g.sub, weights=NA)
      paths <- paths[upper.tri(paths)]
      2 / Nv / (Nv - 1) * sum(1 / paths[paths != 0])
      }, mc.cores=detectCores())
    )
  } else if (type == 'nodal') {
    Nv <- length(degs)
    paths <- shortest.paths(g, weights=NA)
    eff <- simplify2array(mclapply(seq_len(Nv), function(x) {
      path <- paths[x, ]
      eff <- sum(1 / path[path != 0]) / (Nv - 1)
      }, mc.cores=detectCores())
    )
  } else {
    Nv <- length(degs)
    nodes <- V(g)[degs != 0]
    paths <- shortest.paths(g, v=nodes, to=nodes, weights=NA)
    paths <- paths[upper.tri(paths)]

    eff <- 2 / Nv / (Nv - 1) * sum(1 / paths[paths != 0])
  }
  return(eff)
}
