#' Contract graph vertices based on brain lobe and hemisphere
#'
#' This function creates a new graph after merging multiple vertices based on
#' brain lobe and hemisphere membership. The new vertex size is equal to the
#' number of vertices in each lobe. The x- and y- coordinates of the new
#' vertices are equal to the mean of the lobe vertices of the original graph.
#' The new edge weight is equal to the number of inter-lobular connections of
#' the original graph.
#'
#' @param g The graph to contract
#' @export
#'
#' @return A new graph
#'
#' @seealso \code{\link[igraph]{contract.vertices}}

graph.contract.brain <- function(g) {
  g.sub <- contract.vertices(g, V(g)$lobe.hemi)
  E(g.sub)$weight <- 1
  g.sub <- simplify(g.sub)
  V(g.sub)$x <- sapply(1:max(V(g)$lobe.hemi),
                      function(m) mean(V(g)[V(g)$lobe.hemi==m]$x))
  V(g.sub)$y <- sapply(1:max(V(g)$lobe.hemi),
                      function(m) mean(V(g)[V(g)$lobe.hemi==m]$y))
  V(g.sub)$size <- sapply(1:max(V(g)$lobe.hemi),
                      function(m) mean(V(g)[V(g)$lobe.hemi==m]$degree))
  vcols <- unique(V(g)$color.lobe[order(V(g)$lobe)])
  vcols <- rep(vcols, each=2)
  V(g.sub)$color <- vcols
  g.sub
}
