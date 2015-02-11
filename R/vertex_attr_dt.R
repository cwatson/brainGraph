#' Create a data table listing vertex-wise measures for a graph
#'
#' This is just a helper function that creates a data table in which each row is
#' a vertex and each column is a different network measure (degree, centrality,
#' etc.).
#'
#' @param g An igraph graph object
#' @param Group A character string indicating group membership (default:NULL)
#' @export
#'
#' @return A data table with 12 columns and row number equal to the number of
#' vertices in the graph
#' @seealso \code{\link[igraph]{vertex_attr}, \link[igraph]{vertex_attr_names}}

vertex_attr_dt <- function(g, Group=NULL) {
  net.meas <- data.table(region=V(g)$name,
                         lobe=atlas.list$lobe[V(g)$lobe],
                         deg=V(g)$degree,
                         btwn.cent=V(g)$btwn.cent,
                         ev.cent=V(g)$ev.cent,
                         subg.cent=V(g)$subgraph.cent,
                         coreness=V(g)$coreness,
                         trans=V(g)$transitivity,
                         E.local=V(g)$E.local,
                         E.nodal=V(g)$E.nodal,
                         PC=V(g)$PC,
                         z=V(g)$z.score,)

  if (!is.null(Group)) {
    net.meas$Group <- Group
    setkey(net.meas, region, lobe, Group)
  } else {
    setkey(net.meas, region, lobe)
  }

  return(net.meas)
}
