#' Create a data table with graph vertex measures
#'
#' This is just a helper function that creates a data table in which each row is
#' a vertex and each column is a different network measure (degree, centrality,
#' etc.).
#'
#' @param g An igraph graph object
#' @param Group A character string indicating group membership (default:NULL)
#' @export
#'
#' @return A data table with 14 columns and row number equal to the number of
#' vertices in the graph
#' @seealso \code{\link[igraph]{vertex_attr}, \link[igraph]{vertex_attr_names}}

vertex_attr_dt <- function(g, Group=NULL) {
  lobe <- name <- NULL
  atlas.dt <- eval(parse(text=data(list=g$atlas)))
  net.meas <- data.table(density=g$density,
                         region=V(g)$name,
                         lobe=atlas.dt[match(V(g)$name, atlas.dt[, name])][, as.character(lobe)],
                         hemi=V(g)$hemi,
                         deg=V(g)$degree,
                         btwn.cent=V(g)$btwn.cent,
                         ev.cent=V(g)$ev.cent,
                         subg.cent=V(g)$subgraph.cent,
                         lev.cent=V(g)$lev.cent,
                         coreness=V(g)$coreness,
                         trans=V(g)$transitivity,
                         E.local=V(g)$E.local,
                         E.nodal=V(g)$E.nodal,
                         PC=V(g)$PC,
                         z=V(g)$z.score)

  if (!is.null(Group)) {
    net.meas$Group <- Group
    setkey(net.meas, 'region', 'lobe', 'hemi', Group)
  } else {
    setkey(net.meas, 'region', 'lobe', 'hemi')
  }

  return(net.meas)
}
