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
#' @return A data table with 18-19 columns and row number equal to the number of
#' vertices in the graph
#' @seealso \code{\link[igraph]{vertex_attr}, \link[igraph]{vertex_attr_names},
#' \link[igraph]{as_data_frame}}

vertex_attr_dt <- function(g, Group=NULL) {
  lobe <- name <- NULL
  atlas.dt <- eval(parse(text=data(list=g$atlas)))
  net.meas <- data.table(density=g$density,
                         region=V(g)$name,
                         lobe=atlas.dt[, levels(lobe)][V(g)$lobe],
                         hemi=V(g)$hemi,
                         degree=V(g)$degree,
                         knn=V(g)$knn,
                         btwn.cent=V(g)$btwn.cent,
                         hubs=V(g)$hubs,
                         ev.cent=V(g)$ev.cent,
                         subg.cent=V(g)$subgraph.cent,
                         lev.cent=V(g)$lev.cent,
                         coreness=V(g)$coreness,
                         trans=V(g)$transitivity,
                         E.local=V(g)$E.local,
                         E.nodal=V(g)$E.nodal,
                         PC=V(g)$PC,
                         z=V(g)$z.score,
                         vulnerability=V(g)$vulnerability,
                         asymm=V(g)$asymm,
                         eccentricity=V(g)$eccentricity)

  if ('strength' %in% vertex_attr_names(g)) net.meas$strength <- V(g)$strength
  if ('name' %in% graph_attr_names(g)) net.meas$subject <- g$name
  if (!is.null(Group)) {
    net.meas$Group <- Group
    setkey(net.meas, 'region', 'lobe', 'hemi', Group)
  } else {
    setkey(net.meas, 'region', 'lobe', 'hemi')
  }

  return(net.meas)
}
