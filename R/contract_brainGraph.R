#' Contract graph vertices based on brain lobe and hemisphere
#'
#' Create a new graph after merging vertices within specified groups.
#' By default, groups are brain \emph{lobe} and \emph{hemisphere} membership.
#'
#' The \code{size} vertex-level attribute of the resultant graph is equal to the
#' number of vertices in each group. The x-, y-, and z-coordinates of the new
#' graph are equal to the mean coordinates of the vertices per group.
#' The new edge weights are equal to the number of inter-group connections of
#' the original graph.
#'
#' @param g A \code{brainGraph} graph object
#' @param vgroup Character string; the name of the vertex attribute to use when
#'   contracting the graph. Default: \code{'lobe.hemi'}
#' @export
#'
#' @return A new \code{brainGraph} graph object with vertex-level attributes
#'   representing the mean spatial coordinates, and vertex- and edge-level
#'   attributes of color names
#'
#' @seealso \code{\link[igraph]{contract}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

contract_brainGraph <- function(g, vgroup='lobe.hemi') {
  stopifnot(is.brainGraph(g), vgroup %in% vertex_attr_names(g))

  vattrs <- vertex_attr(g, vgroup)
  if (!is.numeric(vattrs)) vattrs <- as.numeric(factor(vattrs))
  g.sub <- contract(g, vattrs, vertex.attr.comb='concat')
  E(g.sub)$count <- 1
  g.sub <- simplify(g.sub, edge.attr.comb=list(count='sum', 'concat'))

  # Simplify coordinate- and color-related vertex and edge attributes
  for (nam in c('x', 'y', 'z')) {
    val <- vapply(vertex_attr(g.sub, nam), mean, numeric(1L))
    g.sub <- delete_vertex_attr(g.sub, nam)
    g.sub <- set_vertex_attr(g.sub, nam, value=val)
  }
  for (nam in grep('color', vertex_attr_names(g.sub), value=TRUE)) {
    vvals <- lapply(vertex_attr(g.sub, nam), unique)
    if (all(lengths(vvals) == 1L)) {
      g.sub <- delete_vertex_attr(g.sub, nam)
      g.sub <- set_vertex_attr(g.sub, nam, value=unlist(vvals))
    }
    evals <- lapply(edge_attr(g.sub, nam), unique)
    if (all(lengths(evals) == 1L)) {
      g.sub <- delete_edge_attr(g.sub, nam)
      g.sub <- set_edge_attr(g.sub, nam, value=unlist(evals))
    }
  }

  V(g.sub)$size <- vapply(V(g.sub), function(v) sum(E(g.sub)[v %--% V(g.sub)]$count), numeric(1L))

  class(g.sub) <- c('brainGraph', class(g.sub))
  return(g.sub)
}
