#' Create a data table with graph global measures
#'
#' This is just a helper function that takes a list of graphs and creates a data
#' table of global measures for each graph, ordered by graph density.
#'
#' @param g.list A list of igraph graph objects
#' @param Group A character string indicating group membership (default:NULL)
#' @export
#'
#' @return A data table with 14-15 columns and row number equal to the number of
#' graphs in the input list
#' @seealso \code{\link[igraph]{graph_attr}, \link[igraph]{graph_attr_names}}

graph_attr_dt <- function(g.list, Group=NULL) {
  glob.meas <- data.table(density=vapply(g.list, function(x) x$density, numeric(1)),
    conn.comp=vapply(g.list, function(x) x$conn.comp[1, 1], numeric(1)),
    num.tri=vapply(g.list, function(x) x$num.tri, numeric(1)),
    clique.num=vapply(g.list, function(x) x$clique.num, numeric(1)),
    diameter=vapply(g.list, function(x) x$diameter, numeric(1)),
    Cp=vapply(g.list, function(x) x$Cp, numeric(1)),
    Lp=vapply(g.list, function(x) x$Lp, numeric(1)),
    assortativity=vapply(g.list, function(x) x$assortativity, numeric(1)),
    assortativity.lobe=vapply(g.list, function(x) x$assortativity.lobe, numeric(1)),
    assortativity.lobe.hemi=vapply(g.list, function(x) x$assortativity.lobe.hemi, numeric(1)),
    E.global=vapply(g.list, function(x) x$E.global, numeric(1)),
    E.local=vapply(g.list, function(x) x$E.local, numeric(1)),
    mod=vapply(g.list, function(x) x$mod, numeric(1)),
    num.hubs=vapply(g.list, function(x) x$num.hubs, numeric(1)),
    asymm=vapply(g.list, function(x) x$asymm, numeric(1)),
    vulnerability=vapply(g.list, function(x) x$vulnerability, numeric(1)))

  if (!is.null(Group)) {
    glob.meas$Group <- Group
    setkey(glob.meas, Group, density)
  } else {
    setkey(glob.meas, density)
  }

  return(glob.meas)
}
