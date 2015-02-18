#' Create a data table with graph global measures
#'
#' This is just a helper function that takes a list of graphs and creates a data
#' table of global measures for each graph, ordered by graph density.
#'
#' @param g.list A list of igraph graph objects
#' @param Group A character string indicating group membership (default:NULL)
#' @export
#'
#' @return A data table with 12 columns and row number equal to the number of
#' graphs in the input list
#' @seealso \code{\link[igraph]{graph_attr}, \link[igraph]{graph_attr_names}}

graph_attr_dt <- function(g.list, Group=NULL) {
  glob.meas <- data.table(density=sapply(g.list, function(x) x$density),
    conn.comp=sapply(g.list, function(x) x$conn.comp[1, 1]),
    num.tri=sapply(g.list, function(x) x$num.tri),
    diameter=sapply(g.list, function(x) x$diameter),
    Cp=sapply(g.list, function(x) x$Cp),
    Lp=sapply(g.list, function(x) x$Lp),
    assortativity=sapply(g.list, function(x) x$assortativity),
    assortativity.lobe=sapply(g.list, function(x) x$assortativity.lobe),
    assortativity.lobe.hemi=sapply(g.list, function(x) x$assortativity.lobe.hemi),
    E.global=sapply(g.list, function(x) x$E.global),
    E.local=sapply(g.list, function(x) x$E.local),
    mod=sapply(g.list, function(x) x$mod))

  if (!is.null(Group)) {
    glob.meas$Group <- Group
    setkey(glob.meas, Group, density)
  } else {
    setkey(glob.meas, density)
  }

  return(glob.meas)
}
