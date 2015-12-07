#' Create a data table with graph global measures
#'
#' This is just a helper function that takes a list of graphs and creates a data
#' table of global measures for each graph, ordered by graph density.
#'
#' @param g.list A list of igraph graph objects
#' @param group A character string indicating group membership (default:NULL)
#' @export
#'
#' @return A data table with several columns (equal to the number of graph
#' attributes) and row number equal to the number of graphs in the input list
#' @seealso \code{\link[igraph]{graph_attr}, \link[igraph]{graph_attr_names}}

graph_attr_dt <- function(g.list, group=NULL) {
  Group <- NULL
  inds <- which(sapply(graph_attr(g.list[[1]]), class) %in% c('numeric', 'integer'))
  g.attr.names <- graph_attr_names(g.list[[1]])[inds]
  g.dt <- as.data.table(sapply(g.attr.names, function(x)
                               sapply(g.list, function(y)
                                      round(graph_attr(y, x), 4))))

  if ('name' %in% graph_attr_names(g.list[[1]])) {
    g.dt$Study.ID <- sapply(g.list, function(x) x$name)
  }
  if ('modality' %in% graph_attr_names(g.list[[1]])) {
    g.dt$modality <- sapply(g.list, function(x) x$modality)
  }
  if (is.null(group)) {
    if ('Group' %in% graph_attr_names(g.list[[1]])) {
      g.dt$Group <- sapply(g.list, function(x) x$Group)
      setkey(g.dt, Group, density)
    }
    setkey(g.dt, density)
  } else {
    g.dt$Group <- group
    setkey(g.dt, Group, density)
  }

  return(g.dt)
}
