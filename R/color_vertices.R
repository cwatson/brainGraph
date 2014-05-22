#' Color the vertices of a graph based on community membership.
#'
#' This function takes the community membership of a given graph, and assigns
#' to the vertices a specific color. Currently only works for up to 5
#' communities (colored: red, green, blue, yellow, magenta). Vertices that
#' belong to small communities (defined by me as degree < 2) are white.
#'
#' @param adj.graph The adjacency graph to get its vertices colored
#' @param comm The community object returned from community.measures
#' @export
#'

color.vertices <- function(adj.graph, comm) {
  cols <- comm$vcolors[comm$community$membership]
}
