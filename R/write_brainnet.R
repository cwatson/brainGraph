#' Write files to be used for visualization with BrainNet Viewer
#'
#' This function will write the \emph{.node} and \emph{.edge} files necessary
#' for visualization with the BrainNet Viewer software (see Reference below).
#'
#' @details For the \emph{.node} file, there are 6 columns:
#' \itemize{
#' \item \emph{Column 1}: x-coordinates
#' \item \emph{Column 2}: y-coordinates
#' \item \emph{Column 3}: z-coordinates
#' \item \emph{Column 4}: Vertex color
#' \item \emph{Column 5}: Vertex size
#' \item \emph{Column 6}: Vertex label
#' }
#' The \emph{.edge} file is the graph's associated adjacency matrix; a weighted
#' adjacency matrix can be returned by using the \code{edge.wt} argument.
#'
#' @param g The \code{igraph} graph object of interest
#' @param node.color Character string indicating whether to color the vertices or
#'   not (default: \code{'none'})
#' @param node.size Character string indicating what size the vertices should be;
#'   can be any vertex-level attribute (default: \code{'constant'})
#' @param edge.wt Character string indicating the edge attribute to use to
#'   return a weighted adjacency matrix (default: \code{NULL})
#' @param file.prefix Character string for the basename of the \emph{.node} and
#'   \emph{.edge} files that are written
#' @export
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Xia, M and Wang, J and He, Y (2013). BrainNet Viewer: a network
#'   visualization tool for human brain connectomics. \emph{PLoS One},
#'   \bold{8(7)}, e68910. \url{https://dx.doi.org/10.1371/journal.pone.0068910}
#' @examples
#' \dontrun{
#' write_brainnet(g, node.color='community', node.size='degree',
#'   edge.wt='t.stat')
#' }

write_brainnet <- function(g, node.color='none', node.size='constant',
                           edge.wt=NULL, file.prefix='') {
  x.mni <- y.mni <- z.mni <- NULL
  stopifnot(is_igraph(g))

  atlas.dt <- get(g$atlas)
  coords.cur <- round(atlas.dt[, matrix(c(x.mni, y.mni, z.mni), ncol=3)])

  if (node.color == 'none') {
    color <- rep(1, vcount(g))
  } else {
    stopifnot(node.color %in% vertex_attr_names(g))
    color <- vertex_attr(g, node.color)
  }

  if (node.size == 'constant') {
    size <- 5
  } else {
    stopifnot(node.size %in% vertex_attr_names(g))
    size <- vertex_attr(g, node.size)
  }

  if (file.prefix == '') {
    nodefile <- paste0(quote(g), '_', node.size, '_', node.color, '.node')
    edgefile <- paste0(quote(g), '.edge')
  } else {
    nodefile <- paste0(file.prefix, '.node')
    edgefile <- paste0(file.prefix, '.edge')
  }
  write.table(cbind(coords.cur, color, size, V(g)$name),
              file=nodefile,
              row.names=F, col.names=F, sep='\t', quote=F)

  write.table(as_adj(g, sparse=F, attr=edge.wt),
              file=edgefile,
              row.names=F, col.names=F, sep='\t', quote=F)
}
