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
#'   not; can be 'none', 'lobe', or 'community'
#' @param node.size Character string indicating what size the vertices should be;
#'   can be any vertex-level attribute (default: 'constant')
#' @param edge.wt Character string indicating the edge attribute to use to
#'   return a weighted adjacency matrix
#' @param file.prefix Character string for the basename of the \emph{.node} and
#'   \emph{.edge} files that are written
#' @export
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Xia M, Wang J, He Y (2013). \emph{BrainNet Viewer: a network
#' visualization tool for human brain connectomics}. PLoS One, 8(7):e68910.
#' @examples
#' \dontrun{
#' write.brainnet(g, node.color='community', node.size='degree',
#'   edge.wt='t.stat')
#' }

write.brainnet <- function(g, node.color=c('none', 'community', 'lobe', 'component'),
                           node.size='constant', edge.wt=NULL, file.prefix='') {
  x.mni <- y.mni <- z.mni <- NULL
  stopifnot(is_igraph(g))

  atlas.dt <- eval(parse(text=g$atlas))
  coords.cur <- round(atlas.dt[, matrix(c(x.mni, y.mni, z.mni), ncol=3)])

  node.color <- match.arg(node.color)
  if (node.color == 'none') {
    color <- rep(1, vcount(g))
  } else if (node.color == 'community') {
    color <- V(g)$comm
  } else if (node.color == 'lobe') {
    color <- V(g)$lobe
  } else if (node.color == 'component') {
    color <- V(g)$comp
  }

  node.size <- match.arg(node.size)
  if (node.size == 'constant') {
    size <- 5
  } else {
    size <- vertex_attr(g, node.size)
  }

  if (file.prefix == '') {
    nodefile <- paste0(quote(g), '_', node.size, '_', node.color, '.node')
    edgefile <- paste0(quote(g), '.edge')
  } else {
    nodefile <- paste0(file.prefix, '.node')
    edgefile <- paste0(file.prefix, '.edge')
  }
  write.table(cbind(coords.cur,
                    color,
                    size,
                    V(g)$name),
              file=nodefile,
              row.names=F, col.names=F, sep='\t', quote=F)

  write.table(as_adj(g, sparse=F, attr=edge.wt),
              file=edgefile,
              row.names=F, col.names=F, sep='\t', quote=F)
}
