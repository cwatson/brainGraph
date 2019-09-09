#' Write files to be used for visualization with BrainNet Viewer
#'
#' Write the \code{.node} and \code{.edge} files necessary for visualization
#' with the BrainNet Viewer software (see Reference below).
#'
#' @details For the \code{.node} file, there are 6 columns:
#' \itemize{
#' \item \emph{Column 1}: x-coordinates
#' \item \emph{Column 2}: y-coordinates
#' \item \emph{Column 3}: z-coordinates
#' \item \emph{Column 4}: Vertex color
#' \item \emph{Column 5}: Vertex size
#' \item \emph{Column 6}: Vertex label
#' }
#' The \code{.edge} file is the graph's associated adjacency matrix; a weighted
#' adjacency matrix can be returned by using the \code{edge.wt} argument.
#'
#' @param g The \code{igraph} graph object of interest
#' @param vcolor Character string indicating how to color the vertices (default:
#'   \code{'none'})
#' @param vsize Character string indicating what size the vertices should be;
#'   can be any vertex-level attribute (default: \code{'constant'})
#' @param edge.wt Character string indicating the edge attribute to use to
#'   return a weighted adjacency matrix (default: \code{NULL})
#' @param file.prefix Character string for the basename of the \code{.node} and
#'   \code{.edge} files that are written
#' @export
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Xia, M. and Wang, J. and He, Y. (2013). BrainNet Viewer: a
#'   network visualization tool for human brain connectomics. \emph{PLoS One},
#'   \bold{8(7)}, e68910. \url{https://dx.doi.org/10.1371/journal.pone.0068910}
#' @examples
#' \dontrun{
#' write_brainnet(g, vcolor='community', vsize='degree', edge.wt='t.stat')
#' }

write_brainnet <- function(g, vcolor='none', vsize='constant', edge.wt=NULL, file.prefix='') {
  x.mni <- y.mni <- z.mni <- NULL
  stopifnot(is.brainGraph(g))

  atlas.dt <- get(g$atlas)
  coords <- round(atlas.dt[, cbind(x.mni, y.mni, z.mni)])

  vnames <- vertex_attr_names(g)
  if (vcolor == 'none') {
    cols <- rep(1, vcount(g))
  } else {
    stopifnot(vcolor %in% vnames)
    cols <- as.numeric(factor(vertex_attr(g, vcolor)))
  }

  if (vsize == 'constant') {
    size <- 5
  } else {
    stopifnot(vsize %in% vnames)
    size <- vertex_attr(g, vsize)
  }

  if (file.prefix == '') {
    nodefile <- paste0(quote(g), '_', vsize, '_', vcolor, '.node')
    edgefile <- paste0(quote(g), '.edge')
  } else {
    nodefile <- paste0(file.prefix, '.node')
    edgefile <- paste0(file.prefix, '.edge')
  }
  fwrite(cbind(coords, cols, size, V(g)$name),
         file=nodefile, quote=FALSE, sep='\t', col.names=FALSE)

  fwrite(as_adj(g, sparse=FALSE, attr=edge.wt),
         file=edgefile, quote=FALSE, sep='\t', col.names=FALSE)
}
