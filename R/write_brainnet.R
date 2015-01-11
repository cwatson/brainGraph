#' Write files to be used for visualization with BrainNet Viewer
#'
#' This function will write the 'node' and 'edge' files necessary for
#' visualization with the BrainNet Viewer software (see Reference below).
#'
#' @param g A graph
#' @param node.color Character string indicating whether to color the nodes or
#' not; can be 'none', 'lobe', or 'community'
#' @param node.size Character string indicating what size the nodes should be;
#' can be any vertex-level attribute or 'constant'
#' @export
#'
#' @seealso \code{\link{vertex_attr}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Xia M, Wang J, He Y (2013). \emph{BrainNet Viewer: a network
#' visualization tool for human brain connectomics}. PLoS One, 8(7):e68910.

write.brainnet <- function(g, node.color=c('none', 'community', 'lobe'),
                           node.size=c('constant', 'degree')) {
  atlas.list <- eval(parse(text=g$atlas))
  coords.cur <- round(atlas.list$brainnet.coords)

  node.color <- match.arg(node.color)
  if (node.color == 'none') {
    color <- rep(1, vcount(g))
  } else if (node.color == 'community') {
    color <- V(g)$comm
  } else if (node.color == 'lobe') {
    color <- V(g)$lobe
  }

  node.size <- match.arg(node.size)
  if (node.size == 'constant') {
    size <- 5
  } else {
    size <- vertex_attr(g, node.size)
  }

  # Write the 'node' and 'edge' files
  #------------------------------------------------------------------------------
  # For the 'node' file, there are 6 columns:
  # Columns 1-3 are coordinates
  # Column 4 is node color
  # Column 5 is node size
  # Column 6 is node label
  write.table(cbind(coords.cur,
                    color,
                    size,
                    V(g)$name),
              file=paste0(group1, '_', node.size, '_', node.color, '.node'),
              row.names=F, col.names=F, sep='\t', quote=F)

  write.table(as_adj(g, sparse=F),
              file=paste0(group1, '.edge'),
              row.names=F, col.names=F, sep='\t', quote=F)

}
