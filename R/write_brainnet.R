#' Write files to be used for visualization with BrainNet Viewer
#'
#' This function will write the "node" and "edge" files necessary for
#' visualization with the BrainNet Viewer software.
#'
#' @param g A graph
#' @param node.size Character string indicating what size the nodes should be;
#' can be any vertex-level attribute (see \link{list.vertex.attributes})
#' @param node.color Character string indicating whether to color the nodes or
#' not; can be "none" or "community"
#' @export
#'

write.brainnet <- function(g, node.color=c('none', 'community'), node.size='degree') {
  coords.cur <- cbind(V(g)$x, V(g)$y, V(g)$z)

  if (node.color == 'none') {
    color <- rep(1, vcount(g))
  } else if (node.color == 'community') {
    color <- V(g)$comm
  }

  size <- get.vertex.attribute(g, node.size)

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
              file=paste0(out.dir, group1, '_', node.size, '.node'),
              row.names=F, col.names=F, sep='\t', quote=F)

  write.table(get.adjacency(g, sparse=F),
              file=paste0(out.dir, group1, '.edge'),
              row.names=F, col.names=F, sep='\t', quote=F)

}
