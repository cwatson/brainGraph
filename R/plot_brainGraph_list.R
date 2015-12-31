#' Write PNG files for a list of graphs
#'
#' This function takes a list of \code{igraph} graph objects and plots them over
#' an axial slice of the brain. A \emph{png} file is written for each element of
#' the list, which can be joined as a \emph{gif} or converted to video using a
#' tool outside of R.
#'
#' You can choose to highlight edge differences between subsequent list
#' elements, and whether to color vertices by \emph{lobe}, \emph{community}
#' membership, or \emph{lightblue} (the default) (or a color of your choosing).
#' By default, the vertex sizes are equal to vertex degree, and max out at 20;
#' however, you may choose other values. Finally, you may choose to plot only a
#' subgraph of vertices based on some criteria (see examples).
#'
#' This function may be particularly useful if the graph list contains graphs of
#' a single subject group at incremental densities, or if the graph list
#' contains graphs of each subject in a group.
#'
#' @param g.list A list of \code{igraph} graph objects
#' @param fname.base A character string specifying the base of the filename for
#'   \emph{png} output
#' @param diffs A logical, indicating whether or not to highlight edge
#'   differences (default: FALSE)
#' @param subgraph A character string specifying an equation for deleting
#'   vertices (default: NULL)
#' @param ... Other parameters (passed to \code{\link{plot_brainGraph}})
#' @export
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' plot_brainGraph_list(g[[1]], 'g1', subgraph='hemi == "R"')
#' plot_brainGraph_list(g[[1]], 'g1', subgraph='degree > 10 | btwn.cent > 50')
#' }

plot_brainGraph_list <- function(g.list, fname.base, diffs=FALSE,
                                 subgraph=NULL, ...) {

  for (i in seq_along(g.list)) {
    png(filename=sprintf('%s_%03d.png', fname.base, i))

    if (!is.null(subgraph)) {
      if (nchar(subgraph) > 0) {
        subs <- strsplit(subgraph, split='\\&|\\|')[[1]]
        if (length(subs) > 1) {
          nchars <- cumsum(sapply(subs, nchar))
          splits <- sapply(seq_along(subs), function(x)
                           substr(subgraph, start=nchars[x]+x, stop=nchars[x]+x))
          subs <- gsub('^\\s+|\\s+$', '', subs) # Remove unnecessary whitespace
          # In case there is a mix of '&' and '|'
          cond.string <- paste(sapply(seq_along(subs), function(x)
                                      paste0('V(g.list[[', i, ']])$', subs[x], splits[x])),
                               collapse='')
        } else {
          cond.string <- paste0('V(g.list[[', i, ']])$', subs)
        }
        cond <- eval(parse(text=cond.string))
        cond <- setdiff(seq_len(vcount(g.list[[i]])), which(cond))
        g.list[[i]] <- delete.vertices(g.list[[i]], cond)
      } else {
        stop(sprintf('%s must be a valid character string', subgraph))
      }
    }

    fargs <- list(...)
    # Choose different vertex colors
    if (hasArg('vertex.color')) {
      if (fargs$vertex.color == 'color.lobe' || fargs$vertex.color == 'color.comm') {
        vcols <- vertex_attr(g.list[[i]], fargs$vertex.color)
        ecols <- edge_attr(g.list[[i]], fargs$vertex.color)
      } else {
        vcols <- fargs$vertex.color
        ecols <- 'red'
      }
    } else {
      vcols <- 'lightblue'
      ecols <- 'red'
    }

    # Choose different vertex sizes
    if (hasArg('vertex.size')) {
      if (is.character(fargs$vertex.size)) {
        vsize <- vertex_attr(g.list[[i]], fargs$vertex.size)
        vsize <- vec.transform(vsize, min(vsize), 20)
      } else {
        vsize <- fargs$vertex.size
      }
    } else {
      vsize <- pmin(V(g.list[[i]])$degree, 20)
    }
    plot_brainGraph_mni('axial')
    plot_brainGraph(g.list[[i]],
                    vertex.size=vsize,
                    vertex.color=vcols, edge.color=ecols,
                    main=g.list[[i]]$Group, ...)

    if (isTRUE(diffs)) {
      if (i > 1) {
        g.diff <- graph.difference(g.list[[i]], g.list[[i-1]])
        if (hasArg('vertex.color')) {
          if (fargs$vertex.color == 'color.lobe' || fargs$vertex.color == 'color.comm') {
            ecols <- edge_attr(g.diff, fargs$vertex.color)
          } else {
            ecols <- 'red'
          }
        } else {
          ecols <- 'deeppink'
        }
        ewidth <- 5
        plot_brainGraph(g.diff, add=T, vertex.label=NA, vertex.shape='none',
                        edge.width=ewidth, edge.color=ecols)
      }
    }
    dev.off()
  }
}
