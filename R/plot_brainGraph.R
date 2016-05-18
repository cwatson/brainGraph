#' Plot a brain graph with a specific spatial layout
#'
#' This function plots a graph when the spatial layout of the nodes is important
#' (e.g. in the brain). The function \code{\link{set.brainGraph.attributes}}
#' needs to be run on the graph, and a valid set of coordinates provided for the
#' vertices. Most of the parameters valid here can be seen in
#' \code{\link{igraph.plotting}}.
#'
#' With the argument \code{subgraph}, you can specify a simple logical equation
#' for which vertices to show. For example, \emph{'degree > 10'} will plot only
#' vertices with a \emph{degree} greater than 10. Combinations of \emph{AND}
#' (i.e., \code{&}) and \emph{OR} (i.e., \code{|}) are allowed.
#'
#' @param g An \code{igraph} graph object
#' @param plane A character string indicating which orientation to plot
#'   (default: 'axial')
#' @param hemi A character string indicating which hemisphere to plot (default:
#'   'both')
#' @param subgraph A character string specifying an equation for deleting
#'   vertices (default: NULL)
#' @param rescale A logical, whether to rescale the coordinates (default: FALSE)
#' @param ylim A vector giving limits for the vertical axis
#'   (default: c(-1.5, 1.5))
#' @param asp A numeric constant for the aspect ratio (default: 0)
#' @param main Character string for the main title (default: NULL)
#' @param ... Other parameters (passed to \code{\link{plot}}).
#' @export
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' plot_brainGraph(g[[1]], hemi='R')
#' plot_brainGraph(g[[1]], subgraph='degree > 10 | btwn.cent > 50')
#' }

plot_brainGraph <- function(g, plane=c('axial', 'sagittal'), hemi=c('both', 'L', 'R'),
                            subgraph=NULL,
                            rescale=FALSE, ylim=c(-1.5, 1.5),
                            asp=0, main=NULL, ...) {
  stopifnot(is_igraph(g))

  plane <- match.arg(plane)

  # Create a subgraph based on user-specified condition
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
                                    paste0('V(g)$', subs[x], splits[x])),
                             collapse='')
      } else {
        cond.string <- paste0('V(g)$', subs)
      }
      cond <- eval(parse(text=cond.string))
      cond <- setdiff(seq_len(vcount(g)), which(cond))
      g <- delete.vertices(g, cond)
    } else {
      stop('"subgraph" must be a valid character string')
    }
  }

  hemi <- match.arg(hemi)
  if (hemi != 'both') {
    Nv <- vcount(g)
    if (plane == 'sagittal') {
      mult <- 100
      ylim.g <- c(-85, 125)
      V(g)$y <- V(g)$z.mni
      if (hemi == 'L') {
        xlim.g <- c(-85, 110)
        V(g)$x <- -V(g)$y.mni
      } else if (hemi == 'R') {
        V(g)$x <- V(g)$y.mni
        xlim.g <- c(-125, 85)
      }
    } else {
      mult <- 1
      ylim.g <- ylim
      xlim.g <- c(-1, 1)
    }
    memb <- which(V(g)$hemi == hemi)
    sg <- induced.subgraph(g, memb)
    for (att in graph_attr_names(sg)) sg <- delete_graph_attr(sg, att)
    for (att in vertex_attr_names(sg)) sg <- delete_vertex_attr(sg, att)
    for (att in edge_attr_names(sg)) sg <- delete_edge_attr(sg, att)
    V(sg)$name <- V(g)$name[memb]
    g <- g %s% sg
    g <- g - vertices(setdiff(seq_len(Nv), memb))

  } else {
    mult <- 1
    xlim.g <- c(-1, 1)
    ylim.g <- ylim
  }

  # Handle extra arguments in case a subgraph was created
  fargs <- list(...)
  if (hasArg('vertex.color')) {
    if (fargs$vertex.color %in% vertex_attr_names(g)) {
      vcols <- vertex_attr(g, fargs$vertex.color)
    } else {
      vcols <- fargs$vertex.color
    }
  } else {
    vcols <- 'lightblue'
  }
  if (hasArg('edge.color')) {
    ecols <- edge_attr(g, fargs$edge.color)
  } else {
    ecols <- 'red'
  }

  # Vertex sizes
  #-------------------------------------
  if (hasArg('vertex.size')) {
    if (is.character(fargs$vertex.size)) {
      vsize <- mult * vertex_attr(g, fargs$vertex.size)
    } else {
      vsize <- mult * fargs$vertex.size
    }
  } else {
    vsize <- mult * 15
  }
  if (hasArg('edge.width')) {
    if (is.character(fargs$edge.width)) {
      ewidth <- edge_attr(g, fargs$edge.width)
      ewidth <- vec.transform(ewidth, min(ewidth), 20)
    } else {
      ewidth <- fargs$edge.width
    }
  } else {
    ewidth <- 1.5
  }

  plot(g, asp=asp, rescale=rescale, xlim=xlim.g, ylim=ylim.g,
       vertex.color=vcols, edge.color=ecols, vertex.size=vsize,
       edge.width=ewidth, ...)

  Nv <- vcount(g)
  Ne <- ecount(g)
  g.density <- round(graph.density(g), digits=3)
  par(new=T, mar=c(5, 0, 3, 0)+0.1)
  subt <- paste('# vertices: ', Nv, '# edges: ', Ne, '\n',
                'Density: ', g.density)
  title(main=main, sub=subt, col.main='white', col.sub='white')
}
