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
#' @param show.legend Logical indicating whether or not to show a legend
#'   (default: FALSE)
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
                            subgraph=NULL, show.legend=FALSE,
                            rescale=FALSE, ylim=c(-1.5, 1.5),
                            asp=0, main=NULL, ...) {
  stopifnot(is_igraph(g))
  lobe <- NULL

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

  plane <- match.arg(plane)
  hemi <- match.arg(hemi)
  Nv <- vcount(g)
  adjust <- 0
  if (hemi != 'both') {
    memb <- which(V(g)$hemi == hemi)
    sg <- induced.subgraph(g, memb)
    for (att in graph_attr_names(sg)) sg <- delete_graph_attr(sg, att)
    for (att in vertex_attr_names(sg)) sg <- delete_vertex_attr(sg, att)
    for (att in edge_attr_names(sg)) sg <- delete_edge_attr(sg, att)
    V(sg)$name <- V(g)$name[memb]
    g <- g %s% sg
    g <- g - vertices(setdiff(seq_len(Nv), memb))

    if (plane == 'sagittal') {
      adjust <- 1
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
  } else {
    mult <- 1
    xlim.g <- c(-1, 1)
    ylim.g <- ylim
  }

  #---------------------------------------------------------
  # Handle extra arguments in case a subgraph was created
  #---------------------------------------------------------
  fargs <- list(...)
  if (hasArg('vertex.color')) {
    if (is.character(fargs$vertex.color) & length(fargs$vertex.color) == 1) {
      if (fargs$vertex.color %in% vertex_attr_names(g)) {
        vcols <- vertex_attr(g, fargs$vertex.color)
        if (adjust == 1) {
          medial <- which(abs(V(g)$x.mni) < 20)
          vcols[medial] <- adjustcolor(vcols[medial], alpha.f=0.4)
        }
      }
    } else {
      vcols <- fargs$vertex.color
    }
  } else {
    vcols <- 'lightblue'
  }
  if (hasArg('edge.color') && length(fargs$edge.color) == 1) {
    ecols <- edge_attr(g, fargs$edge.color)
    if (adjust == 1) {
      medial.es <- as.numeric(E(g)[medial %--% medial])
      ecols[medial.es] <- adjustcolor(ecols[medial.es], alpha.f=0.4)
    }
  } else {
    ecols <- 'red'
  }

  # Vertex sizes
  #-------------------------------------
  if (hasArg('vertex.size')) {
    if (is.character(fargs$vertex.size)) {
      vsize <- mult * vertex_attr(g, fargs$vertex.size)
      if (any(!is.finite(vsize))) g <- delete.vertices(g, which(!is.finite(vsize)))
      vsize <- vec.transform(vsize, 0, 15)
    } else {
      vsize <- mult * fargs$vertex.size
    }
  } else {
    vsize <- mult * 15
  }
  if (hasArg('edge.width')) {
    if (is.character(fargs$edge.width)) {
      ewidth <- edge_attr(g, fargs$edge.width)
    } else {
      ewidth <- fargs$edge.width
    }
  } else {
    ewidth <- 1.5
  }
  # Vertex labels
  #-------------------------------------
  if (hasArg('vertex.label')) {
    if (is.na(fargs$vertex.label)) {
      vlabel <- vlabel.cex <- vlabel.dist <- vlabel.col <- vlabel.font <- NA
    } else {
      vlabel <- V(g)$name
      vlabel.cex <- 0.75
      vlabel.dist <- ifelse(vsize >= 10, 0, 0.75)
      vlabel.col <- ifelse(vcols %in% c('red', 'blue'), 'white', 'blue')
      vlabel.font <- 2
    }
  }

  plot(g, asp=asp, rescale=rescale, xlim=xlim.g, ylim=ylim.g,
       vertex.color=vcols, edge.color=ecols, vertex.size=vsize,
       edge.width=ewidth, vertex.label=vlabel, vertex.label.cex=vlabel.cex,
       vertex.label.dist=vlabel.dist, vertex.label.font=vlabel.font,
       vertex.label.color=vlabel.col, ...)

  Nv <- vcount(g)
  Ne <- ecount(g)
  g.density <- round(graph.density(g), digits=3)
  par(new=T, mar=c(5, 0, 3, 0)+0.1)
  subt <- paste('# vertices: ', Nv, '# edges: ', Ne, '\n',
                'Density: ', g.density)
  title(main=main, sub=subt, col.main='white', col.sub='white', cex.main=2)

  if (isTRUE(show.legend)) {
    if (hasArg('vertex.color')) {
      atlas.dt <- eval(parse(text=g$atlas))
      if (fargs$vertex.color == 'color.lobe') {
        lobes <- sort(unique(V(g)$lobe))
        classnames <- levels(atlas.dt[, lobe])[lobes]
        total <- unname(atlas.dt[, table(lobe)])[lobes]
        classnames <- paste0(classnames, ': ', table(V(g)$lobe),
                             ' / ', total)
        cols <- unique(V(g)$color.lobe[order(V(g)$lobe)])
        cex <- 1
      } else if (fargs$vertex.color == 'color.class') {
        classnames <- levels(atlas.dt[, class])
        cols <- c('red', 'green', 'blue')
      } else if (fargs$vertex.color == 'color.rich') {
        classnames <- c('Rich-club', 'Feeder', 'Local')
        cols <- c('red', 'orange', 'green')
        cex <- 1.5
      }
      legend('topleft',
             classnames,
             fill=cols,
             bg='black',
             text.col='white',
             cex=cex)
    }
  }
}
