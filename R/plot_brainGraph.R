#' Plot a brain graph with a specific spatial layout
#'
#' \code{plot.brainGraph} plots a graph in which the spatial layout of the nodes
#' is important. The network itself is plotted over a brain MRI slice from the
#' MNI152 template if \code{mni=TRUE}.
#'
#' With the argument \code{subgraph}, you can specify a simple logical equation
#' for which vertices to show. For example, \emph{'degree > 10'} will plot only
#' vertices with a \emph{degree} greater than 10. Combinations of \emph{AND}
#' (i.e., \code{&}) and \emph{OR} (i.e., \code{|}) are allowed.
#'
#' To remove the subtitle at the bottom, simply specify \code{subt=NULL}.
#'
#' @param x A \code{brainGraph} graph object
#' @param plane Character string indicating which orientation to plot
#'   (default: \code{'axial'})
#' @param hemi Character string indicating which hemisphere to plot (default:
#'   \code{'both'})
#' @param subgraph Character string specifying an equation for vertices to plot
#'   (default: \code{NULL})
#' @param show.legend Logical indicating whether or not to show a legend
#'   (default: \code{FALSE})
#' @param rescale Logical, whether to rescale the coordinates (default:
#'   \code{FALSE})
#' @param asp Numeric constant; the aspect ratio (default: 0)
#' @param main Character string; the main title (default: \code{NULL})
#' @param subt Character string; the subtitle (default: \code{default})
#' @param mni Logical indicating whether or not to plot over a slice of the
#'   brain (default: \code{TRUE})
#' @param ... Other parameters (passed to \code{\link[igraph]{plot.igraph}}).
#'   See \code{\link[igraph]{igraph.plotting}} for details.
#' @export
#' @importFrom oro.nifti image nifti
#'
#' @family Plotting functions
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' plot(g[[1]], hemi='R')
#' plot(g[[1]], subgraph='degree > 10 | btwn.cent > 50')
#' }

plot.brainGraph <- function(x, plane=c('axial', 'sagittal', 'circular'),
                            hemi=c('both', 'L', 'R'),
                            subgraph=NULL, show.legend=FALSE,
                            rescale=FALSE, asp=0, main=NULL, subt='default', mni=TRUE, ...) {
  stopifnot(is.brainGraph(x), 'atlas' %in% graph_attr_names(x))
  lobe <- network <- NULL

  plane <- match.arg(plane)
  hemi <- match.arg(hemi)
  if (plane == 'sagittal') hemi <- 'L'  # Choose other than "both" if argument not provided
  mult <- 100

  if (isTRUE(mni) & plane != 'circular') {
    if (plane == 'axial') {
      slice <- 46
      X <- mni152
      slicemax <- max(X[, , slice])
    } else {
      if (hemi == 'R') {
        X <- mni152
      } else if (hemi == 'L') {
        tmp <- mni152@.Data
        X <- nifti(tmp[rev(seq_len(nrow(tmp))), rev(seq_len(ncol(tmp))), ])
      }
      slice <- 30
      slicemax <- max(X[slice, , ])
    }
    image(X, plot.type='single', plane=plane, z=slice, zlim=c(3500, slicemax))
    par(new=TRUE, mai=c(0, 0, 0, 0), mar=c(0, 0, 0, 0))
  }

  #---------------------------------------------------------
  # Create a subgraph based on user-specified condition
  #---------------------------------------------------------
  if (!is.null(subgraph)) x <- subset_graph(x, subgraph)$g
  if (is.null(x)) return()  # No vertices met criteria

  atlas <- graph_attr(x, 'atlas')
  Nv <- vcount(x)
  adjust <- 0
  #---------------------------------------------------------
  # Handle plotting of individual hemispheres
  #---------------------------------------------------------
  if (hemi != 'both') {
    memb <- which(V(x)$hemi == hemi)
    sg <- induced.subgraph(x, memb)
    sg <- delete_all_attr(sg)
    V(sg)$name <- V(x)$name[memb]
    x <- x %s% sg
    x <- x - vertices(setdiff(seq_len(Nv), memb))

    if (plane == 'sagittal') {
      adjust <- 1
      V(x)$y <- V(x)$z.mni

      if (hemi == 'L') {
        V(x)$x <- -V(x)$y.mni
        xlim.g <- switch(atlas,
                         destrieux=, destrieux.scgm=, dosenbach160=c(-85, 120),
                         hoa112=c(-85, 125),
                         lpba40=c(-75, 120),
                         c(-85, 110))
        ylim.g <- switch(atlas,
                         destrieux=, destrieux.scgm=, dosenbach160=c(-85, 120),
                         hoa112=c(-85, 130),
                         lpba40=c(-80, 100),
                         c(-85, 125))
      } else if (hemi == 'R') {
        V(x)$x <- V(x)$y.mni
        xlim.g <- switch(atlas,
                         hoa112=c(-125, 85),
                         lpba40=c(-125, 75),
                         dosenbach160=c(-120, 85),
                         c(-115, 85))
        ylim.g <- switch(atlas,
                         aal90=, aal116=, aal2.94=, aal2.120=, brainsuite=c(-85, 125),
                         hoa112=c(-85, 130),
                         lpba40=c(-87, 107),
                         c(-85, 120))
      }
    } else if (plane == 'axial') {
      xlim.g <- switch(atlas,
                       dk=, dk.scgm=, dkt=, dkt.scgm=c(-105, 105),
                       hoa112=, lpba40=c(-100, 100),
                       dosenbach160=c(-95, 95),
                       c(-92, 95))
      ylim.g <- switch(atlas,
                       aal90=, aal116=, aal2.94=, aal2.120=c(-110, 70),
                       brainsuite=, dosenbach160=c(-115, 85),
                       destrieux=, destrieux.scgm=c(-120, 83),
                       hoa112=c(-122, 77),
                       lpba40=c(-115, 70),
                       c(-125, 85))
    }
  } else {
    if (plane != 'circular') {
      xlim.g <- switch(atlas,
                       dk=, dk.scgm=, dkt=, dkt.scgm=c(-105, 105),
                       hoa112=, lpba40=c(-100, 100),
                       dosenbach160=c(-95, 95),
                       c(-92, 95))
      ylim.g <- switch(atlas,
                       aal90=, aal116=, aal2.94=, aal2.120=c(-110, 70),
                       brainsuite=, dosenbach160=c(-115, 85),
                       destrieux=, destrieux.scgm=c(-120, 83),
                       hoa112=c(-122, 77),
                       lpba40=c(-115, 70),
                       c(-125, 85))
    }
  }
  if (plane == 'circular') {
    par(bg='black')
    mult <- 1
    layout.x <- rotation(layout.circle(x, order=V(x)$circle.layout), -pi/2 - pi/Nv)
    V(x)$x <- layout.x[, 1]
    V(x)$y <- layout.x[, 2]
    xlim.g <- c(-1.5, 1.2)
    ylim.g <- c(-1, 1.25)
  }

  #---------------------------------------------------------
  # Handle extra arguments in case a subgraph was created
  #---------------------------------------------------------
  fargs <- list(...)
  vcols <- 'lightblue'
  if (hasArg('vertex.color')) {
    if (is.character(fargs$vertex.color) & length(fargs$vertex.color) == 1) {
      stopifnot(fargs$vertex.color %in% vertex_attr_names(x))
      vcols <- vertex_attr(x, fargs$vertex.color)
      if (adjust == 1) {
        medial <- which(abs(V(x)$x.mni) < 20)
        vcols[medial] <- adjustcolor(vcols[medial], alpha.f=0.4)
      }
    } else {
      vcols <- fargs$vertex.color
    }
  }
  ecols <- 'red'
  if (hasArg('edge.color')) {
    if (is.character(fargs$edge.color) & length(fargs$edge.color) == 1) {
      stopifnot(fargs$edge.color %in% edge_attr_names(x))
      ecols <- edge_attr(x, fargs$edge.color)
      if (adjust == 1) {
        medial.es <- as.numeric(E(x)[medial %--% medial])
        ecols[medial.es] <- adjustcolor(ecols[medial.es], alpha.f=0.4)
      }
    } else {
      ecols <- fargs$edge.color
    }
  }

  # Vertex sizes
  #-------------------------------------
  vsize <- mult * 10
  if (hasArg('vertex.size')) {
    if (is.character(fargs$vertex.size)) {
      stopifnot(fargs$vertex.size %in% vertex_attr_names(x))
      vsize <- mult * vertex_attr(x, fargs$vertex.size)
      if (any(!is.finite(vsize))) x <- delete.vertices(x, which(!is.finite(vsize)))
      vsize <- mult * vec.transform(vsize, as.numeric(min(vsize) > 0), 15)
    } else {
      vsize <- mult * fargs$vertex.size
    }
  }
  ewidth <- 1.5
  if (hasArg('edge.width')) {
    if (is.character(fargs$edge.width)) {
      stopifnot(fargs$edge.width %in% edge_attr_names(x))
      ewidth <- edge_attr(x, fargs$edge.width)
    } else {
      ewidth <- fargs$edge.width
    }
  }
  # Vertex labels
  #-------------------------------------
  vlabel <- V(x)$name
  vlabel.dist <- ifelse(vsize >= mult * 10, 0, 5)
  vlabel.col <- ifelse(vcols %in% c('red', 'blue', 'magenta'), 'white', 'blue')
  vlabel.font <- 2
  if (hasArg('vertex.label')) {
    if (is.na(fargs$vertex.label)) {
      vlabel <- vlabel.cex <- vlabel.dist <- vlabel.col <- vlabel.font <- NA
    } else {
      vlabel <- vertex_attr(x, fargs$vertex.label)
    }
  }
  if (hasArg('vertex.label.cex')) {
    vlabel.cex <- fargs$vertex.label.cex
  } else {
    vlabel.cex <- 1
  }

  NextMethod(generic='plot', object=x, asp=asp, rescale=rescale, xlim=xlim.g, ylim=ylim.g,
       vertex.color=vcols, edge.color=ecols, vertex.size=vsize,
       edge.width=ewidth, vertex.label=vlabel, vertex.label.cex=vlabel.cex,
       vertex.label.dist=vlabel.dist, vertex.label.font=vlabel.font,
       vertex.label.color=vlabel.col, ...)

  if (!is.null(subt)) {
    if (subt == 'default') {
      Ne <- ecount(x)
      g.density <- round(graph.density(x), digits=3)
      par(new=TRUE, mar=c(5, 0, 3, 0)+0.1)
      subt <- paste('# vertices: ', Nv, '# edges: ', Ne, '\n',
                   'Density: ', g.density)
    }
  }
  if (hasArg('cex.main')) {
    cex.main <- fargs$cex.main
  } else {
    cex.main <- 2.5
  }
  title(main=main, sub=subt, col.main='white', col.sub='white', cex.main=cex.main)

  if (isTRUE(show.legend)) {
    if (hasArg('vertex.color')) {
      atlas.dt <- get(atlas)
      if (fargs$vertex.color == 'color.lobe') {
        lobes.g <- unique(V(x)$lobe)
        classnames <- intersect(levels(atlas.dt$lobe), lobes.g)
        total <- unname(atlas.dt[, table(lobe)][classnames])
        classnames <- paste0(classnames, ': ', table(V(x)$lobe)[classnames],
                             ' / ', total)
        cols <- group.cols[which(levels(atlas.dt$lobe) %in% lobes.g)]
        cex <- vlabel.cex
      } else if (fargs$vertex.color == 'color.class') {
        classnames <- levels(atlas.dt[, class])
        cols <- c('red', 'green', 'blue')
        cex <- 1.5
      } else if (fargs$vertex.color == 'color.rich') {
        classnames <- c('Rich-club', 'Feeder', 'Local')
        cols <- c('red', 'orange', 'green')
        cex <- 1.5
      } else if (fargs$vertex.color == 'color.network') {
        networks.g <- unique(V(x)$network)
        classnames <- intersect(levels(atlas.dt$network), networks.g)
        total <- unname(atlas.dt[, table(network)][classnames])
        classnames <- paste0(classnames, ': ', table(V(x)$network)[classnames],
                             ' / ', total)
        cols <- group.cols[which(levels(atlas.dt$network) %in% networks.g)]
        cex <- vlabel.cex
      } else if (fargs$vertex.color == 'color.comm') {
        memb <- V(x)$comm
        tab <- table(memb)
        tab <- tab[tab >= 3]
        group.nums <- as.integer(names(tab))
        classnames <- paste('Community', group.nums)
        cols <- group.cols[group.nums]
        cex <- vlabel.cex
      }
      l.inset <- ifelse(plane == 'circular', c(-0.1, 0), c(-0.01, 0))
      legend('topleft',
             classnames,
             fill=cols,
             bg='black',
             text.col='white',
             cex=cex,
             inset=l.inset)
    }
  }
}

#' @inheritParams plot.brainGraph
#' @export
#' @rdname plot.brainGraph

plot_brainGraph <- function(x, plane=c('axial', 'sagittal', 'circular'),
                            hemi=c('both', 'L', 'R'),
                            subgraph=NULL, show.legend=FALSE,
                            rescale=FALSE, asp=0, main=NULL, subt='default', mni=TRUE, ...) {
  plot.brainGraph(x, plane, hemi, subgraph, show.legend, rescale, asp, main, subt, mni, ...)
}

#' Plot a graph with results from the network-based statistic
#'
#' This is a convenience function for plotting a graph based on results from
#' \code{\link{NBS}}. There are several default arguments that are set:
#' vertex/edge colors will correspond to connected component membership, and
#' only those vertices in which \code{V(g)$p.nbs > 1 - alpha} will be shown.
#' Finally, vertex names will be omitted.
#'
#' @param x A \code{brainGraph_NBS} graph object (from
#'   \code{\link{make_nbs_brainGraph}})
#' @param alpha Numeric; the significance level (default: 0.05)
#' @param subgraph Character string specifying the condition for subsetting the
#'   graph. By default, it will show only the vertices which are members of
#'   components determined to be significant based on \code{alpha}.
#' @param vertex.label Character vector of the vertex labels to be displayed.
#'   Default behavior is to omit them.
#' @param vertex.color Character string specifying the vertex attribute to color
#'   the vertices by (default: \code{color.comp}, which groups vertices by
#'   connected component)
#' @param edge.color Character string specifying the edge attribute to color
#'   the edges by (default: \code{color.comp}, which groups edges by connected
#'   component)
#' @param cex.main Numeric; the scaling factor for text size; see
#'   \code{\link[graphics]{par}} (default: 2)
#' @param ... Other arguments passed to \code{\link{plot.brainGraph}}
#' @inheritParams plot.brainGraph
#' @export
#' @method plot brainGraph_NBS
#' @family Plotting functions

plot.brainGraph_NBS <- function(x, alpha=0.05, subgraph=paste('p.nbs >', 1 - alpha),
                                vertex.label=NA,
                                vertex.color='color.comp', edge.color='color.comp', subt=NULL,
                                main=paste0('\n\nNBS: ', x$name), cex.main=2, ...) {
  NextMethod(generic='plot', object=x, subgraph=subgraph,
             vertex.label=vertex.label, vertex.color=vertex.color, edge.color=edge.color,
             subt=subt, main=main, cex.main=cex.main, ...)
}

#' Plot a graph with results from brainGraph_GLM
#'
#' This is a convenience function for plotting a graph based on results from
#' \code{\link{brainGraph_GLM}}. There are a few argument defaults: to plot only those
#' vertices for which \eqn{p < \alpha}; a plot title with the outcome
#' measure and contrast name, and to omit the plot subtitle.
#'
#' @param x A \code{brainGraph_GLM} graph object (from
#'   \code{\link{make_glm_brainGraph}})
#' @param p.sig Character string indicating which p-value to use for determining
#'   significance (default: \code{p})
#' @param cex.main Numeric indicating the scaling for plot title size (see
#'   \code{\link[graphics]{par}}.
#' @inheritParams plot.brainGraph
#' @export
#' @method plot brainGraph_GLM
#' @family Plotting functions

plot.brainGraph_GLM <- function(x, p.sig=c('p', 'p.fdr', 'p.perm'),
                                subgraph=NULL,
                                main=paste0('\n\n', x$outcome, ': ', x$name),
                                subt=NULL, cex.main=2, ...) {
  p.sig <- match.arg(p.sig)
  if (is.null(subgraph)) subgraph <- paste0(p.sig, ' > 1 -', x$alpha)
  NextMethod(generic='plot', object=x, subgraph=subgraph,
             main=main, subt=subt, cex.main=cex.main, ...)
}

#' Plot a graph with results from MTPC
#'
#' This is a convenience function for plotting a graph based on results from
#' \code{\link{mtpc}}. There are a few argument defaults: to plot only those
#' vertices for which \eqn{A_{mtpc} > A_{crit}}; a plot title with the outcome
#' measure and contrast name, and to omit the plot subtitle.
#'
#' @param x A \code{brainGraph_mtpc} graph object (from
#'   \code{\link{make_glm_brainGraph}})
#' @param cex.main Numeric indicating the scaling for plot title size (see
#'   \code{\link[graphics]{par}}.
#' @inheritParams plot.brainGraph
#' @export
#' @method plot brainGraph_mtpc
#' @family Plotting functions

plot.brainGraph_mtpc <- function(x, subgraph='sig == 1',
                                 main=paste0('\n\n', x$outcome, ': ', x$name),
                                 subt=NULL, cex.main=2, ...) {
  NextMethod(generic='plot', object=x, subgraph=subgraph,
             main=main, subt=subt, cex.main=cex.main, ...)
}

#' Plot a graph with results from a mediation analysis
#'
#' @param x A \code{brainGraph_mediate} graph object (from
#'   \code{\link{make_mediate_brainGraph}})
#' @param cex.main Numeric indicating the scaling for plot title size (see
#'   \code{\link[graphics]{par}}.
#' @inheritParams plot.brainGraph
#' @export
#' @method plot brainGraph_mediate
#' @family Plotting functions

plot.brainGraph_mediate <- function(x, subgraph='p.acme > 0.95',
                                    main=sprintf('\n\n\nEffect of "%s" on\n"%s"\nmediated by "%s"', x$treat, x$outcome, x$mediator),
                                    subt=NULL, cex.main=1, ...) {
  NextMethod(generic='plot', object=x, subgraph=subgraph, main=main, subt=subt, cex.main=cex.main, ...)
}
