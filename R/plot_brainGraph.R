#' Plot a brain graph with a specific spatial layout
#'
#' \code{plot.brainGraph} plots a graph in which the spatial layout of the nodes
#' is important. The network itself is plotted over a brain MRI slice from the
#' MNI152 template by default (when \code{mni=TRUE}).
#'
#' @section Selecting specific vertices to display:
#' With the argument \code{subgraph}, you can supply a simple logical expression
#' specifying which vertices to show. For example, \emph{'degree > 10'} will
#' plot only vertices with a \emph{degree} greater than 10. Combinations of
#' \emph{AND} (i.e., \code{&}) and \emph{OR} (i.e., \code{|}) are allowed. This
#' requires that any vertex attribute in the expression must be present in the
#' graph; e.g., \code{V(g)$degree} must exist.
#'
#' @section Title, subtitle, and label:
#' By default, a \emph{title} (i.e., text displayed at the top of the figure) is
#' not included. You can include one by passing a character string to
#' \code{main}, and control the size with \code{cex.main}. A \emph{subtitle}
#' (i.e., text at the bottom), is included by default and displays the number of
#' vertices and edges along with the graph density. To exclude this, specify
#' \code{subtitle=NULL}. A \dQuote{label} can be included in one corner of the
#' figure (for publications). For example, you can choose \code{label='A.'} or
#' \code{label='a)'}. Arguments controlling the location and appearance can be
#' changed, but the default values are optimal for bottom-left placement. See
#' \code{\link[graphics]{mtext}} for more details. The label-specific arguments
#' are:
#' \describe{
#'   \item{side}{The location. \code{1} is for bottom placement.}
#'   \item{line}{If \code{side=1} (bottom), a negative number places the text
#'     \emph{above} the bottom of the figure; a higher number could result in
#'     the bottom part of the text to be missing. This can differ if
#'     \code{plane='circular'}, in which case you may want to specify a positive
#'     number.}
#'   \item{adj}{Seems to be the percentage away from the margin. So, for
#'     example, \code{adj=0.1} would place the text closer to the center than
#'     the default value, and \code{adj=0.5} places it in the center.}
#'   \item{cex}{The degree of \dQuote{character expansion}. A value of 1 would
#'     not increase the text size.}
#'   \item{col}{The text color.}
#'   \item{font}{The font type. The default \code{font=2} is bold face. See
#'     \code{\link[graphics]{par}} for details.}
#' }
#'
#' @param x A \code{brainGraph} graph object
#' @param plane Character string indicating which orientation to plot.
#'   Default: \code{'axial'}
#' @param hemi Character string indicating which hemisphere to plot. Default:
#'   \code{'both'}
#' @param mni Logical indicating whether or not to plot over a slice of the
#'   brain. Default: \code{TRUE}
#' @param subgraph Character string specifying a logical condition for vertices
#'   to plot. Default: \code{NULL}
#' @param main Character string; the main title. Default: \code{NULL}
#' @param subtitle Character string; the subtitle. Default: \code{'default'}
#' @param label Character string specifying text to display in one corner of the
#'   plot (e.g., \code{'A.'}). Default: NULL
#' @param side Label placement. Default: \code{1} (bottom)
#' @param line Which margin line to place the text.
#' @param adj If \code{side=1}, a value closer to 0 places the text closer to
#'   the left margin. Default: \code{0.025}
#' @param cex Amount of character expansion of the label text. Default:
#'   \code{2.5}
#' @param col Label font color. Default: \code{'white'}
#' @param font Integer specifying the font type. Default: \code{2} (bold face)
#' @param show.legend Logical indicating whether or not to show a legend.
#'   Default: \code{FALSE}
#' @param rescale Logical, whether to rescale the coordinates. Default:
#'   \code{FALSE}
#' @param asp Numeric constant; the aspect ratio. Default: 0
#' @param ... Other parameters (passed to \code{\link[igraph]{plot.igraph}}).
#'   See \code{\link[igraph]{plot.common}} for details.
#' @export
#' @importFrom graphics legend mtext par title
#' @importFrom methods hasArg
#'
#' @family Plotting functions
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' plot(g[[1]], hemi='R')
#' plot(g[[1]], subgraph='degree > 10 | btwn.cent > 50')
#'
#' ## Place label in upper-left
#' plot(g.ex, label='A)', side=3, line=-2.5)
#' }

plot.brainGraph <- function(x, plane=c('axial', 'sagittal', 'circular'),
                            hemi=c('both', 'L', 'R'), mni=TRUE, subgraph=NULL,
                            main=NULL, subtitle='default', label=NULL, side=1,
                            line=-2, adj=0.025, cex=2.5, col='white', font=2,
                            show.legend=FALSE, rescale=FALSE, asp=0, ...) {
  stopifnot(is.brainGraph(x))

  plane <- match.arg(plane)
  hemi <- match.arg(hemi)
  if (plane == 'sagittal' && hemi == 'both') hemi <- 'L'
  if (isTRUE(mni) && plane != 'circular') plot_mni(plane, hemi)

  mult <- 100
  adjust <- 0L
  atlas <- graph_attr(x, 'atlas')
  Nv <- vcount(x)
  keep <- seq_len(Nv)
  if (plane == 'circular') {
    layout.circ <- rotation(layout.circle(x, order=V(x)$circle.layout), -pi/2 - pi/Nv)
  }
  #---------------------------------------------------------
  # Set figure window limits and handle plotting of individual hemispheres
  #---------------------------------------------------------
  if (hemi != 'both') subgraph <- paste0(subgraph, " & hemi == '", hemi, "'")
  if (!is.null(subgraph)) {
    if (substr(subgraph, 1L, 3L) == ' & ') subgraph <- substr(subgraph, 4L, nchar(subgraph))
    sg <- subset_graph(x, subgraph)
    x <- sg$g
    keep <- sg$inds
    Nv <- vcount(x)
  }

  if (plane == 'sagittal') {
    adjust <- 1L
    V(x)$x <- -V(x)$y.mni
    V(x)$y <- V(x)$z.mni
    xlim.g <- switch(atlas,
                     destrieux=, destrieux.scgm=, dosenbach160=c(-85, 120),
                     hcp_mmp1.0=c(-90, 125), hoa112=c(-85, 125), lpba40=c(-75, 120),
                     brainnetome=, gordon333=, power264=c(-90, 120),
                     c(-85, 110))
    ylim.g <- switch(atlas,
                     destrieux=, destrieux.scgm=, dosenbach160=c(-85, 120),
                     hoa112=c(-85, 130), lpba40=c(-80, 100), c(-85, 125))

    if (hemi == 'R') {
      V(x)$x <- -V(x)$x
      xlim.g <- -rev(xlim.g)
    }
  } else if (plane == 'axial') {
    xlim.g <- switch(atlas,
                     dk=, dk.scgm=, dkt=, dkt.scgm=c(-105, 105),
                     brainnetome=, hcp_mmp1.0=, hoa112=, lpba40=, power264=c(-100, 100),
                     dosenbach160=c(-95, 95), c(-92, 95))
    ylim.g <- switch(atlas,
                     aal90=, aal116=, aal2.94=, aal2.120=c(-110, 70),
                     brainsuite=, dosenbach160=c(-115, 85),
                     destrieux=, destrieux.scgm=, hcp_mmp1.0=c(-120, 83),
                     hoa112=c(-122, 77), lpba40=c(-115, 70), c(-125, 85))
  } else if (plane == 'circular') {
    par(new=TRUE, bg='black', mar=c(5, 0, 2, 0)+0.1)
    mult <- 1
    V(x)$x <- layout.circ[keep, 1L]
    V(x)$y <- layout.circ[keep, 2L]
    xlim.g <- c(-1.5, 1.5)
    ylim.g <- c(-1.00, 1.5)
  }

  #---------------------------------------------------------
  # Handle extra arguments in case a subgraph was created
  #---------------------------------------------------------
  # Default values
  #-------------------------------------
  medial <- which(abs(V(x)$x.mni) < 20)
  vcols <- 'lightblue'
  ecols <- 'red'
  vsize <- mult * 10
  ewidth <- 1.5
  vlabel <- V(x)$name
  fargs <- list(...)

  # Vertex and edge colors
  #-------------------------------------
  if (hasArg('vertex.color')) {
    vcols <- fargs$vertex.color
    if (is.character(vcols) && length(vcols) == 1L) {
      stopifnot(vcols %in% vertex_attr_names(x))
      vcols <- vertex_attr(x, vcols)
      if (adjust == 1L) vcols[medial] <- adjustcolor(vcols[medial], alpha.f=0.4)
    }
  }
  if (hasArg('edge.color')) {
    ecols <- fargs$edge.color
    if (is.character(ecols) && length(ecols) == 1L) {
      stopifnot(ecols %in% edge_attr_names(x))
      ecols <- edge_attr(x, ecols)
      if (adjust == 1L) {
        medial.es <- as.numeric(E(x)[medial %--% medial])
        ecols[medial.es] <- adjustcolor(ecols[medial.es], alpha.f=0.4)
      }
    }
  }

  # Vertex sizes & edge widths
  #-------------------------------------
  if (hasArg('vertex.size')) {
    if (is.character(fargs$vertex.size)) {
      stopifnot(fargs$vertex.size %in% vertex_attr_names(x))
      vsize <- mult * vertex_attr(x, fargs$vertex.size)
      if (any(!is.finite(vsize))) x <- delete_vertices(x, which(!is.finite(vsize)))
      vsize <- mult * vec.transform(vsize, as.numeric(min(vsize) > 0), 15)
    } else {
      vsize <- mult * fargs$vertex.size
    }
  }
  if (hasArg('edge.width')) {
    ewidth <- fargs$edge.width
    if (is.character(ewidth)) {
      stopifnot(ewidth %in% edge_attr_names(x))
      ewidth <- edge_attr(x, ewidth)
    }
  }
  # Vertex labels
  #-------------------------------------
  vlabel.dist <- ifelse(vsize >= mult * 10, 0, 5)
  vlabel.col <- ifelse(vcols %in% c('red', 'blue', 'magenta'), 'white', 'blue')
  vlabel.font <- 2
  if (hasArg('vertex.label')) {
    vlabel <- if (!is.na(fargs$vertex.label)) vertex_attr(x, fargs$vertex.label) else NA
  }
  vlabel.cex <- if (hasArg('vertex.label.cex')) fargs$vertex.label.cex else 1

  NextMethod(generic='plot', object=x, asp=asp, rescale=rescale, xlim=xlim.g, ylim=ylim.g,
       vertex.color=vcols, edge.color=ecols, vertex.size=vsize,
       edge.width=ewidth, vertex.label=vlabel, vertex.label.cex=vlabel.cex,
       vertex.label.dist=vlabel.dist, vertex.label.font=vlabel.font,
       vertex.label.color=vlabel.col, ...)

  if (!is.null(label)) mtext(label, side=side, line=line, adj=adj, cex=cex, col=col, font=font)

  if (!is.null(subtitle) && subtitle == 'default') {
    Ne <- ecount(x)
    g.density <- round(graph.density(x), digits=3L)
    subtitle <- paste('# vertices: ', Nv, '\t# edges: ', Ne, '\n', 'Density: ', g.density)
  }
  par(new=TRUE, mar=c(5, 0, 2, 0)+0.1)
  cex.main <- if (hasArg('cex.main')) fargs$cex.main else 2.5
  title(main=main, sub=subtitle, col.main='white', col.sub='white', cex.main=cex.main)

  if (isTRUE(show.legend) && hasArg('vertex.color')) {
    atlas.dt <- get(atlas)
    vattr <- strsplit(fargs$vertex.color, '.', fixed=TRUE)[[1L]][2L]
    if (vattr %in% c('lobe', 'network', 'class', 'area', 'Yeo_7network', 'Yeo_17network')) {
      vattr_all <- vertex_attr(x, vattr)
      vattr_graph <- unique(vattr_all)
      atlas_vals <- atlas.dt[, get(vattr)]
      cnames <- intersect(levels(atlas_vals), vattr_graph)
      totals <- unname(table(atlas_vals)[cnames])
      if (vattr != 'class') {
        cnames <- paste0(cnames, ': ', table(vattr_all)[cnames], ' / ', totals)
      }
      group.nums <- which(levels(atlas_vals) %in% vattr_graph)
    } else if (vattr == 'rich') {
      group.nums <- 1L:3L
      cnames <- c('Rich-club', 'Feeder', 'Local')
    } else {
      group.nums <- as.integer(which(table(vertex_attr(x, vattr)) > 2))
      cnames <- paste0(simpleCap(vattr), '. ', group.nums)
    }
    nGroups <- length(group.nums)
    cex <- ifelse(nGroups > 10L, 0.75, ifelse(nGroups > 5L, 0.85, 1))
    cols <- group.cols[group.nums]
    l.inset <- c(-0.01, 0)  # Increasing 1st moves it toward the center (L-R), 2nd toward center (U-D)
    loc <- 'topleft'
    if (!is.null(label) && side == 3) loc <- 'topright'
    legend(loc, cnames, fill=cols, bg='black', text.col='white', cex=cex, inset=l.inset)
  }
}

#' Plot a slice of the MNI152 template
#'
#' Helper function that plots a single slice from the MNI152 template.
#' @inheritParams plot.brainGraph
#' @noRd

plot_mni <- function(plane, hemi) {
  if (plane == 'axial') {
    slice <- 46
    X <- mni152_mat
  } else {
    if (hemi == 'L') {
      X <- mni152_L
      slice <- 62
    } else if (hemi == 'R') {
      X <- mni152_R
      slice <- 30
    }
  }
  dims <- dim(X)
  x <- seq_len(dims[1L])
  y <- seq_len(dims[2L])
  zlim <- c(3500, max(X[, , slice]))
  imcol <- gray(0:64/64)
  breaks <- seq.int(min(zlim), max(zlim), length.out=length(imcol) + 1L)

  par(mai=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), bg='black', bty='o')
  graphics::image(x, y, X[, , slice], col=imcol, breaks=breaks, asp=1)
  par(new=TRUE)
}

#' Plot a graph with results from GLM-based analyses
#'
#' These methods are convenience functions for plotting a graph based on results
#' from GLM-based analyses (i.e., \code{\link{brainGraph_GLM},
#' \link{brainGraph_mediate}, \link{mtpc}, \link{NBS}}). There are several
#' default arguments which differ depending on the input object.
#'
#' The default arguments are specified so that the user only needs to type
#' \code{plot(x)} at the console, if desired. For all methods, the plot's
#' \emph{subtitle} will be omitted.
#'
#' @section NBS:
#' By default, a subgraph will be plotted consisting of only those vertices
#' which are part of a significant connected component. Vertex/edge colors will
#' correspond to connected component membership. Vertex names will be omitted.
#' Finally, the plot title will contain the contrast name.
#'
#' @section brainGraph_GLM:
#' By default, a subgraph will be plotted consisting of only those vertices for
#' which \eqn{p < \alpha}. It will also include a plot title with the outcome
#' measure and contrast name.
#'
#' @section mtpc:
#' By default, a subgraph will be plotted consisting of only those vertices for
#' which \eqn{A_{mtpc} > A_{crit}}. It will also include a plot title with the
#' outcome measure and contrast name.
#'
#' @section brainGraph_mediate:
#' By default, a subgraph will be plotted consisting of only those vertices for
#' which \eqn{P_{acme} < \alpha}. It will also include a plot title with the
#' treatment, mediator, and outcome variable names.
#'
#' @param x A \code{brainGraph_GLM}, \code{brainGraph_mtpc},
#'   \code{brainGraph_mediate}, or \code{brainGraph_NBS} object
#' @param alpha Numeric; the significance level. Default: \code{0.05}
#' @param subgraph Character string specifying the condition for subsetting the
#'   graph.
#' @param vertex.label Character vector of the vertex labels to be displayed.
#' @param vertex.color Character string specifying the vertex attribute to color
#'   the vertices by.
#' @param edge.color Character string specifying the edge attribute to color
#'   the edges by.
#' @param cex.main Numeric; the scaling factor for text size; see
#'   \code{\link[graphics]{par}}
#' @param ... Other arguments passed to \code{\link{plot.brainGraph}}
#' @inheritParams plot.brainGraph
#' @export
#' @name Plotting GLM graphs
#' @rdname glm_graph_plots
#' @family Plotting functions

plot.brainGraph_NBS <- function(x, alpha=0.05, subgraph=paste('p.nbs >', 1 - alpha),
                                vertex.label=NA,
                                vertex.color='color.comp', edge.color='color.comp', subtitle=NULL,
                                main=paste0('NBS: ', x$name), cex.main=2, ...) {
  NextMethod(generic='plot', object=x, subgraph=subgraph,
             vertex.label=vertex.label, vertex.color=vertex.color, edge.color=edge.color,
             main=main, subtitle=subtitle, cex.main=cex.main, ...)
}

#' @param p.sig Character string indicating which p-value to use for determining
#'   significance (default: \code{p})
#' @export
#' @rdname glm_graph_plots

plot.brainGraph_GLM <- function(x, p.sig=c('p', 'p.fdr', 'p.perm'),
                                subgraph=NULL,
                                main=paste0(x$outcome, ': ', x$name),
                                subtitle=NULL, cex.main=2, ...) {
  p.sig <- match.arg(p.sig)
  if (is.null(subgraph)) subgraph <- paste0(p.sig, ' > 1 -', x$alpha)
  NextMethod(generic='plot', object=x, subgraph=subgraph,
             main=main, subtitle=subtitle, cex.main=cex.main, ...)
}

#' @export
#' @rdname glm_graph_plots

plot.brainGraph_mtpc <- function(x, subgraph='sig == 1',
                                 main=paste0(x$outcome, ': ', x$name),
                                 subtitle=NULL, cex.main=2, ...) {
  NextMethod(generic='plot', object=x, subgraph=subgraph,
             main=main, subtitle=subtitle, cex.main=cex.main, ...)
}

#' @export
#' @rdname glm_graph_plots

plot.brainGraph_mediate <- function(x, subgraph='p.acme > 0.95',
          main=sprintf('Effect of "%s" on\n"%s"\nmediated by "%s"', x$treat, x$outcome, x$mediator),
          subtitle=NULL, cex.main=1, ...) {
  NextMethod(generic='plot', object=x, subgraph=subgraph,
             main=main, subtitle=subtitle, cex.main=cex.main, ...)
}
