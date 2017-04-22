#' Plot a brain graph with a specific spatial layout
#'
#' This function plots a graph when the spatial layout of the nodes is important
#' (e.g. in the brain). The function \code{\link{set_brainGraph_attr}}
#' needs to be run on the graph, and a valid set of coordinates provided for the
#' vertices. Most of the parameters valid here can be seen in
#' \code{\link{igraph.plotting}}.
#'
#' With the argument \code{subgraph}, you can specify a simple logical equation
#' for which vertices to show. For example, \emph{'degree > 10'} will plot only
#' vertices with a \emph{degree} greater than 10. Combinations of \emph{AND}
#' (i.e., \code{&}) and \emph{OR} (i.e., \code{|}) are allowed.
#'
#' To remove the subtitle at the bottom, simply specify \code{sub=NULL}.
#'
#' @param g An \code{igraph} graph object
#' @param plane A character string indicating which orientation to plot
#'   (default: \code{'axial'})
#' @param hemi A character string indicating which hemisphere to plot (default:
#'   \code{'both'})
#' @param subgraph A character string specifying an equation for deleting
#'   vertices (default: \code{NULL})
#' @param show.legend Logical indicating whether or not to show a legend
#'   (default: \code{FALSE})
#' @param rescale A logical, whether to rescale the coordinates (default:
#'   \code{FALSE})
#' @param asp A numeric constant for the aspect ratio (default: 0)
#' @param main Character string for the main title (default: \code{NULL})
#' @param sub Character string for the subtitle (default: \code{default})
#' @param ... Other parameters (passed to \code{\link{plot}}).
#' @export
#'
#' @family Plotting functions
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' plot_brainGraph(g[[1]], hemi='R')
#' plot_brainGraph(g[[1]], subgraph='degree > 10 | btwn.cent > 50')
#' }

plot_brainGraph <- function(g, plane=c('axial', 'sagittal', 'circular'),
                            hemi=c('both', 'L', 'R'),
                            subgraph=NULL, show.legend=FALSE,
                            rescale=FALSE, asp=0, main=NULL, sub='default', ...) {
  stopifnot(is_igraph(g), 'atlas' %in% graph_attr_names(g))
  lobe <- network <- NULL

  #---------------------------------------------------------
  # Create a subgraph based on user-specified condition
  #---------------------------------------------------------
  if (!is.null(subgraph)) {
    stopifnot(nchar(subgraph) > 0)
    # Function for creating the condition string to later subset the graph
    get_cond_string <- function(orig) {
      substrings <- strsplit(orig, split='\\s\\&\\s|\\s\\|\\s')[[1]]
      if (length(substrings) > 1) {  # Multiple conditions
        if (!isTRUE(grepl('\\s\\&\\s|\\s\\|\\s', orig))) {
          stop('Logical operators must be surrounded by spaces!')
        }
        nchars <- cumsum(sapply(substrings, nchar))
        splits <- sapply(seq_along(substrings), function(x)
                         substr(orig, start=nchars[x]+(3*x-1), stop=nchars[x]+(3*x-1)))
        substrings <- gsub('^\\s+|\\s+$', '', substrings) # Remove unnecessary whitespace

        cond.string <- paste(sapply(seq_along(substrings), function(x)
                                    paste0('V(g)$', substrings[x], splits[x])),
                             collapse='')
      } else {
        cond.string <- paste0('V(g)$', substrings)
      }
      return(cond.string)
    }

    # Handle when logical expressions are separated by parentheses
    if (isTRUE(grepl('\\(.*\\)', subgraph))) {
      subs <- strsplit(subgraph, split='\\)\\s\\&\\s\\(')[[1]]
      subs <- as.list(gsub('^\\(|\\)$|^\\s+|\\s+$', '', subs))
      cond.strings <- sapply(subs, get_cond_string)
      cond.string <- paste0('(', cond.strings[1], ') & (', cond.strings[2], ')')
    } else {
      cond.string <- get_cond_string(subgraph)
    }

    cond <- eval(parse(text=cond.string))
    cond <- setdiff(seq_len(vcount(g)), which(cond))
    g <- delete.vertices(g, cond)
  }

  atlas <- graph_attr(g, 'atlas')
  plane <- match.arg(plane)
  hemi <- match.arg(hemi)
  Nv <- vcount(g)
  adjust <- 0
  mult <- 100
  #---------------------------------------------------------
  # Handle plotting of individual hemispheres
  #---------------------------------------------------------
  if (hemi != 'both') {
    memb <- which(V(g)$hemi == hemi)
    sg <- induced.subgraph(g, memb)
    sg <- delete_all_attr(sg)
    V(sg)$name <- V(g)$name[memb]
    g <- g %s% sg
    g <- g - vertices(setdiff(seq_len(Nv), memb))

    if (plane == 'sagittal') {
      adjust <- 1
      V(g)$y <- V(g)$z.mni

      if (hemi == 'L') {
        V(g)$x <- -V(g)$y.mni
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
        V(g)$x <- V(g)$y.mni
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
    } else {
      mult <- 1
      xlim.g <- ylim.g <- c(-1.25, 1.25)
    }
  } else {
    if (plane == 'circular') {
      mult <- 1
      xlim.g <- ylim.g <- c(-1.25, 1.25)
    } else {
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

  #---------------------------------------------------------
  # Handle extra arguments in case a subgraph was created
  #---------------------------------------------------------
  fargs <- list(...)
  vcols <- 'lightblue'
  if (hasArg('vertex.color')) {
    if (is.character(fargs$vertex.color) & length(fargs$vertex.color) == 1) {
      stopifnot(fargs$vertex.color %in% vertex_attr_names(g))
      vcols <- vertex_attr(g, fargs$vertex.color)
      if (adjust == 1) {
        medial <- which(abs(V(g)$x.mni) < 20)
        vcols[medial] <- adjustcolor(vcols[medial], alpha.f=0.4)
      }
    } else {
      vcols <- fargs$vertex.color
    }
  }
  ecols <- 'red'
  if (hasArg('edge.color')) {
    if (is.character(fargs$edge.color) & length(fargs$edge.color) == 1) {
      stopifnot(fargs$edge.color %in% edge_attr_names(g))
      ecols <- edge_attr(g, fargs$edge.color)
      if (adjust == 1) {
        medial.es <- as.numeric(E(g)[medial %--% medial])
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
      stopifnot(fargs$vertex.size %in% vertex_attr_names(g))
      vsize <- mult * vertex_attr(g, fargs$vertex.size)
      if (any(!is.finite(vsize))) g <- delete.vertices(g, which(!is.finite(vsize)))
      vsize <- mult * vec.transform(vsize, as.numeric(min(vsize) > 0), 15)
    } else {
      vsize <- mult * fargs$vertex.size
    }
  }
  ewidth <- 1.5
  if (hasArg('edge.width')) {
    if (is.character(fargs$edge.width)) {
      stopifnot(fargs$edge.width %in% edge_attr_names(g))
      ewidth <- edge_attr(g, fargs$edge.width)
    } else {
      ewidth <- fargs$edge.width
    }
  }
  # Vertex labels
  #-------------------------------------
  vlabel <- V(g)$name
  vlabel.dist <- ifelse(vsize >= mult * 10, 0, 5)
  vlabel.col <- ifelse(vcols %in% c('red', 'blue', 'magenta'), 'white', 'blue')
  vlabel.font <- 2
  if (hasArg('vertex.label')) {
    if (is.na(fargs$vertex.label)) {
      vlabel <- vlabel.cex <- vlabel.dist <- vlabel.col <- vlabel.font <- NA
    } else {
      vlabel <- vertex_attr(g, fargs$vertex.label)
    }
  }
  if (hasArg('vertex.label.cex')) {
    vlabel.cex <- fargs$vertex.label.cex
  } else {
    vlabel.cex <- 1
  }

  plot(g, asp=asp, rescale=rescale, xlim=xlim.g, ylim=ylim.g,
       vertex.color=vcols, edge.color=ecols, vertex.size=vsize,
       edge.width=ewidth, vertex.label=vlabel, vertex.label.cex=vlabel.cex,
       vertex.label.dist=vlabel.dist, vertex.label.font=vlabel.font,
       vertex.label.color=vlabel.col, ...)

  if (!is.null(sub)) {
    if (sub == 'default') {
      Ne <- ecount(g)
      g.density <- round(graph.density(g), digits=3)
      par(new=TRUE, mar=c(5, 0, 3, 0)+0.1)
      sub <- paste('# vertices: ', Nv, '# edges: ', Ne, '\n',
                   'Density: ', g.density)
    }
  }
  if (hasArg('cex.main')) {
    cex.main <- fargs$cex.main
  } else {
    cex.main <- 2.5
  }
  title(main=main, sub=sub, col.main='white', col.sub='white', cex.main=cex.main)

  if (isTRUE(show.legend)) {
    if (hasArg('vertex.color')) {
      atlas.dt <- eval(parse(text=g$atlas))
      if (fargs$vertex.color == 'color.lobe') {
        lobes.g <- unique(V(g)$lobe)
        classnames <- intersect(levels(atlas.dt$lobe), lobes.g)
        total <- unname(atlas.dt[, table(lobe)][classnames])
        classnames <- paste0(classnames, ': ', table(V(g)$lobe)[classnames],
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
        networks.g <- unique(V(g)$network)
        classnames <- intersect(levels(atlas.dt$network), networks.g)
        total <- unname(atlas.dt[, table(network)][classnames])
        classnames <- paste0(classnames, ': ', table(V(g)$network)[classnames],
                             ' / ', total)
        cols <- group.cols[which(levels(atlas.dt$network) %in% networks.g)]
        cex <- vlabel.cex
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
