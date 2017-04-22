#' Function to dynamically plot a graph
#'
#' This function is called by \link{plot_brainGraph_gui} to update the GUI.
#'
#' @param plotDev A Cairo device for the plotting area
#' @param g An \code{igraph} graph object for the first plotting area
#' @param g2 An \code{igraph} graph object for the second plotting area
#' @param plotFunc A function specifying which type of plot to use
#' @param vsize.measure Character string of the name of the attribute for vertex
#'   scaling
#' @param ewidth.measure Character string of the name of the attribute for edge
#'   scaling
#' @param vertColor A GTK combo box for changing vertex colors
#' @param hemi A GTK combo box for plotting individual hemispheres
#' @param lobe Character string of the lobe name to select (comma-separated if
#'   multiple lobes are listed)
#' @param orient A GTK combo box for plotting a specific orientation
#' @param vertSize.min A GTK spin button for minimum vertex size
#' @param edgeWidth.min A GTK spin button for minimum edge width
#' @param edgeWidth.max A GTK spin button for maximum edge width
#' @param vertSize.const A GTK entry for constant vertex size
#' @param edgeWidth.const A GTK entry for constant width
#' @param vertLabels A GTK check button for showing vertex labels
#' @param showLegend A GTK check button for showing a legend
#' @param comms Integer vector of the communities to plot
#' @param neighb Character vector of vertex names for getting neighborhoods
#' @param neighbInd Integer vector of vertex indices of \code{neighb}
#' @param slider A GTK horizontal slider widget for changing edge curvature
#' @param vertSize.other A GTK entry for vertex size (other attributes)
#' @param edgeWidth.other A GTK entry for edge width (other attributes)
#' @param vertSize.eqn A GTK entry for equations to exclude vertices
#' @param showDiameter A GTK check button for showing the graph's diameter
#' @param edgeDiffs A GTK check button for showing edge diffs between graphs
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

update_brainGraph_gui <- function(plotDev, g, g2, plotFunc, vsize.measure, ewidth.measure,
                                  vertColor, hemi, lobe, orient, vertSize.min,
                                  edgeWidth.min, edgeWidth.max, vertSize.const,
                                  edgeWidth.const, vertLabels, showLegend, comms,
                                  neighb, neighbInd=NULL, slider, vertSize.other,
                                  edgeWidth.other, vertSize.eqn, showDiameter,
                                  edgeDiffs) {
  subgraph <- v.attr <- e.attr <- NULL
  dev.set(plotDev)
  Nv <- vcount(g)
  #=======================================================
  # Lobe to plot
  #=======================================================
  if (lobe != 'All') {
    subs <- strsplit(lobe, split=', ')[[1]]
    subgraph <- paste(paste0('lobe == "', subs, '"'), collapse=' | ')
  }

  #=======================================================
  # Hemisphere to plot
  #=======================================================
  plotHemi <- switch(hemi$getActive() + 1, 'both', 'L', 'R', rep('both', 6))
  if (hemi$getActive() > 2) {
    groups <- switch(hemi$getActive() - 2,
                     seq_len(Nv),
                     seq_len(Nv),
                     seq_len(max(V(g)$comm)),
                     seq_len(max(V(g)$comm)),
                     sort(unique(V(g)$lobe)),
                     sort(unique(V(g)$lobe)))
    eids <- switch(hemi$getActive() - 2,
                   E(g)[which(V(g)$hemi == 'L') %--% which(V(g)$hemi == 'R')],
                   count_homologous(g),
                   unique(unlist(sapply(groups, function(x)
                                        as.numeric(E(g)[which(V(g)$comm == x) %--%
                                                   which(V(g)$comm %in% groups[-x])])))),
                   unique(unlist(sapply(groups, function(x)
                                        as.numeric(E(g)[which(V(g)$comm == x) %--%
                                                   which(V(g)$comm == x)])))),
                   unique(unlist(sapply(groups, function(x)
                                        as.numeric(E(g)[which(V(g)$lobe == x) %--%
                                                   which(V(g)$lobe %in% groups[-x])])))),
                   unique(unlist(sapply(groups, function(x)
                                        as.numeric(E(g)[which(V(g)$lobe == x) %--%
                                                   which(V(g)$lobe == x)])))))
    sg.hemi <- subgraph.edges(g, eids)
    memb.hemi <- which(V(g)$name %in% V(sg.hemi)$name)
    sg.hemi <- delete_all_attr(sg.hemi)
    V(sg.hemi)$name <- V(g)$name[memb.hemi]
    g <- sg.hemi %s% g
  }

  #====================================================
  # Orientation of plots
  #====================================================
  plane <- switch(orient$getActive() + 1, 'axial', 'sagittal', 'sagittal', 'circular')
  plotHemi <- switch(orient$getActive() + 1, plotHemi, 'L', 'R', plotHemi)
  imSlice <- switch(orient$getActive() + 1, 46, 30, 30, NULL)
  if (orient$getActive() == 3) { # CIRCULAR LAYOUT
    par(bg='black')
    circ <- V(g)$circle.layout
    layout.g <- rotation(layout.circle(g, order=circ), -pi/2 - pi/Nv)
    V(g)$x <- layout.g[, 1]
    V(g)$y <- layout.g[, 2]
  } else {
    plot_brainGraph_mni(plane=plane, slice=imSlice, hemi=plotHemi)
  }
  par(pty='s', mar=rep(0, 4))

  #=======================================================
  # Vertex neighborhoods, if applicable
  if (!is.function(plotFunc) && plotFunc == 'plot_neighborhood') {
    if (length(neighbInd) == 1 && V(g)[neighbInd]$degree < 2) {
      g.sub <- make_ego_graph(g, order=1, nodes=neighbInd)[[1]]
    } else {
      g.sub <- make_ego_brainGraph(g, neighbInd)
      g.sub <- set_brainGraph_attr(g.sub, g$atlas)
    }
    g <- g.sub
    neighbInd <- which(V(g)$name %in% neighb)
    Nv <- vcount(g)
  }

  # Community number, if applicable
  if (!is.function(plotFunc) && plotFunc == 'plot_community') {
    subgraph <- paste(paste0('comm == ', comms), collapse=' | ')
  }
  plotFunc <- plot_brainGraph
  #=======================================================

  #-----------------------------------
  # Vertex sizes
  #-----------------------------------
  v.min <- vertSize.min$getValue()
  if (vertSize.const$getSensitive()) {
    vsize <- eval(parse(text=vertSize.const$getText()))
  } else {
    if (vsize.measure == 'other') {
      vsize <- vertSize.other$getText()
      if (!is.null(subgraph)) {
        subgraph <- paste0('(', subgraph, ') & (', vsize, ' >= ', v.min, ')')
      } else {
        subgraph <- paste0(vsize, ' >= ', v.min)
      }
    } else if (vsize.measure == 'eqn') {
      if (!is.null(subgraph)) {
        subgraph <- paste0('(', subgraph, ') & (', vertSize.eqn$getText(), ')')
      } else {
        subgraph <- vertSize.eqn$getText()
      }
      vsize <- 7.5
    } else {
      if (vsize.measure == 'hub.score' && !is_directed(g)) vsize.measure <- 'ev.cent'
      if (!is.null(subgraph)) {
        subgraph <- paste0('(', subgraph, ') & (', vsize.measure, ' >= ', v.min, ')')
      } else {
        subgraph <- paste0(vsize.measure, ' >= ', v.min)
      }
      vsize <- vsize.measure
    }
  }

  #-----------------------------------
  # Edge width
  #-----------------------------------
  e.min <- edgeWidth.min$getValue()
  e.max <- edgeWidth.max$getValue()

  if (edgeWidth.const$getSensitive()) {
    E(g)$width <- eval(parse(text=edgeWidth.const$getText()))
  } else {
    if (ewidth.measure == 'other') {
      e.attr <- edgeWidth.other$getText()
      e.vals <- edge_attr(g, e.attr)
      # If all are negative, take the absolute value
      if (sum(e.vals > 0) == 0) edge_attr(g, e.attr) <- abs(e.vals)
      g <- delete.edges(g, which(edge_attr(g, e.attr) < e.min))
      g <- delete.edges(g, which(edge_attr(g, e.attr) > e.max))
      E(g)$width <- vec.transform(edge_attr(g, e.attr), 0, 5)
    } else {
      g <- delete.edges(g, which(edge_attr(g, ewidth.measure) < e.min))
      g <- delete.edges(g, which(edge_attr(g, ewidth.measure) > e.max))
      E(g)$width <- switch(ewidth.measure,
                           btwn=log1p(E(g)$btwn),
                           dist=vec.transform(E(g)$dist, 0.1, 5),
                           weight=vec.transform(E(g)$weight, 1, 5))
    }
  }

  # Vertex & edge colors
  if (vertColor$getActive() == 0) {
    vertex.color <- rep('lightblue', Nv)
    vertex.color[neighbInd] <- 'yellow'
    edge.color <- rep('red', ecount(g))
  } else {
    edge.color <- vertex.color <- switch(vertColor$getActive(), 'color.comm', 'color.lobe',
                           'color.comp', 'color.comm.wt', 'color.class',
                           'color.network', 'color.nbhood')
    if (vertex.color == 'color.nbhood') {
      nbs.l <- lapply(neighbInd, function(x) neighbors(g, x))
      vertex.color <- rep('lightblue', Nv)
      for (i in seq_along(nbs.l)) {
        vertex.color[c(neighbInd[i], nbs.l[[i]])] <- group.cols[i]
      }
      edge.color <- set_edge_color(g, as.numeric(factor(vertex.color, levels=group.cols)))
    }
    #V(g)$color <- vertex_attr(g, vertex.color)
    #E(g)$color <- edge_attr(g, vertex.color)
  }

  # Show vertex labels?
  vlabel <- ifelse(vertLabels$active == FALSE, NA, 'name')

  # Slider for curvature of edges in circle plots
  curv <- ifelse(length(class(slider)) > 1, slider$getValue(), 0)

  main <- g$Group
  cex.main <- 2.5
  if (!is.null(neighbInd)) {
    main <- paste0('Neighborhoods of: ', paste(neighb, collapse=', '))
    cex.main <- 1
  }

  # Show a legend for vertex colors
  show.legend <- FALSE
  if (vertColor$getActive() %in% c(2, 5, 6) & showLegend$active == TRUE) {
    show.legend <- TRUE
  }

  plotFunc(g, plane=plane, hemi=plotHemi, subgraph=subgraph,
           vertex.label=vlabel, vertex.size=vsize, edge.width='width',
           vertex.color=vertex.color, edge.color=edge.color,
           edge.curved=curv, main=main, cex.main=cex.main, show.legend=show.legend)

  if (!is.null(showDiameter) && (showDiameter$active == TRUE | edgeDiffs$active == TRUE)) {
    if (showDiameter$active == TRUE) {
      inds <- as.numeric(get.diameter(g))
      es <- get.edge.ids(g, combn(inds, 2))
      es <- es[which(es != 0)]
      g.sub <- subgraph.edges(g, es)
    }
    if (edgeDiffs$active == TRUE) g.sub <- graph.difference(g, g2)
    plotFunc(g.sub, add=TRUE, vertex.label=NA,
             vertex.shape='none', edge.width=5,
             vertex.color='deeppink', edge.color='deeppink',
             edge.curved=curv, sub=NULL)
  }
}
