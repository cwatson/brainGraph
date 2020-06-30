#' Function to dynamically plot a graph
#'
#' This function is called by \link{plot_brainGraph_gui} to update the GUI.
#'
#' @param g,g2 Graph objects for the plotting areas
#' @param orient A GTK combo box for plotting a specific orientation
#' @param hemi A GTK combo box for plotting individual hemispheres
#' @param vlabels A GTK check button for showing vertex labels
#' @param showLegend A GTK check button for showing a legend
#' @param mainTitle,subTitle A GTK check button for showing the main and
#'   subtitle, respectively
#' @param vcolor A GTK combo box for changing vertex colors
#' @param vsizeAttr Character string; attribute name for vertex scaling
#' @param vsizeConst A GTK entry for constant vertex size
#' @param vsizeMin A GTK spin button for minimum vertex size
#' @param vsizeEqn A GTK entry for equations to exclude vertices
#' @param ewidthAttr Character string; attribute name for edge scaling
#' @param ewidthConst A GTK entry for constant width
#' @param ewidthMin A GTK spin button for minimum edge width
#' @param ewidthMax A GTK spin button for maximum edge width
#' @param vGroupAttr A character string specifying which vertex groups to select
#'   from (e.g., \code{lobe}, \code{network}, etc.)
#' @param vgroups Character vector of the vertex groups to select; either
#'   lobe names, vertex names (for neighborhoods), communities, functional
#'   networks, cortical areas, gyri, or Yeo networks
#' @param slider A GTK horizontal slider widget for changing edge curvature
#' @param showDiameter A GTK check button for showing the graph's diameter
#' @param edgeDiffs A GTK check button for showing edge diffs between graphs
#' @noRd

update_brainGraph_gui <- function(g, orient, hemi,
                                  vlabels, showLegend, mainTitle, subTitle, vcolor,
                                  vsizeAttr, vsizeConst, vsizeMin, vsizeEqn,
                                  ewidthAttr, ewidthConst, ewidthMin, ewidthMax,
                                  vGroupAttr, vgroups=NULL,
                                  slider=NULL, showDiameter=NULL, edgeDiffs=NULL, g2=NULL) {
  subgraph <- NULL
  Nv <- vcount(g)
  main <- g$Group; cex.main <- 2.5
  #=======================================================
  # Hemisphere/vgroup to plot
  #=======================================================
  hemiInd <- hemi$getActive() + 1L
  hemiIndChar <- as.character(hemiInd)
  plotHemi <- switch(hemiIndChar, '2'='L', '3'='R', 'both')
  if (hemiInd > 3L) {
    vHemi <- switch(hemiIndChar,
                    '4'=, '5'=seq_len(Nv),                  # Inter-hemi & homologous
                    '6'=, '7'=seq_len(max(V(g)$comm)),      # Inter- & intra-comm
                    '8'=, '9'=sort(unique(V(g)$lobe)),      # Inter- & intra-lobe
                    '10'=, '11'=sort(unique(V(g)$network)), # Inter- & intra-network
                    '12'=, '13'=sort(unique(V(g)$area)),    # Inter- & intra-area (HCP)
                    '14'=, '15'=sort(unique(V(g)$Yeo_7network)),
                    '16'=, '17'=sort(unique(V(g)$Yeo_17network)))

    edges_intra <- function(g, vattr, vgroups) {
      vattrs <- vertex_attr(g, vattr)
      ids <- sapply(vgroups, function(x)
                    as.numeric(E(g)[which(vattrs == x) %--% which(vattrs == x)]))
      unique(unlist(ids))
    }
    edges_inter <- function(g, vattr, vgroups) {
      vattrs <- vertex_attr(g, vattr)
      ids <- sapply(vgroups, function(x)
                    as.numeric(E(g)[which(vattrs == x) %--% which(vattrs != x)]))
      unique(unlist(ids))
    }
    eids <- switch(hemiIndChar,
                   '4'=E(g)[which(V(g)$hemi == 'L') %--% which(V(g)$hemi == 'R')],
                   '5'=count_homologous(g),
                   '6'=unique(unlist(sapply(vHemi, function(x)
                                            as.numeric(E(g)[which(V(g)$comm == x) %--%
                                                       which(V(g)$comm %in% vHemi[-x])])))),
                   '7'=edges_intra(g, 'comm', vHemi),
                   '8'=edges_inter(g, 'lobe', vHemi),
                   '9'=edges_intra(g, 'lobe', vHemi),
                   '10'=edges_inter(g, 'network', vHemi),
                   '11'=edges_intra(g, 'network', vHemi),
                   '12'=edges_inter(g, 'area', vHemi),
                   '13'=edges_intra(g, 'area', vHemi),
                   '14'=edges_inter(g, 'Yeo_7network', vHemi),
                   '15'=edges_intra(g, 'Yeo_7network', vHemi),
                   '16'=edges_inter(g, 'Yeo_17network', vHemi),
                   '17'=edges_intra(g, 'Yeo_17network', vHemi))
    sg.hemi <- subgraph.edges(g, eids)
    memb.hemi <- which(V(g)$name %in% V(sg.hemi)$name)
    sg.hemi <- delete_all_attr(sg.hemi)
    V(sg.hemi)$name <- V(g)$name[memb.hemi]
    g <- sg.hemi %s% g
    class(g) <- c('brainGraph', class(g))
  }

  # Orientation of plots
  #====================================================
  plane <- switch(orient$getActive() + 1, 'axial', 'sagittal', 'sagittal', 'circular')
  plotHemi <- switch(orient$getActive() + 1, plotHemi, 'L', 'R', plotHemi)
  par(pty='s', mar=rep(0, 4))

  # Update "subgraph" for lobe, neighborhood, community, network, etc. plots
  #=======================================================
  if (vGroupAttr == 'neighborhood') {
    if (length(vgroups) == 1L && V(g)[vgroups]$degree < 2) {
      g.sub <- make_ego_graph(g, order=1, nodes=vgroups)[[1L]]
    } else {
      g.sub <- make_ego_brainGraph(g, vgroups)
    }
    class(g.sub) <- 'igraph'
    g <- make_brainGraph(g.sub, g$atlas)
    neighbInd <- which(V(g)$name %in% vgroups)
    Nv <- vcount(g)
    main <- paste0('Neighborhoods of: ', paste(vgroups, collapse=', '))
    cex.main <- 1

  } else if (vGroupAttr == 'comm') {
    subgraph <- paste0(vGroupAttr, ' == ', vgroups, collapse=' | ')
  } else {
    if (vgroups != 'All') {
      subs <- strsplit(vgroups, split=', ')[[1]]
      subgraph <- paste0(vGroupAttr, ' == "', subs, '"', collapse=' | ')
      main <- split_string(paste('\n', vGroupAttr, ':', paste(subs, collapse=', ')),
                           max_len=50L, delim=',')
      cex.main <- 1
    }
  }

  #-----------------------------------
  # Vertex sizes
  #-----------------------------------
  if (vsizeConst$getSensitive()) {
    vsize <- eval(parse(text=vsizeConst$getText()))
  } else {
    if (vsizeAttr == 'eqn') {
      vsize <- 7.5
      eqn <- vsizeEqn$getText()
    } else {
      v.min <- vsizeMin$getValue()
      v.attr <- vsize <- vsizeAttr
      eqn <- paste0(v.attr, ' >= ', v.min)
    }
    subgraph <- if (is.null(subgraph)) eqn else paste0('(', subgraph, ') & (', eqn, ')')
  }

  #-----------------------------------
  # Edge width
  #-----------------------------------
  if (ewidthConst$getSensitive()) {
    E(g)$width <- eval(parse(text=ewidthConst$getText()))
  } else {
    e.min <- ewidthMin$getValue()
    e.max <- ewidthMax$getValue()
    e.attr <- ewidthAttr
    e.vals <- edge_attr(g, e.attr)
    if (all(e.vals <= 0)) e.vals <- abs(e.vals)
    toDelete <- union(which(e.vals < e.min), which(e.vals > e.max))
    g <- delete_edges(g, toDelete)
    class(g) <- c('brainGraph', class(g))
    E(g)$width <- switch(ewidthAttr,
                         btwn=log1p(edge_attr(g, e.attr)),
                         dist=vec.transform(edge_attr(g, e.attr), 0.1, 5),
                         weight=vec.transform(edge_attr(g, e.attr), 1, 5),
                         vec.transform(edge_attr(g, e.attr), 0, 5))
  }

  # Vertex & edge colors
  edge.color <- rep('red', ecount(g))
  vertex.color <- rep('lightblue', Nv)
  if (vGroupAttr == 'neighborhood') vertex.color[neighbInd] <- 'yellow'
  if (vcolor$getActive() > 0) {
    edge.color <- vertex.color <-
      switch(vcolor$getActive(), 'color.comm', 'color.lobe', 'color.comp',
             'color.comm.wt', 'color.class', 'color.network', 'color.area',
             'color.Yeo_7network', 'color.Yeo_17network', 'color.nbhood')
    if (vertex.color == 'color.nbhood') {
      nbs.l <- lapply(neighbInd, function(x) neighbors(g, x))
      nbs.l <- nbs.l[lengths(nbs.l) > 0L]
      tab <- table(unlist(nbs.l))
      V(g)$nbhood <- 1
      for (i in seq_along(nbs.l)) V(g)[nbs.l[[i]]]$nbhood <- i + 1
      V(g)[as.numeric(names(tab[tab > 1]))]$nbhood <- i + 2
      V(g)[neighbInd]$nbhood <- 1 + 1:length(neighbInd)
      V(g)$nbhood <- V(g)$nbhood - 1
      g <- set_graph_colors(g, 'color.nbhood', V(g)$nbhood)
    }
  }

  # Show vertex labels?
  vlabel <- if (isTRUE(vlabels$active)) 'name' else NA

  # Slider for curvature of edges in circle plots
  curv <- if (length(class(slider)) > 1L) slider$getValue() else 0

  # Show a legend for vertex colors
  show.legend <- FALSE
  if (vcolor$getActive() %in% c(2, 5:10) && showLegend$active == TRUE) {
    show.legend <- TRUE
  }

  if (isFALSE(mainTitle$active)) main <- ''
  subTitle <- if (isTRUE(subTitle$active)) 'default' else NULL
  plot(g, plane=plane, hemi=plotHemi, subgraph=subgraph,
       vertex.label=vlabel, vertex.size=vsize, edge.width='width',
       vertex.color=vertex.color, edge.color=edge.color,
       edge.curved=curv, main=main, cex.main=cex.main, subtitle=subTitle, show.legend=show.legend)

  if (!is.null(showDiameter) && (showDiameter$active == TRUE || edgeDiffs$active == TRUE)) {
    if (showDiameter$active == TRUE) {
      inds <- as.numeric(get_diameter(g))
      es <- get.edge.ids(g, combn(inds, 2))
      es <- es[which(es != 0)]
      g.sub <- subgraph.edges(g, es)
    }
    if (edgeDiffs$active == TRUE) g.sub <- difference(g, g2)
    class(g.sub) <- c('brainGraph', class(g.sub))
    plot(g.sub, add=TRUE, vertex.label=NA, vertex.shape='none',
         vertex.size=vsize, vertex.color=rep('deeppink', vcount(g.sub)),
         edge.width=5, edge.color=rep('deeppink', ecount(g.sub)),
         edge.curved=curv, subtitle=NULL, mni=FALSE)
  }
}
