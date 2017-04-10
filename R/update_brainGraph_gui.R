#' Function to dynamically plot a graph
#'
#' This function is called by \link{plot_brainGraph_gui} to update a plot
#' on-the-fly. It updates by calling the helper function \code{make.plot}.
#'
#' @param plotDev A Cairo device for the plotting area
#' @param graph1 An igraph graph object for the first plotting area
#' @param graph2 An igraph graph object for the  second plotting area
#' @param plotFunc A function specifying which type of plot to use
#' @param vertSize A GTK combo box for scaling vertex size
#' @param edgeWidth A GTK entry for changing edge width
#' @param vertColor A GTK combo box for changing vertex colors
#' @param hemi A GTK combo box for plotting individual hemispheres
#' @param lobe A GTK combo box for plotting individual lobes
#' @param orient A GTK combo box for plotting a specific orientation
#' @param vertSize.min A GTK spin button for minimum vertex size
#' @param edgeWidth.min A GTK spin button for minimum edge width
#' @param edgeWidth.max A GTK spin button for maximum edge width
#' @param vertSize.const A GTK entry for constant vertex size
#' @param edgeWidth.const A GTK entry for constant width
#' @param vertLabels A GTK check button for showing vertex labels
#' @param showLegend A GTK check button for showing a legend
#' @param comm A GTK combo box for plotting individual communities
#' @param kNumComms Integer indicating the number of total communities (optional)
#' @param neighb A GTK combo box for plotting individual neighborhoods
#' @param neighbMult A GTK entry for joint neighborhoods of multiple vertices
#' @param slider A GTK horizontal slider widget for changing edge curvature
#' @param vertSize.other A GTK entry for vertex size (other attributes)
#' @param edgeWidth.other A GTK entry for edge width (other attributes)
#' @param vertSize.eqn A GTK entry for equations to exclude vertices
#' @param showDiameter A GTK check button for showing the graph's diameter
#' @param edgeDiffs A GTK check button for showing edge diffs between graphs
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

update_brainGraph_gui <- function(plotDev, graph1, graph2, plotFunc, vertSize,
                       edgeWidth, vertColor, hemi, lobe, orient, vertSize.min,
                       edgeWidth.min, edgeWidth.max, vertSize.const=NULL,
                       edgeWidth.const=NULL, vertLabels=NULL, showLegend=NULL,
                       comm=NULL, kNumComms=NULL, neighb=NULL, neighbMult=NULL,
                       slider=NULL, vertSize.other=NULL, edgeWidth.other=NULL,
                       vertSize.eqn=NULL, showDiameter=NULL, edgeDiffs=NULL) {

  #===========================================================================
  #===========================================================================
  # Function that does all the work
  #===========================================================================
  #===========================================================================
  make.plot <- function(dev, g, g2, orient, vertSize.min, edgeWidth.min,
                        edgeWidth.max, vertSize.const=NULL, vertSize.eqn=NULL,
                        the.slider=NULL, kNumComms=NULL, comm=NULL, ...) {
    dev.set(dev)
    atlas.dt <- eval(parse(text=g$atlas))
    Nv <- vcount(g)
    #=======================================================
    # Lobe to plot
    #=======================================================
    memb.lobe <- rep(TRUE, Nv)
    sg.lobe <- g
    if (lobe$getActive() > 0) {
      kNumLobes <- nlevels(atlas.dt[, lobe])
      combos <- lapply(seq_len(kNumLobes - 1),
                       function(x) combn(seq_along(levels(atlas.dt[, lobe])), x))
      n <- lobe$getActive()
      if (n <= kNumLobes) {
        ind1 <- 1
        ind2 <- n
      } else if (n > kNumLobes && n < 2*kNumLobes) {
        ind1 <- 2
        ind2 <- n - kNumLobes
      } else {
        ind1 <- which(cumsum(vapply(combos, ncol, integer(1))) %/% n >= 1)[1]
        ind2 <- ncol(combos[[ind1]]) - (cumsum(vapply(combos, ncol, integer(1))) %% n)[ind1]
      }
      memb.lobe <- apply(vapply(combos[[ind1]][, ind2],
                                function(x) V(g)$lobe == x, logical(Nv)), 1, any)
      sg.lobe <- induced.subgraph(g, memb.lobe)
    }
    sg.lobe <- delete_all_attr(sg.lobe)
    V(sg.lobe)$name <- V(g)$name[memb.lobe]

    #=======================================================
    # Hemisphere to plot
    #=======================================================
    plotHemi <- 'both'
    if (hemi$getActive() == 1) {  # LH only
      plotHemi <- 'L'
    } else if (hemi$getActive() == 2) {  # RH only
      plotHemi <- 'R'

    } else if (hemi$getActive() > 2) {
      groups <- switch(hemi$getActive() - 2,
                       seq_len(Nv),
                       seq_len(Nv),
                       seq_len(max(V(g)$comm)),
                       seq_len(max(V(g)$comm)),
                       sort(unique(V(g)$lobe)),
                       sort(unique(V(g)$lobe)))
      eids <- switch(hemi$getActive() - 2,
                     E(g)[which(atlas.dt$hemi == 'L') %--% which(atlas.dt$hemi == 'R')],
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
      g <- g %s% (sg.lobe %s% sg.hemi)
    }

    if (lobe$getActive() > 0) {
      g <- g %s% sg.lobe
#    g <- g %s% (sg.lobe %s% sg.hemi)
    if (orient$getActive() != 3) {
      g <- g - vertices(setdiff(seq_len(Nv), which(memb.lobe)))#intersect(which(memb.lobe), memb.hemi)))
    }
    }

    #====================================================
    # Orientation of plots
    #====================================================
    if (orient$getActive() == 3) { # CIRCULAR LAYOUT
      plane <- 'circular'
      par(bg='black')
      circ <- V(g)$circle.layout

      layout.g <- rotation(layout.circle(g, order=circ), -pi/2 - pi/Nv)
      V(g)$x <- layout.g[, 1]
      V(g)$y <- layout.g[, 2]
      xlim.g <- ylim.g <- c(-1.25, 1.25)
    } else {

      if (orient$getActive() == 0) {
        plane <- 'axial'
        imSlice <- 46
      } else {
        plane <- 'sagittal'
        imSlice <- 30
        if (orient$getActive() == 1) {
          plotHemi <- 'L'
        } else if (orient$getActive() == 2) {
          plotHemi <- 'R'
        }
      }
      plot_brainGraph_mni(plane=plane, slice=imSlice, hemi=plotHemi)
    }

    #=======================================================
    #=======================================================
    par(pty='s', mar=rep(0, 4))

    # Vertex neighborhoods, if applicable
    if (!is.function(plotFunc) && plotFunc == 'plot_neighborhood') {
      n <- neighb$getActive()
      if (n == 0) { # Multiple vertices
        splits <- strsplit(neighbMult, split=', ')[[1]]
        if (any(grep('^[[:digit:]]*$', splits))) {  # numeric
          verts <- as.numeric(splits)
          vnames <- V(g)[verts]$name
        } else {
          vnames <- verts <- c(splits)
        }
      } else {
        verts <- n
        vnames <- V(g)[verts]$name
      }
      if (n > 0 && V(g)[n]$degree < 2) {
        g.sub <- make_ego_graph(g, order=1, nodes=n)[[1]]
        n <- 0
      } else {
        g.sub <- make_ego_brainGraph(g, verts)
        g.sub <- set_brainGraph_attr(g.sub, g$atlas)
      }
      g <- g.sub
      verts <- which(V(g)$name %in% vnames)
      Nv <- vcount(g)
      plotFunc <- plot_brainGraph
      n <- 0
    } else {
      n <- 1
      verts <- NULL
    }

    # Community number, if applicable
    if (!is.function(plotFunc) && plotFunc == 'plot_community') {
      cNum <- comm$getActive() + 1
      combos <- lapply(seq_len(kNumComms), function(x)
                       combn(seq_len(kNumComms), x))
      if (cNum <= kNumComms) {
        ind1 <- 1
        ind2 <- cNum
      } else if (cNum > kNumComms && cNum < 2*kNumComms) {
        ind1 <- 2
        ind2 <- cNum - kNumComms
      } else {
        ind1 <- which(cumsum(vapply(combos, ncol, integer(1))) %/% cNum >= 1)[1]
        ind2 <- ncol(combos[[ind1]]) -
            (cumsum(vapply(combos, ncol, integer(1))) %% cNum)[ind1]
      }
      cNums <- combos[[ind1]][, ind2]
      comms <- seq_len(max(V(g)$comm))[cNums]
      memb <- which(V(g)$comm %in% comms)
      g.sub <- induced.subgraph(g, memb)
      vcomms <- V(g.sub)$comm
      eids <- lapply(seq_along(comms), function(x)
                     as.numeric(E(g.sub)[which(vcomms == comms[x]) %--%
                                         which(vcomms == comms[x])]))
      if (length(cNums) == 1) {
        ecols <- rep(group.cols[cNums], length=ecount(g.sub))
      } else {
        listcols <- as.list(group.cols[cNums])
        ecols <- rep('gray50', length=ecount(g.sub))
        for (i in seq_along(eids)) ecols[eids[[i]]] <- listcols[[i]]
      }
      g <- g.sub
      E(g)$color.comm <- ecols
      plotFunc <- plot_brainGraph
    }

    #-----------------------------------
    # Vertex sizes
    #-----------------------------------
    subgraph <- const <- v.attr <- e.const <- e.attr <- NULL
    if (vertSize.const$getSensitive()) const <- eval(parse(text=vertSize.const$getText()))
    if (vertSize.other$getSensitive()) v.attr <- vertSize.other$getText()
    v.min <- vertSize.min$getValue()

    i <- vertSize$getActive()
    vsize.opts <- c('const', 'degree', 'ev.cent', 'btwn.cent',
                    'coreness', 'transitivity', 'PC', 'E.local', 'E.nodal',
                    'z.score', 'hub.score', 'vulnerability', 'knn', 'asymm',
                    'eccentricity', 'dist', 'dist.strength', 'Lp', 'other',
                    'eqn', 'strength', 'knn.wt', 'E.local.wt', 'E.nodal.wt')
    if (i == 0) {
      vsize <- const
    } else {
      if (i == 18) {  # Other
        subgraph <- paste0(v.attr, ' >= ', v.min)
        vsize <- v.attr
      } else if (i == 19) {  # equation
        subgraph <- vertSize.eqn$getText()
        vsize <- 7.5
      } else {
        if (i == 11 && !is_directed(g)) i <- 2
        subgraph <- paste0(vsize.opts[i+1], ' >= ', v.min)
        vsize <- vsize.opts[i+1]
      }
    }

    #-----------------------------------
    # Edge width
    #-----------------------------------
    if (edgeWidth.const$getSensitive()) e.const <- eval(parse(text=edgeWidth.const$getText()))
    if (edgeWidth.other$getSensitive()) e.attr <- edgeWidth.other$getText()
    e.min <- edgeWidth.min$getValue()
    e.max <- edgeWidth.max$getValue()

    j <- edgeWidth$getActive()
    ewidth.opts <- c('const', 'btwn', 'dist', 'other', 'weight')
    ewidth <- e.const
    if (j > 0 & j == 3) {
      e.vals <- edge_attr(g, e.attr)
      # If all are negative, take the absolute value
      if (sum(e.vals > 0) == 0) edge_attr(g, e.attr) <- abs(e.vals)
      g <- delete.edges(g, which(edge_attr(g, e.attr) < e.min))
      g <- delete.edges(g, which(edge_attr(g, e.attr) > e.max))
      ewidth <- vec.transform(edge_attr(g, e.attr), 0, 5)
    } else {
      g <- delete.edges(g, which(edge_attr(g, ewidth.opts[j+1]) < e.min))
      g <- delete.edges(g, which(edge_attr(g, ewidth.opts[j+1]) > e.max))

      if (j == 1) ewidth <- log1p(E(g)$btwn)
      if (j == 2) ewidth <- vec.transform(E(g)$dist, 0.1, 5)
      if (j == 4) ewidth <- vec.transform(E(g)$weight, min(E(g)$weight), 5)
    }

    # Vertex & edge colors
    if (vertColor$getActive() == 0) {
      V(g)$color <- rep('lightblue', Nv)
      V(g)$color[verts] <- 'yellow'
      E(g)$color <- rep('red', ecount(g))
    } else {
      vertex.color <- switch(vertColor$getActive(), 'color.comm', 'color.lobe',
                             'color.comp', 'color.comm.wt', 'color.class',
                             'color.network')
      V(g)$color <- vertex_attr(g, vertex.color)
      E(g)$color <- edge_attr(g, vertex.color)
    }

    # Show vertex labels?
    vlabel <- ifelse(vertLabels$active == FALSE, NA, 'name')

    # Slider for curvature of edges in circle plots
    curv <- ifelse(length(class(the.slider)) > 1, the.slider$getValue(), 0)

    if (n == 0) {
      main <- paste0('Neighborhoods of: ', paste(vnames, collapse=', '))
      cex.main <- 1
    } else {
      main <- g$Group
      cex.main <- 2.5
    }

    # Show a legend for vertex colors
    show.legend <- FALSE
    if (vertColor$getActive() %in% c(2, 5, 6) & showLegend$active == TRUE) {
      show.legend <- TRUE
    }

    plotFunc(g, plane=plane, hemi=plotHemi, subgraph=subgraph,
             vertex.label=vlabel, vertex.size=vsize, edge.width=ewidth,
             vertex.color='color', edge.color='color',
             edge.curved=curv, main=main, cex.main=cex.main, show.legend=show.legend)



    if (showDiameter$active == TRUE | edgeDiffs$active == TRUE) {
      # Show the diameter of each graph?
      if (showDiameter$active == TRUE) {
        inds <- as.numeric(get.diameter(g))
        es <- get.edge.ids(g, combn(inds, 2))
        es <- es[which(es != 0)]
        g.sub <- subgraph.edges(g, es)
      }
      # Show edge differences?
      if (edgeDiffs$active == TRUE) g.sub <- graph.difference(g, g2)
      plotFunc(g.sub, add=T, vertex.label=NA,
               vertex.shape='none', edge.width=5,
               vertex.color='deeppink', edge.color='deeppink',
               xlim=xlim.g, ylim=ylim.g,
               edge.curved=curv, sub=NULL)
    }
  }
  #=============================================================================
  make.plot(dev=plotDev, g=graph1, g2=graph2, orient, vertSize.min,
            edgeWidth.min, edgeWidth.max, vertLabels, vertSize, edgeWidth,
            vertColor, showLegend, vertSize.const=vertSize.const,
            vertSize.eqn=vertSize.eqn, hemi, the.slider=slider,
            kNumComms=kNumComms, comm=comm)
}
