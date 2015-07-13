#' Function to dynamically plot a graph
#'
#' This function is called by \link{plot_brainGraph_gui} to update a plot
#' on-the-fly. It updates both of the plots by calling the helper function
#' "make.plot".
#'
#' @param graphname1 A GTK entry for the first plotting area
#' @param graphname2 A GTK entry for the second plotting area
#' @param vertLabels A GTK check button for showing vertex labels
#' @param vertSize A GTK combo box for scaling vertex size
#' @param edgeWidth A GTK entry for changing edge width
#' @param edgeDiffs A GTK check button for showing edge differences between
#' graphs
#' @param vertColor A GTK combo box for changing vertex colors
#' @param firstplot A Cairo device for the first plotting area
#' @param secondplot A Cairo device for the second plotting area
#' @param vertSize.const A GTK entry for constant vertex size
#' @param edgeWidth.const A GTK entry for constant width
#' @param plotFunc A function specifying which type of plot to use
#' @param comm A GTK combo box for plotting individual communities
#' @param comboNeighb A GTK combo box for plotting individual neighborhoods
#' @param hemi A GTK combo box for plotting individual hemispheres
#' @param lobe A GTK combo box for plotting individual lobes
#' @param orient1 A GTK combo box for plotting a specific orientation
#' @param orient2 A GTK combo box for plotting a specific orientation
#' @param showDiameter A GTK check button for showing the graph's diameter
#' @param slider1 A GTK horizontal slider widget for changing edge curvature
#' @param slider2 A GTK horizontal slider widget for changing edge curvature
#' @param vertSize.other A GTK entry for vertex size (other attributes)
#' @param vertSize.min1 A GTK entry for minimum vertex size, group 1
#' @param vertSize.min2 A GTK entry for minimum vertex size, group 2
#' @param edgeWidth.min A GTK entry for minimum edge width
#' @param kNumComms Integer indicating the number of total communities (optional)
#' @param comboLobeMult A GTK entry for joint neighborhoods of multiple vertices
#'
#' @export

update_brainGraph_gui <- function(graphname1, graphname2, vertLabels, vertSize,
                       edgeWidth, edgeDiffs, vertColor, firstplot, secondplot,
                       vertSize.const=NULL, edgeWidth.const=NULL, plotFunc,
                       comm=NULL, comboNeighb, hemi, lobe, orient1, orient2,
                       showDiameter, slider1, slider2, vertSize.other=NULL,
                       vertSize.min1, vertSize.min2, edgeWidth.min,
                       kNumComms=NULL, comboLobeMult=NULL) {
  #===========================================================================
  #===========================================================================
  # Function that does all the work
  #===========================================================================
  #===========================================================================
  make.plot <- function(dev, g, g2, orient, comm=NULL, the.slider, kNumComms, vertSize.min, ...) {
    dev.set(dev)
    atlas <- g$atlas
    atlas.dt <- eval(parse(text=data(list=atlas)))
    Nv <- vcount(g)

    #=======================================================
    # Lobe to plot
    #=======================================================
    if (lobe$getActive() > 0) {
      kNumLobes <- nlevels(atlas.dt[, lobe])
      combos <- sapply(seq_len(kNumLobes - 1),
                       function(x) combn(seq_along(levels(atlas.dt[, lobe])), x))
      n <- lobe$getActive()
      if (n <= kNumLobes) {
        ind1 <- 1
        ind2 <- n
      } else if (n > kNumLobes && n < 2*kNumLobes) {
        ind1 <- 2
        ind2 <- n - kNumLobes
      } else {
        ind1 <- which(cumsum(sapply(combos, ncol)) %/% n >= 1)[1]
        ind2 <- ncol(combos[[ind1]]) - (cumsum(sapply(combos, ncol)) %% n)[ind1]
      }
      memb.lobe <- apply(sapply(combos[[ind1]][, ind2],
                           function(x) V(g)$lobe == x), 1, any)
      sg.lobe <- induced.subgraph(g, memb.lobe)
    } else {
      memb.lobe <- rep(TRUE, Nv)
      sg.lobe <- g
    }

    for (att in graph_attr_names(sg.lobe)) sg.lobe <- delete_graph_attr(sg.lobe, att)
    for (att in vertex_attr_names(sg.lobe)) sg.lobe <- delete_vertex_attr(sg.lobe, att)
    for (att in edge_attr_names(sg.lobe)) sg.lobe <- delete_edge_attr(sg.lobe, att)
    V(sg.lobe)$name <- V(g)$name[memb.lobe]
    V(sg.lobe)$hemi <- V(g)$hemi[memb.lobe]

    #=======================================================
    # Hemisphere to plot
    #=======================================================
    if (hemi$getActive() == 0) {
      # Both hemispheres
      sg.hemi <- g
      memb.hemi <- seq_len(Nv)

    } else if (hemi$getActive() == 1) {
      # LH only
      memb.hemi <- which(V(g)$hemi == 'L')
      sg.hemi <- induced.subgraph(g, memb.hemi)

    } else if (hemi$getActive() == 2) {
      # RH only
      memb.hemi <- which(V(g)$hemi == 'R')
      sg.hemi <- induced.subgraph(g, memb.hemi)

    } else if (hemi$getActive() == 3) {
      # Interhemispheric only
      if (atlas %in% c('aal90', 'lpba40', 'hoa112')) {
        sg.hemi <- g - E(g) + subgraph.edges(g, E(g)[seq(1, Nv, 2) %--% seq(2, Nv, 2)])
      } else {
        sg.hemi <- g - E(g) + subgraph.edges(g, E(g)[1:(Nv/2) %--% (Nv/2 + 1):Nv])
      }
      memb.hemi <- seq_len(Nv)
    } else if (hemi$getActive() == 4) {
      # Homologous connections only
      eids <- count_homologous(g)
      sg.hemi <- subgraph.edges(g, eids)
      memb.hemi <- which(V(g)$name %in% V(sg.hemi)$name)
    }

    for (att in graph_attr_names(sg.hemi)) sg.hemi <- delete_graph_attr(sg.hemi, att)
    for (att in vertex_attr_names(sg.hemi)) sg.hemi <- delete_vertex_attr(sg.hemi, att)
    for (att in edge_attr_names(sg.hemi)) sg.hemi <- delete_edge_attr(sg.hemi, att)
    V(sg.hemi)$name <- V(g)$name[memb.hemi]
    V(sg.lobe)$hemi <- V(g)$hemi[memb.hemi]

    g <- g %s% (sg.lobe %s% sg.hemi)
    if (orient$getActive() != 3) {
      g <- g - vertices(setdiff(seq_len(Nv), intersect(which(memb.lobe), memb.hemi)))
    }

    #====================================================
    #====================================================
    # Orientation of plots
    #====================================================
    #====================================================
    if (orient$getActive() == 0) {
      plot_brainGraph_mni('axial')
      xlim.g <- c(-1, 1)
      ylim.g <- c(-1.5, 1.5)
      mult <- 1

      adjust <- 0

    # Plot the left sagittal only
    #---------------------------------
    } else if (orient$getActive() == 1) {
      plot_brainGraph_mni('sagittal', slice=30, hemi='left')
      V(g)$x <- -V(g)$y.mni
      V(g)$y <- V(g)$z.mni
      xlim.g <- c(-85, 110)
      ylim.g <- c(-85, 125)
      mult <- 100

      adjust <- 1

    # Plot the right sagittal only
    #---------------------------------
    } else if (orient$getActive() == 2) {
      plot_brainGraph_mni('sagittal', slice=30, hemi='right')
      V(g)$x <- V(g)$y.mni
      V(g)$y <- V(g)$z.mni
      xlim.g <- c(-125, 85)
      ylim.g <- c(-85, 125)
      mult <- 100

      adjust <- 1

    # Plot in a circular layout
    #---------------------------------
    } else if (orient$getActive() == 3) {
      par(bg='black')
      circ <- V(g)$circle.layout

      layout.g <- rotation(layout.circle(g, order=circ), -pi/2 - pi/Nv)
      V(g)$x <- layout.g[, 1]
      V(g)$y <- layout.g[, 2]
      xlim.g <- c(-1.25, 1.25)
      ylim.g <- c(-1.25, 1.25)
      mult <- 1
      adjust <- 0
    }
    #=======================================================
    #=======================================================
    par(pty='s', mar=rep(0, 4))


    # Multiple vertex neighborhoods, if applicable
    if (identical(plotFunc, plot_neighborhood)) {
      n <- comboNeighb$getActive()
      if (n == 0) {
        splits <- strsplit(comboLobeMult, split=', ')[[1]]
        if (any(grep('^[[:digit:]]*$', splits))) {  # numeric
          verts <- as.numeric(splits)
          vnames <- V(g)[verts]$name
        } else {
          verts <- c(splits)
          vnames <- verts
        }
        g.sub <- graph_neighborhood_multiple(g, verts)
        g.sub <- set.brainGraph.attributes(g.sub, atlas)
        g <- g.sub
        verts <- which(V(g)$name %in% vnames)
        Nv <- vcount(g)
        plotFunc <- plot_brainGraph
      }
    } else {
      n <- 1
      verts <- NULL
    }

    # Vertex sizes
    if (!vertSize.const$getSensitive()) {
      const <- NULL
      if (!vertSize.other$getSensitive()) {
        v.attr <- NULL
      } else {
        v.attr <- vertSize.other$getText()
      }
    } else {
      const <- eval(parse(text=vertSize.const$getText()))
      v.attr <- NULL
    }
    v.min <- eval(parse(text=vertSize.min$getText()))

    i <- vertSize$getActive()
    vsize.opts <- c('const', 'degree', 'ev.cent', 'btwn.cent', 'subgraph.cent',
                    'coreness', 'transitivity', 'PC', 'E.local', 'E.nodal',
                    'z.score', 'hub.score', 'vulnerability')
    if (i == 0) {
      vsize <- mult * const
    } else {
      if (i < 13) {
        if (i == 11 && !is.directed(g)) {
          i <- 2
        }
        g <- delete.vertices(g, which(vertex_attr(g,
                                vsize.opts[i + 1]) < v.min))
        vnum <- vertex_attr(g, vsize.opts[i + 1])
        if (vcount(g) %in% c(0, 1)) {
          vsize <- vnum
        } else {
          vsize <- mult * vec.transform(vnum, 0, 15)
        }

    } else if (vertSize$getActive() == 13) {
      g <- delete.vertices(g, which(vertex_attr(g, v.attr) < v.min))
      if (v.attr %in% c('p', 'p.adj', 'p.perm')) {
        vsize <- mult * 15 * vertex_attr(g, v.attr)
      } else {
        vsize <- mult * vertex_attr(g, v.attr)
      }
      #TODO: add lev.cent; same problem as w/ z.score though
    }
    }

    # Edge width
    if (!edgeWidth.const$getSensitive()) {
      e.const <- NULL
    } else {
      e.const <- eval(parse(text=edgeWidth.const$getText()))
    }
    e.min <- eval(parse(text=edgeWidth.min$getText()))

    if (edgeWidth$getActive() == 0) {
      ewidth <- e.const
    } else if (edgeWidth$getActive() == 1) {
      g <- delete.edges(g, which(E(g)$btwn < e.min))
      ewidth <- log1p(E(g)$btwn)
    } else if (edgeWidth$getActive() == 2) {
      g <- delete.edges(g, which(E(g)$dist < e.min))
      ewidth <- log1p(E(g)$dist)
    }

    # Vertex & edge colors
    vertex.color <- switch(vertColor$getActive() + 1,
                           {tmp <- rep('lightblue', Nv); tmp[verts] <- 'yellow'; tmp},
                           V(g)$color.comm,
                           V(g)$color.lobe,
                           V(g)$color.comp,
                           V(g)$color.class
    )
    edge.color <- switch(vertColor$getActive() + 1,
                         'red',
                         E(g)$color.comm,
                         E(g)$color.lobe,
                         E(g)$color.comp,
                         E(g)$color.class
    )

    # Make most medial colors more transparent if in sagittal view
    if (adjust == 1) {
      medial <- which(abs(V(g)$x.mni) < 20)
      vertex.color[medial] <- adjustcolor(vertex.color[medial], alpha.f=0.4)
      medial.es <- as.numeric(E(g)[medial %--% medial])
      edge.color[medial.es] <- adjustcolor(edge.color[medial.es], alpha.f=0.4)
    }

    # Show vertex labels?
    if (vertLabels$active == FALSE) {
      vlabel <- vlabel.cex <- vlabel.dist <- NA
      vlabel.color <- vlabel.font <- NA
    } else {
      vlabel <- V(g)$name
      vlabel.cex <- 0.75
      vlabel.dist <- ifelse(vsize >= 10, 0, 0.75)
      vlabel.color <- ifelse(vertex.color %in% c('red', 'blue'), 'white', 'blue')
      vlabel.font <- 2
    }


    # Community number, if applicable
    if (identical(plotFunc, plot_community)) {
      cNum <- comm$getActive() + 1
      combos <- sapply(seq_len(kNumComms), function(x)
                       combn(seq_len(kNumComms), x))
      if (cNum <= kNumComms) {
        ind1 <- 1
        ind2 <- cNum
      } else if (cNum > kNumComms && cNum < 2*kNumComms) {
        ind1 <- 2
        ind2 <- cNum - kNumComms
      } else {
        ind1 <- which(cumsum(sapply(combos, ncol)) %/% cNum >= 1)[1]
        ind2 <- ncol(combos[[ind1]]) - (cumsum(sapply(combos, ncol)) %% cNum)[ind1]
      }
      cNum <- combos[[ind1]][, ind2]
      plotFunc(g, n=cNum,
               vertex.label=vlabel, vertex.label.cex=vlabel.cex,
               vertex.label.dist=vlabel.dist, vertex.label.font=vlabel.font,
               vertex.label.color=vlabel.color,
               vertex.size=vsize, edge.width=ewidth,
               vertex.color=vertex.color, edge.color=edge.color)
    }

    # Vertex neighborhood, if applicable
    if (identical(plotFunc, plot_neighborhood)) {
      n <- comboNeighb$getActive()
      if (vertColor$getActive() == 0) {
        vertex.color <- V(g)$color <- rep('lightblue', Nv)
        vertex.color[n] <- V(g)[n]$color <- 'yellow'
        edge.color <- E(g)$color <- 'red'
      }
      plotFunc(g, n=n,
               vertex.label=vlabel, vertex.label.cex=vlabel.cex,
               vertex.label.dist=vlabel.dist, vertex.label.font=vlabel.font,
               vertex.label.color=vlabel.color,
               vertex.size=vsize, edge.width=ewidth,
               vertex.color=vertex.color, edge.color=edge.color)
    }

    if (identical(plotFunc, plot_brainGraph)) {
      if (orient$getActive() == 3) {
        plotFunc <- plot
      }

      # Slider for curvature of edges in circle plots
      if (length(class(the.slider)) > 1) {
        curv <- the.slider$getValue()
      } else {
        curv <- 0
      }

      if (n == 0) {
        main <- paste0('Neighborhoods of: ', paste(vnames, collapse=', '))
      } else {
        main <- NULL
      }
      plotFunc(g,
               vertex.label=vlabel, vertex.label.cex=vlabel.cex,
               vertex.label.dist=vlabel.dist, vertex.label.font=vlabel.font,
               vertex.label.color=vlabel.color,
               vertex.size=vsize, edge.width=ewidth,
               vertex.color=vertex.color, edge.color=edge.color,
               xlim=xlim.g, ylim=ylim.g,
               edge.curved=curv, main=main)

      if (orient$getActive() == 3) {
        g.density <- round(graph.density(g), digits=3)
        par(new=T, mar=c(5, 0, 3, 0) + 0.1)
        subt <- paste('# vertices: ', vcount(g), '# edges: ', ecount(g), '\n',
                      'Density: ', g.density)
        title(main=main, sub=subt, col.sub='white')
      }
    }

    # Show a legend for lobe colors
    if (vertColor$getActive() + 1 == 3) {
      lobes <- sort(unique(V(g)$lobe))
      lobe.names <- levels(atlas.dt[, lobe])[lobes]
      nonzero <- V(g)$degree > 0
      total <- unname(atlas.dt[, table(lobe)])[lobes]
      lobe.names <- paste0(lobe.names, ': ', table(V(g)$lobe[nonzero]),
                           ' / ', total)
      lobe.cols <- unique(V(g)$color.lobe[order(V(g)$lobe)])
      legend('topleft',
             lobe.names,
             fill=lobe.cols,
             bg='black',
             text.col='white')
    } else if (vertColor$getActive() + 1 == 5) {
      classes <- levels(atlas.dt[, class])
      cols <- c('red', 'green', 'blue')
      legend('topleft',
             classes,
             fill=cols,
             bg='black',
             text.col='white')
    }

    # Show the diameter of each graph?
    if (showDiameter$active == TRUE) {
      inds <- as.numeric(get.diameter(g))
      es <- get.edge.ids(g, combn(inds, 2))
      es <- es[which(es != 0)]
      g.sub <- subgraph.edges(g, es)
      #g.sub <- induced.subgraph(g, inds)
      plotFunc(g.sub, add=T, vertex.label=NA,
               vertex.size=10, edge.width=5,
               vertex.color='deeppink', edge.color='deeppink',
               xlim=xlim.g, ylim=ylim.g,
               edge.curved=curv)
    }

    # Show edge differences?
    if (edgeDiffs$active == TRUE) {
      g.diff12 <- graph.difference(g, g2)
      plotFunc(g.diff12, add=T,
               vertex.label=NA, vertex.size=degree(g.diff12),
               vertex.color='deeppink', edge.color='deeppink',
               edge.width=5,
               xlim=xlim.g, ylim=ylim.g,
               edge.curved=curv)
    }
  }
  #=============================================================================
  graph1 <- eval(parse(text=graphname1$getText()))
  graph2 <- eval(parse(text=graphname2$getText()))
  if (!is.igraph(graph1)) {
    stop(sprintf('%s is not a graph object.', graphname1$getText()))
  }

  make.plot(dev=firstplot, g=graph1, g2=graph2, orient=orient1,
            comm, vertLabels, vertSize, edgeWidth,
            vertColor, vertSize.const=NULL, hemi, the.slider=slider1,
            kNumComms=kNumComms, vertSize.min=vertSize.min1)

  if (nchar(graphname2$getText()) > 0) {
    if (!is.igraph(graph2)) {
      stop(sprintf('%s is not a graph object.', graphname2$getText()))
    }
    make.plot(dev=secondplot, g=graph2, g2=graph1, orient=orient2,
              comm, vertLabels, vertSize, edgeWidth,
              vertColor, vertSize.const=NULL, hemi, the.slider=slider2,
              kNumComms=kNumComms, vertSize.min=vertSize.min2)
  }
}
