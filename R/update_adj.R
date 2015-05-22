#' Function to dynamically plot a graph
#'
#' This function is called by \link{plot.adj.gui} to update a plot on-the-fly.
#' It updates both of the plots by calling the helper function "make.plot".
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
#'
#' @export

update_adj <- function(graphname1, graphname2, vertLabels, vertSize,
                       edgeWidth, edgeDiffs, vertColor, firstplot, secondplot,
                       vertSize.const=NULL, edgeWidth.const=NULL, plotFunc,
                       comm=NULL, comboNeighb, hemi, lobe, orient1, orient2,
                       showDiameter, slider1, slider2, vertSize.other=NULL,
                       vertSize.min1, vertSize.min2, edgeWidth.min, kNumComms=NULL) {
  #===========================================================================
  #===========================================================================
  # Function that does all the work
  #===========================================================================
  #===========================================================================
  make.plot <- function(dev, g, g2, orient, comm=NULL, the.slider, kNumComms, vertSize.min, ...) {
    dev.set(dev)
    atlas <- g$atlas
    atlas.list <- eval(parse(text=atlas))
    Nv <- vcount(g)

    #=======================================================
    # Lobe to plot
    #=======================================================
    if (lobe$getActive() > 0) {
      kNumLobes <- nlevels(atlas.list[, lobe])
      combos <- sapply(seq_len(kNumLobes - 1),
                       function(x) combn(seq_along(levels(atlas.list[, lobe])), x))
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
      plot.over.brain.axial(0)
      #layout.g <- matrix(c(atlas.list$brainnet.coords[, 1],
      #                     atlas.list$brainnet.coords[, 2]),
      #                   ncol=2, byrow=F)
      #V(g)$x <- layout.g[, 1]
      #V(g)$y <- layout.g[, 2]
      #xlim.g <- c(-100, 100)
      #ylim.g <- c(-125, 85)
      #mult <- 100
      xlim.g <- c(-1, 1)
      ylim.g <- c(-1.5, 1.5)
      mult <- 1

    # Plot the left sagittal only
    #---------------------------------
    } else if (orient$getActive() == 1) {
      plot.over.brain.sagittal(0, hemi='left', z=30)
      cur.vertices <- intersect(which(memb.lobe), memb.hemi)
      mni.coords <- as.matrix(atlas.list[, c('x.mni', 'y.mni', 'z.mni'), with=F])
      layout.g <- matrix(c(-mni.coords[cur.vertices, 2],
                           mni.coords[cur.vertices, 3]),
                         ncol=2, byrow=F)
      V(g)$x <- layout.g[, 1]
      V(g)$y <- layout.g[, 2]
      xlim.g <- c(-85, 110)
      ylim.g <- c(-85, 125)
      mult <- 100

    # Plot the right sagittal only
    #---------------------------------
    } else if (orient$getActive() == 2) {
      plot.over.brain.sagittal(0, hemi='right', z=30)
      cur.vertices <- intersect(which(memb.lobe), memb.hemi)
      mni.coords <- as.matrix(atlas.list[, c('x.mni', 'y.mni', 'z.mni'), with=F])
      layout.g <- matrix(c(mni.coords[cur.vertices, 2],
                           mni.coords[cur.vertices, 3]),
                         ncol=2, byrow=F)
      V(g)$x <- layout.g[, 1]
      V(g)$y <- layout.g[, 2]
      xlim.g <- c(-125, 85)
      ylim.g <- c(-85, 125)
      mult <- 100

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
    }
    #=======================================================
    #=======================================================
    par(pty='s', mar=rep(0, 4))

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

    if (vertSize$getActive() == 0) {
      vsize <- mult * const
    } else if (vertSize$getActive() == 1) {
      vsize <- mult * vec.transform(V(g)$degree, 0, 15)
      g <- delete.vertices(g, V(g)$degree < v.min)
    } else if (vertSize$getActive() == 2) {
      g <- delete.vertices(g, V(g)$ev.cent < v.min)
      vsize <- mult * 15 * V(g)$ev.cent
    } else if (vertSize$getActive() == 3) {
      g <- delete.vertices(g, V(g)$btwn.cent < v.min)
      vsize <- mult * 3 * log1p(V(g)$btwn.cent)
    } else if (vertSize$getActive() == 4) {
      g <- delete.vertices(g, V(g)$subgraph.cent < v.min)
      vsize <- mult * 2 * log1p(V(g)$subgraph.cent)
    } else if (vertSize$getActive() == 5) {
      g <- delete.vertices(g, V(g)$coreness < v.min)
      vsize <- mult * V(g)$coreness
    } else if (vertSize$getActive() == 6) {
      g <- delete.vertices(g, V(g)$transitivity < v.min)
      vsize <- mult * 20 * V(g)$transitivity
    } else if (vertSize$getActive() == 7) {
      g <- delete.vertices(g, V(g)$PC < v.min)
      vsize <- mult * 25 * V(g)$PC
    } else if (vertSize$getActive() == 8) {
      g <- delete.vertices(g, V(g)$E.local < v.min)
      vsize <- mult * 15 * V(g)$E.local
    } else if (vertSize$getActive() == 9) {
      g <- delete.vertices(g, V(g)$E.nodal < v.min)
      vsize <- mult * 30 * V(g)$E.nodal
    } else if (vertSize$getActive() == 10) {
      g <- delete.vertices(g, V(g)$z.score < v.min)
      vsize <- mult*ifelse(V(g)$z.score == 0, 0,
                           vec.transform(V(g)$z.score, 0, 15))
    } else if (vertSize$getActive() == 11) {
      if (is.directed(g)) {
        g <- delete.vertices(g, V(g)$hub.score < v.min)
        vsize <- mult * 10 * sqrt(V(g)$hub.score)
      } else {
        g <- delete.vertices(g, V(g)$ev.cent < v.min)
        vsize <- mult * 15 * V(g)$ev.cent
      }
    } else if (vertSize$getActive() == 12) {
      g <- delete.vertices(g, V(g)$vulnerability < v.min)
      vsize <- mult * vec.transform(V(g)$vulnerability, 0, 15)
    } else if (vertSize$getActive() == 13) {
      g <- delete.vertices(g, which(vertex_attr(g, v.attr) < v.min))
      if (v.attr %in% c('p', 'p.adj', 'p.perm')) {
        vsize <- mult * 15 * vertex_attr(g, v.attr)
      } else {
        vsize <- mult * vertex_attr(g, v.attr)
      }
      #TODO: add lev.cent; same problem as w/ z.score though 
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
      ewidth <- 2
    }

    # Vertex & edge colors
    vertex.color <- switch(vertColor$getActive() + 1,
                           'lightblue',
                           V(g)$color.comm,
                           V(g)$color.lobe,
                           V(g)$color.comp
    )
    edge.color <- switch(vertColor$getActive() + 1,
                         'red',
                         E(g)$color.comm,
                         E(g)$color.lobe,
                         E(g)$color.comp
    )

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
      n <- comboNeighb$getActive() + 1
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

    if (identical(plotFunc, plot.adj)) {
      if (orient$getActive() == 3) {
        plotFunc <- plot
      }

      # Slider for curvature of edges in circle plots
      if (length(class(the.slider)) > 1) {
        curv <- the.slider$getValue()
      } else {
        curv <- 0
      }

      plotFunc(g,
               vertex.label=vlabel, vertex.label.cex=vlabel.cex,
               vertex.label.dist=vlabel.dist, vertex.label.font=vlabel.font,
               vertex.label.color=vlabel.color,
               vertex.size=vsize, edge.width=ewidth,
               vertex.color=vertex.color, edge.color=edge.color,
               xlim=xlim.g, ylim=ylim.g,
               #layout=layout.g,
               edge.curved=curv)
      if (orient$getActive() == 3) {
        g.density <- round(graph.density(g), digits=3)
        par(new=T, mar=c(5, 0, 3, 0) + 0.1)
        subt <- paste('# vertices: ', vcount(g), '# edges: ', ecount(g), '\n',
                      'Density: ', g.density)
        title(sub=subt, col.sub='white')
      }
    }

    # Show a legend for lobe colors
    if (vertColor$getActive() + 1 == 3) {
      lobes <- sort(unique(V(g)$lobe))
      lobe.names <- levels(atlas.list[, lobe])
      lobe.cols <- unique(V(g)$color.lobe[order(V(g)$lobe)])
      kNumLobes <- max(V(g)$lobe)
      legend('topleft',
             lobe.names,
             fill=lobe.cols,
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
               #layout=layout.g[inds, ],
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
               #layout=layout.g,
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
