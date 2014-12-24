#' Function to dynamically plot a graph
#'
#' This function is called by \link{plot.adj.gui} to update a plot on-the-fly.
#' It updates both of the plots by calling the helper function "make.plot".

update.adj <- function(graphname1, graphname2, vertLabels, vertSize,
                       edgeWidth, edgeDiffs, vertColor, firstplot, secondplot,
                       vertSize.const=NULL, plotFunc, comm=NULL, comboNeighb, hemi,
                       lobe, orient1, orient2, showDiameter,
                       slider1, slider2) {
  #===========================================================================
  #===========================================================================
  # Function that does all the work
  #===========================================================================
  #===========================================================================
  make.plot <- function(dev, g, g2, orient, comm=NULL, the.slider, ...) {
    dev.set(dev)
    atlas <- g$atlas
    atlas.list <- eval(parse(text=atlas))

    #=================================
    # Lobe to plot
    #=================================
    if (lobe$getActive() > 0) {
      kNumLobes <- length(atlas.list$lobe.names)
      combos <- sapply(seq_len(kNumLobes - 1),
                       function(x) combn(seq_along(atlas.list$lobe.names), x))
#      inds.all <- lapply(seq_along(combos),
#             function(z) sapply(seq_len(ncol(combos[[z]])),
#                    function(x) apply(sapply(seq_len(z),
#                             function(y) V(g)$lobe == t(combos[[z]])[x, y]), 1, any)))
#      inds.all <- do.call('cbind', inds.all)
#      g <- induced.subgraph(g, inds.all[, lobe$getActive()])

      n <- lobe$getActive()
      if (n <= kNumLobes) {
        FIRST <- 1
        SECOND <- n
      } else if (n > kNumLobes && n < 2*kNumLobes) {
        FIRST <- 2
        SECOND <- n - kNumLobes
      } else {
        FIRST <- which(cumsum(sapply(combos, ncol)) %/% n >= 1)[1]
        SECOND <- ncol(combos[[FIRST]]) - (cumsum(sapply(combos, ncol)) %% n)[FIRST]
      }
      memb.alt <- apply(sapply(combos[[FIRST]][, SECOND], function(x) V(g)$lobe == x), 1, any)
      g <- induced.subgraph(g, memb.alt)
    }
    Nv <- vcount(g)

    #====================================================
    #====================================================
    # Orientation of plots
    #====================================================
    #====================================================
    if (orient$getActive() == 0) {
      plot.over.brain.axial(0)
      # Hemisphere to plot
      if (atlas == 'aal90' || atlas == 'lpba40' || atlas == 'hoa112') {
        g <- switch(hemi$getActive()+1,
          g,
          induced.subgraph(g, seq(1, Nv, 2)),
          induced.subgraph(g, seq(2, Nv, 2)),
          subgraph.edges(g, E(g)[seq(1, Nv, 2) %--% seq(2, Nv, 2)])
        )
      } else {
        g <- switch(hemi$getActive()+1,
          g,
          induced.subgraph(g, 1:(Nv/2)),
          induced.subgraph(g, (Nv/2 + 1):Nv),
          subgraph.edges(g, E(g)[1:(Nv/2) %--% (Nv/2 + 1):Nv])
        )
      }
      layout.g <- cbind(V(g)$x, V(g)$y)
      xlim.g <- c(-1, 1)
      ylim.g <- c(-1.5, 1.5)
      mult <- 1

    #---------------------------------
    # Plot the left sagittal only
    #---------------------------------
    } else if (orient$getActive() == 1) {
      plot.over.brain.sagittal(0, hemi='left', z=30)
      if (atlas == 'aal90' || atlas == 'lpba40' || atlas == 'hoa112') {
        g <- induced.subgraph(g, seq(1, Nv, 2))
        layout.g <- matrix(c(-atlas.list$brainnet.coords[, 2],
                             atlas.list$brainnet.coords[, 3]),
                           ncol=2, byrow=F)
      } else {
        g <- induced.subgraph(g, 1:(Nv/2))
        layout.g <- matrix(c(-atlas.list$brainnet.coords[, 2],
                             atlas.list$brainnet.coords[, 3]),
                           ncol=2, byrow=F)
      }
      xlim.g <- c(-85, 110)
      ylim.g <- c(-85, 125)
      mult <- 100

    #---------------------------------
    # Plot the right sagittal only
    #---------------------------------
    } else if (orient$getActive() == 2) {
      plot.over.brain.sagittal(0, hemi='right', z=30)
      if (atlas == 'aal90' || atlas == 'lpba40' || atlas == 'hoa112') {
        g <- induced.subgraph(g, seq(2, Nv, 2))
        layout.g <- matrix(c(atlas.list$brainnet.coords[, 2],
                             atlas.list$brainnet.coords[, 3]),
                           ncol=2, byrow=F)
      } else {
        g <- induced.subgraph(g, ((Nv/2 + 1):Nv))
        layout.g <- matrix(c(atlas.list$brainnet.coords[, 2],
                             atlas.list$brainnet.coords[, 3]),
                           ncol=2, byrow=F)
      }
      xlim.g <- c(-125, 85)
      ylim.g <- c(-85, 125)
      mult <- 100

    #---------------------------------
    # Plot in a circular layout
    #---------------------------------
    } else if (orient$getActive() == 3) {
      par(bg='white')
      # Hemisphere to plot
      circ <- V(g)$circle.layout
      lobe <- V(g)$lobe

      # Both hemispheres
      if (hemi$getActive()+1 == 1) {
        sg <- g

      } else if (hemi$getActive()+1 == 2) {
        # LH only
        if (atlas == 'aal90' || atlas == 'lpba40' || atlas == 'hoa112') {
          sg <- g - E(g) + subgraph.edges(g, E(g)[seq(1, Nv, 2) %--% seq(1, Nv, 2)])
        } else {
          sg <- g - E(g) + subgraph.edges(g, E(g)[1:(Nv/2) %--% 1:(Nv/2)])
        }

      } else if (hemi$getActive()+1 == 3) {
        # RH only
        if (atlas == 'aal90' || atlas == 'lpba40' || atlas == 'hoa112') {
          sg <- g - E(g) + subgraph.edges(g, E(g)[seq(2, Nv, 2) %--% seq(2, Nv, 2)])
        } else {
          sg <- g - E(g) + subgraph.edges(g, E(g)[(Nv/2 + 1):Nv %--% (Nv/2 + 1):Nv])
        }

      } else if (hemi$getActive()+1 == 4) {
        # Interhemispheric only
        if (atlas == 'aal90' || atlas == 'lpba40' || atlas == 'hoa112') {
          sg <- g - E(g) + subgraph.edges(g, E(g)[seq(1, Nv, 2) %--% seq(2, Nv, 2)])
        } else {
          sg <- g - E(g) + subgraph.edges(g, E(g)[1:(Nv/2) %--% (Nv/2 + 1):Nv])
        }
      }

      for (att in list.graph.attributes(sg)) {
        sg <- remove.graph.attribute(sg, att)
      }
      for (att in list.vertex.attributes(sg)) {
        sg <- remove.vertex.attribute(sg, att)
      }
      for (att in list.edge.attributes(sg)) {
        sg <- remove.edge.attribute(sg, att)
      }

      V(sg)$name <- V(g)$name
      g <- graph.intersection(g, sg)
      layout.g <- rotation(layout.circle(g, order=circ), -pi/2)
      mult <- 1
      xlim.g <- c(-1.25, 1.25)
      ylim.g <- c(-1.25, 1.25)
    }
    #====================================================
    #====================================================
    par(pty='s', mar=rep(0, 4))

    # Show vertex labels?
    if (vertLabels$active == FALSE) {
      vertex.label <- NA
      vertex.label.cex <- NA
    } else {
      vertex.label <- V(g)$name
      vertex.label.cex <- 0.75
    }
    # Vertex & edge colors
    vertex.color <- switch(vertColor$getActive() + 1,
                           'lightblue3',
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

    # Vertex sizes
    if (!vertSize.const$getSensitive()) {
      V <- NULL
    } else {
      V <- eval(parse(text=vertSize.const$getText()))
    }
    vsize <- switch(vertSize$getActive()+1,
      mult*V,
      mult*V(g)$degree,
      mult*25*V(g)$ev.cent,
      mult*3*log1p(V(g)$btwn.cent)+.05,
      mult*range.transform(V(g)$subgraph.cent, 0, 15),
      mult*V(g)$coreness,
      mult*20*V(g)$transitivity,
      mult*range.transform(V(g)$PC, 0, 15),
      mult*range.transform(V(g)$l.eff, 0, 15),
      mult*range.transform(abs(V(g)$z.score), 0, 15),
      mult*10*sqrt(V(g)$hub.score)
    )

    # Edge width
    ewidth <- eval(parse(text=edgeWidth$getText()))

    # Community number, if applicable
    if (identical(plotFunc, plot.community)) {
      cNum <- comm$getActive() + 1
      plotFunc(g, n=cNum,
               vertex.label=vertex.label, vertex.label.cex=vertex.label.cex,
               vertex.size=vsize,
               edge.width=ewidth,
               vertex.color=vertex.color,
               edge.color=edge.color)
    }

    # Vertex neighborhood, if applicable
    if (identical(plotFunc, plot.neighborhood)) {
      n <- comboNeighb$getActive() + 1
      if (vertColor$getActive() + 1 == 1) {
        vertex.color <- V(g)$color <- rep('lightblue3', Nv)
        vertex.color[n] <- V(g)[n]$color <- 'yellow'
        edge.color <- E(g)$color <- 'red'
      }
      plotFunc(g, n=n,
               vertex.label=vertex.label,
               vertex.label.cex=vertex.label.cex,
               vertex.size=vsize,
               edge.width=ewidth,
               vertex.color=vertex.color,
               edge.color=edge.color)
    }

    if (identical(plotFunc, plot.adj)) {
      if (orient$getActive() == 3) {
        plotFunc <- plot
      }

      # Slider for curvature of edges in circle plots
      #if (!is.null(the.slider)) {
      if (length(class(the.slider)) > 1) {
        curv <- the.slider$getValue()
      } else {
        curv <- 0
      }

      plotFunc(g,
               vertex.label=vertex.label, vertex.label.cex=vertex.label.cex,
               vertex.size=vsize,
               edge.width=ewidth,
               vertex.color=vertex.color,
               edge.color=edge.color,
               xlim=xlim.g,
               ylim=ylim.g,
               layout=layout.g,
               edge.curved=curv)
    }

    # Show a legend for lobe colors
    if (vertColor$getActive() + 1 == 3) {
      lobes <- sort(unique(V(g)$lobe))
      lobe.names <- atlas.list$lobe.names[lobes]
      lobe.cols <- unique(V(g)$color.lobe[order(V(g)$lobe)])
      kNumLobes <- max(V(g)$lobe)
      if (orient$getActive() == 3) {
        legend.text.col <- 'black'
      } else {
        legend.text.col <- 'white'
      }
      legend('topleft',
             lobe.names,
             fill=lobe.cols,
             text.col=legend.text.col)
    }

    # Show the diameter of each graph?
    if (showDiameter$active == TRUE) {
      g.sub <- induced.subgraph(g, get.diameter(g))
      V(g.sub)$color <- 'deeppink'
      E(g.sub)$color <- 'deeppink'
      E(g.sub)$width <- 5
      V(g.sub)$size <- 10
      #browser()
      plotFunc(g.sub, add=T, vertex.label=NA,
               xlim=xlim.g,
               ylim=ylim.g)
    }

    # Show edge differences?
    if (edgeDiffs$active == TRUE) {
      g.diff12 <- graph.difference(g, g2)
      plot.adj(g.diff12, add=T,
               vertex.label=NA, vertex.size=degree(g.diff12),
               vertex.color='deeppink',edge.color='deeppink', 
               edge.width=5,
               layout=layout.g)
    }
  }

#browser()
  graph1 <- eval(parse(text=graphname1$getText()))
  graph2 <- eval(parse(text=graphname2$getText()))
  if (!is.igraph(graph1)) {
    stop(sprintf('%s is not a graph object.', graphname1$getText()))
  }
  
  make.plot(dev=firstplot, g=graph1, g2=graph2, orient=orient1,
            comm, vertLabels, vertSize, edgeWidth,
            vertColor, vertSize.const=NULL, hemi, the.slider=slider1)

  if (nchar(graphname2$getText()) > 0) {
    if (!is.igraph(graph2)) {
      stop(sprintf('%s is not a graph object.', graphname2$getText()))
    }
    make.plot(dev=secondplot, g=graph2, g2=graph1, orient=orient2,
              comm, vertLabels, vertSize, edgeWidth,
              vertColor, vertSize.const=NULL, hemi, the.slider=slider2)
  }
}
