#' Function to dynamically plot a graph
#'
#' This function is called by \link{plot.adj.gui} to update a plot on-the-fly.
#' It updates both of the plots by calling the helper function "make.plot".

update.adj <- function(graphname1, graphname2, vertLabels, vertSize,
                       edgeWidth, edgeDiffs, vertColor, firstplot, secondplot,
                       vertSize.const=NULL, plotFunc, comm=NULL, comboNeighb, hemi,
                       lobe, orient1, orient2, showDiameter) {
  #===========================================================================
  #===========================================================================
  # Function that does all the work
  #===========================================================================
  #===========================================================================
  make.plot <- function(dev, g, g2, orient, comm=NULL, ...) {
    dev.set(dev)
    atlas <- g$atlas
    atlas.list <- eval(parse(text=g$atlas))

    #=================================
    # Lobe to plot
    #=================================
    g <- switch(lobe$getActive()+1,
                g,
                induced.subgraph(g, V(g)$lobe==1),
                induced.subgraph(g, V(g)$lobe==2),
                induced.subgraph(g, V(g)$lobe==3),
                induced.subgraph(g, V(g)$lobe==4),
                induced.subgraph(g, V(g)$lobe==1 | V(g)$lobe==2),
                induced.subgraph(g, V(g)$lobe==1 | V(g)$lobe==3),
                induced.subgraph(g, V(g)$lobe==1 | V(g)$lobe==4),
                induced.subgraph(g, V(g)$lobe==2 | V(g)$lobe==3),
                induced.subgraph(g, V(g)$lobe==2 | V(g)$lobe==4),
                induced.subgraph(g, V(g)$lobe==3 | V(g)$lobe==4)
    )               
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
      lobe.color <- V(g)$lobe.color

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
                           V(g)$color,
                           V(g)$lobe.color
    )
    edge.color <- switch(vertColor$getActive() + 1,
                         'red',
                         E(g)$color,
                         E(g)$lobe.color
    )

    # Vertex sizes
    if (!vertSize.const$getSensitive()) {
      V <- NULL
    } else {
      V <- eval(parse(text=vertSize.const$getText()))
    }
    vsize <- switch(vertSize$getActive()+1,
      mult*V,
      #mult*range.transform(V(g)$degree, 0, 15),
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

    if (is.null(comm) && is.null(comboNeighb)) {
      if (orient$getActive() == 3) {
        plotFunc <- plot
      }
      plotFunc(g,
               vertex.label=vertex.label, vertex.label.cex=vertex.label.cex,
               vertex.size=vsize,
               edge.width=ewidth,
               vertex.color=vertex.color,
               edge.color=edge.color,
               xlim=xlim.g,
               ylim=ylim.g,
               layout=layout.g)
    }

    # Show a legend for lobe colors
    if (vertColor$getActive() + 1 == 3) {
      lobe.cols <- unique(V(g)$lobe.color[order(V(g)$lobe)])
      kNumLobes <- max(V(g)$lobe)
      if (orient$getActive() == 3) {
        legend.text.col <- 'black'
      } else {
        legend.text.col <- 'white'
      }
      legend('topleft',
             atlas.list$lobe.names[1:kNumLobes],
             fill=lobe.cols[1:kNumLobes],
             text.col=legend.text.col)
    }

    # Show the diameter of each graph?
    if (showDiameter$active == TRUE) {
      d <- get.diameter(g)
      E(g, path=d)$color <- 'deeppink'
      vertex.color[d] <- 'deeppink'
      E(g, path=d)$width <- 5
      vsize[d] <- 10
      g.sub <- induced.subgraph(g, d)
      plot.adj(g.sub, add=T,
               vertex.label=NA, vertex.size=vsize[d],
               vertex.color=vertex.color[d])
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

  graph1 <- eval(parse(text=graphname1$getText()))
  graph2 <- eval(parse(text=graphname2$getText()))
  if (!is.igraph(graph1)) {
    stop(sprintf('%s is not a graph object.', graphname1$getText()))
  }
  
  make.plot(dev=firstplot, g=graph1, g2=graph2, orient=orient1,
            comm, vertLabels, vertSize, edgeWidth,
            vertColor, vertSize.const=NULL, hemi)
  if (nchar(graphname2$getText()) > 0) {
    if (!is.igraph(graph2)) {
      stop(sprintf('%s is not a graph object.', graphname2$getText()))
    }
    make.plot(dev=secondplot, g=graph2, g2=graph1, orient=orient2,
              comm, vertLabels, vertSize, edgeWidth,
              vertColor, vertSize.const=NULL, hemi)
  }
}
