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
#'
#' @export

update.adj <- function(graphname1, graphname2, vertLabels, vertSize,
                       edgeWidth, edgeDiffs, vertColor, firstplot, secondplot,
                       vertSize.const=NULL, plotFunc, comm=NULL, comboNeighb,
                       hemi, lobe, orient1, orient2, showDiameter,
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
    Nv <- vcount(g)

    #=================================
    # Lobe to plot
    #=================================
    if (lobe$getActive() > 0) {
      kNumLobes <- length(atlas.list$lobe.names)
      combos <- sapply(seq_len(kNumLobes - 1),
                       function(x) combn(seq_along(atlas.list$lobe.names), x))
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

    #====================================================
    # Hemisphere to plot
    #====================================================
    # Both hemispheres
    if (hemi$getActive() == 0) {
      sg.hemi <- g
      memb.hemi <- seq_len(Nv)

    } else if (hemi$getActive() == 1) {
      # LH only
      if (atlas == 'aal90' || atlas == 'lpba40' || atlas == 'hoa112') {
        memb.hemi <- seq(1, Nv, 2)
      } else {
        memb.hemi <- 1:(Nv/2)
      }
      sg.hemi <- induced.subgraph(g, memb.hemi)

    } else if (hemi$getActive() == 2) {
      # RH only
      if (atlas == 'aal90' || atlas == 'lpba40' || atlas == 'hoa112') {
        memb.hemi <- seq(2, Nv, 2)
      } else {
        memb.hemi <- (Nv/2  + 1):Nv
      }
      sg.hemi <- induced.subgraph(g, memb.hemi)

    } else if (hemi$getActive() == 3) {
      # Interhemispheric only
      if (atlas == 'aal90' || atlas == 'lpba40' || atlas == 'hoa112') {
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
      layout.g <- cbind(V(g)$x, V(g)$y)
      xlim.g <- c(-1, 1)
      ylim.g <- c(-1.5, 1.5)
      mult <- 1

    #---------------------------------
    # Plot the left sagittal only
    #---------------------------------
    } else if (orient$getActive() == 1) {
      plot.over.brain.sagittal(0, hemi='left', z=30)
      cur.vertices <- intersect(which(memb.lobe), memb.hemi)
      layout.g <- matrix(c(-atlas.list$brainnet.coords[cur.vertices, 2],
                           atlas.list$brainnet.coords[cur.vertices, 3]),
                         ncol=2, byrow=F)
      xlim.g <- c(-85, 110)
      ylim.g <- c(-85, 125)
      mult <- 100

    #---------------------------------
    # Plot the right sagittal only
    #---------------------------------
    } else if (orient$getActive() == 2) {
      plot.over.brain.sagittal(0, hemi='right', z=30)
      cur.vertices <- intersect(which(memb.lobe), memb.hemi)
      layout.g <- matrix(c(atlas.list$brainnet.coords[cur.vertices, 2],
                           atlas.list$brainnet.coords[cur.vertices, 3]),
                         ncol=2, byrow=F)
      xlim.g <- c(-125, 85)
      ylim.g <- c(-85, 125)
      mult <- 100

    #---------------------------------
    # Plot in a circular layout
    #---------------------------------
    } else if (orient$getActive() == 3) {
      par(bg='white')
      circ <- V(g)$circle.layout

      layout.g <- rotation(layout.circle(g, order=circ), -pi/2)
      xlim.g <- c(-1.25, 1.25)
      ylim.g <- c(-1.25, 1.25)
      mult <- 1
    }
    #=======================================================
    #=======================================================
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
      mult*ifelse(V(g)$z.score == 0, 0, range.transform(V(g)$z.score, 0, 15)),
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
      if (vertColor$getActive() == 0) {
        vertex.color <- V(g)$color <- rep('lightblue3', Nv)
        vertex.color[n] <- V(g)[n]$color <- 'yellow'
        edge.color <- E(g)$color <- 'red'
      }
      plotFunc(g, n=n,
               vertex.label=vertex.label, vertex.label.cex=vertex.label.cex,
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
               vertex.label=vertex.label, vertex.label.cex=vertex.label.cex,
               vertex.size=vsize, edge.width=ewidth,
               vertex.color=vertex.color, edge.color=edge.color,
               xlim=xlim.g, ylim=ylim.g,
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
        legend.bg <- 'white'
        legend.text.col <- 'black'
      } else {
        legend.bg <- 'black'
        legend.text.col <- 'white'
      }
      legend('topleft',
             lobe.names,
             fill=lobe.cols,
             bg=legend.bg,
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
  #=============================================================================
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
