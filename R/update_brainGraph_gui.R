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
#' @param vertSize.const A GTK entry for constant vertex size
#' @param edgeWidth.const A GTK entry for constant width
#' @param vertLabels A GTK check button for showing vertex labels
#' @param comm A GTK combo box for plotting individual communities
#' @param kNumComms Integer indicating the number of total communities (optional)
#' @param neighb A GTK combo box for plotting individual neighborhoods
#' @param neighbMult A GTK entry for joint neighborhoods of multiple vertices
#' @param slider A GTK horizontal slider widget for changing edge curvature
#' @param vertSize.other A GTK entry for vertex size (other attributes)
#' @param vertSize.eqn A GTK entry for equations to exclude vertices
#' @param showDiameter A GTK check button for showing the graph's diameter
#' @param edgeDiffs A GTK check button for showing edge diffs between graphs
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

update_brainGraph_gui <- function(plotDev, graph1, graph2, plotFunc, vertSize,
                       edgeWidth, vertColor, hemi, lobe, orient, vertSize.min,
                       edgeWidth.min, vertSize.const=NULL, edgeWidth.const=NULL,
                       vertLabels=NULL, comm=NULL, kNumComms=NULL, neighb=NULL,
                       neighbMult=NULL, slider=NULL, vertSize.other=NULL,
                       vertSize.eqn=NULL, showDiameter=NULL, edgeDiffs=NULL) {

  #===========================================================================
  #===========================================================================
  # Function that does all the work
  #===========================================================================
  #===========================================================================
  make.plot <- function(dev, g, g2, orient, vertSize.min, edgeWidth.min,
                        vertSize.const=NULL, vertSize.eqn=NULL, the.slider=NULL,
                        kNumComms=NULL, comm=NULL, ...) {
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


    # Vertex neighborhoods, if applicable
    if (!is.function(plotFunc) && plotFunc == 'plot_neighborhood') {
      n <- neighb$getActive()
      if (n == 0) {
        splits <- strsplit(neighbMult, split=', ')[[1]]
        if (any(grep('^[[:digit:]]*$', splits))) {  # numeric
          verts <- as.numeric(splits)
          vnames <- V(g)[verts]$name
        } else {
          verts <- c(splits)
          vnames <- verts
        }
      } else {
        verts <- n
        vnames <- V(g)[verts]$name
      }
        g.sub <- graph_neighborhood_multiple(g, verts)
        g.sub <- set.brainGraph.attributes(g.sub, atlas)
        g <- g.sub
        verts <- which(V(g)$name %in% vnames)
        Nv <- vcount(g)
        plotFunc <- plot_brainGraph
    } else {
      n <- 1
      verts <- NULL
    }

    # Community number, if applicable
    if (!is.function(plotFunc) && plotFunc == 'plot_community') {
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
      cNums <- combos[[ind1]][, ind2]
      comms <- as.numeric(names(rev(sort(table(V(g)$comm)))[cNums]))
      memb <- which(V(g)$comm %in% comms)
      g.sub <- induced.subgraph(g, memb)
      lobe.cols <- c('red', 'green', 'blue', 'magenta', 'yellow', 'orange',
                     'lightgreen', 'lightblue', 'lightyellow')
      vcomms <- V(g.sub)$comm
      eids <- sapply(seq_along(comms),
        function(x) as.numeric(E(g.sub)[which(vcomms == comms[x]) %--% which(vcomms == comms[x])]))
      if (length(cNums) == 1) {
        ecols <- rep(lobe.cols[cNums], length=ecount(g.sub))
      } else {
        listcols <- as.list(lobe.cols[cNums])
        ecols <- rep('gray50', length=ecount(g.sub))
        for (i in seq_along(eids)) {
          ecols[eids[[i]]] <- listcols[[i]]
        }
      }
      g <- g.sub
      E(g)$color.comm <- ecols
      plotFunc <- plot_brainGraph
    }

    #-----------------------------------
    # Vertex sizes
    #-----------------------------------
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
    v.min <- vertSize.min$getValue()

    i <- vertSize$getActive()
    vsize.opts <- c('const', 'degree', 'ev.cent', 'btwn.cent', 'subgraph.cent',
                    'coreness', 'transitivity', 'PC', 'E.local', 'E.nodal',
                    'z.score', 'hub.score', 'vulnerability', 'knn')
    if (i == 0) {
      vsize <- mult * const
    } else {
      if (i < 14) {
        if (i == 11 && !is.directed(g)) {
          i <- 2
        }
        vnum <- vertex_attr(g, vsize.opts[i + 1])
        removed <- union(which(vnum < v.min), which(!is.finite(vnum)))
        g <- delete.vertices(g, removed)
        if (vcount(g) %in% c(0, 1)) {
          vsize <- vnum[-removed]
        } else {
          vsize <- mult * vec.transform(vnum, 0, 15)
          if (length(removed) > 0) vsize <- vsize[-removed]
          vsize <- vsize[!is.nan(vsize)]
        }

      } else if (vertSize$getActive() == 14) {
        g <- delete.vertices(g, which(vertex_attr(g, v.attr) < v.min))
        if (v.attr %in% c('p', 'p.adj', 'p.perm')) {
          vsize <- mult * 15 * vertex_attr(g, v.attr)
        } else {
          vsize <- mult * vertex_attr(g, v.attr)
        }
      } else if (vertSize$getActive() == 15) {  # equation
        x <- vertSize.eqn$getText()
        if (nchar(x) > 0) {
          subs <- strsplit(x, split='&')[[1]]
          subs <- gsub('^\\s+|\\s+$', '', subs)
          cond <- eval(parse(text=paste0('V(g)$', subs, collapse='&')))
          g <- delete.vertices(g, setdiff(seq_len(vcount(g)), which(cond)))
        }
        vsize <- mult * 7.5
      }
    }

    # Edge width
    if (!edgeWidth.const$getSensitive()) {
      e.const <- NULL
    } else {
      e.const <- eval(parse(text=edgeWidth.const$getText()))
    }

    if (edgeWidth$getActive() == 0) {
      ewidth <- e.const
    } else if (edgeWidth$getActive() == 1) {
      e.min <- edgeWidth.min$getValue()
      g <- delete.edges(g, which(E(g)$btwn < e.min))
      ewidth <- log1p(E(g)$btwn)
    } else if (edgeWidth$getActive() == 2) {
      e.min <- edgeWidth.min$getValue()
      g <- delete.edges(g, which(E(g)$dist < e.min))
      ewidth <- vec.transform(E(g)$dist, 0.1, 5)
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
  if (!is.igraph(graph1)) {
    stop(sprintf('%s is not a graph object.', deparse(substitute(graph1))))
  }

  make.plot(dev=plotDev, g=graph1, g2=graph2, orient, vertSize.min,
            edgeWidth.min, vertLabels, vertSize, edgeWidth,
            vertColor, vertSize.const=vertSize.const, vertSize.eqn=vertSize.eqn,
            hemi, the.slider=slider, kNumComms=kNumComms, comm=comm)
}
