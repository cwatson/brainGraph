#' GUI for plotting graphs overlaid on MNI axial image.
#'
#' This function creates a GUI for plotting the graph over an axial image from
#' the MNI template. It gives the user control over several plotting parameters.
#'
#' @export

plot.adj.gui <- function() {
  window <- gtkWindow('toplevel')
  window['title'] <- 'brainGraph'

  #=============================================================================
  # Add a menubar for different plotting functions
  #=============================================================================
  plotFunc <- plot.adj

  menubar <- gtkMenuBar()
  plot_menu <- gtkMenu()
  plot_item <- gtkMenuItemNewWithMnemonic(label='_Plot')
  plot_item$setSubmenu(plot_menu)
  menubar$append(plot_item)

  # Plot the entire graph (default)
  #---------------------------------------------------------
  adj_item <- gtkMenuItemNewWithMnemonic('_Entire graph')
  gSignalConnect(adj_item, 'activate', function(item) {
                 commNum$setSensitive(F)
                 neighbVert$setSensitive(F)
                 plotFunc <<- plot.adj
  })
  plot_menu$append(adj_item)

  # Plot just a neighborhood
  #---------------------------------------------------------
  neighb_item <- gtkMenuItemNewWithMnemonic('_Neighborhood graph')
  gSignalConnect(neighb_item, 'activate', function(item) {
                 commNum$setSensitive(F)
                 neighbVert$setSensitive(T)
                 comboVcolor$setActive(0)
                 plotFunc <<- plot.neighborhood
  })
  plot_menu$append(neighb_item)

  # Plot just a community
  #---------------------------------------------------------
  community_item <- gtkMenuItemNewWithMnemonic('_Community graph')
  gSignalConnect(community_item, 'activate', function(item) {
                 commNum$setSensitive(T)
                 neighbVert$setSensitive(F)
                 comboVcolor$setActive(0)
                 plotFunc <<- plot.community
  })
  plot_menu$append(community_item)

  vboxMainMenu <- gtkVBoxNew(F, 8)
  window$add(vboxMainMenu)
  vboxMainMenu$packStart(menubar, F, F)
  #=============================================================================


  # Create main (horizontal) container
  hboxMain <- gtkHBoxNew(F, 8)
  vboxMainMenu$add(hboxMain)

  #=============================================================================
  # Create main (left side) vertical container
  #=============================================================================
  vboxMain <- gtkVBoxNew(F, 8)
  hboxMain$add(vboxMain)

  #---------------------------------------------------------
  # Create frame + vert container for plot parameters
  #---------------------------------------------------------
  frame <- gtkFrameNew('Specify plotting parameters')
  vboxMain$add(frame)
  
  vbox <- gtkVBoxNew(F, 8)
  vbox$setBorderWidth(10)
  frame$add(vbox)

  # Create a frame for each plot area
  #---------------------------------------------------------
  frameG <- list(length=2)
  vboxG <- list(length=2)
  hboxOrient <- list(length=2)
  comboOrient <- list(length=2)
  graphname <- list(length=2)
  for (i in 1:2) {
    frameG[[i]] <- gtkFrameNew(sprintf('Graph %i', i))
    vbox$add(frameG[[i]])

    vboxG[[i]] <- gtkVBoxNew(F, 8)
    vboxG[[i]]$setBorderWidth(10)
    frameG[[i]]$add(vboxG[[i]])

    hbox <- gtkHBoxNew(F, 8)
    vboxG[[i]]$packStart(hbox, expand=F, fill=F, padding=0)

    label <- gtkLabelNewWithMnemonic('Name:')
    hbox$packStart(label, F, F, 0)
    graphname[[i]] <- gtkEntryNew()
    graphname[[i]]$setWidthChars(20)
    label$setMnemonicWidget(graphname[[i]])
    hbox$packStart(graphname[[i]], F, F, 0)

    hboxOrient[[i]] <- gtkHBoxNew(F, 8)
    vboxG[[i]]$add(hboxOrient[[i]])
    choices <- c('Axial', 'Sagittal (left)', 'Sagittal (right)', 'Circular')
    comboOrient[[i]] <- gtkComboBoxNewText()
    comboOrient[[i]]$show()
    for (choice in choices) comboOrient[[i]]$appendText(choice)
    comboOrient[[i]]$setActive(0)

    label <- gtkLabelNewWithMnemonic('Orientation')
    hboxOrient[[i]]$packStart(label, F, F, 0)
    hboxOrient[[i]]$add(comboOrient[[i]])
  }
  #---------------------------------------------------------


  # Add horizontal container to specify options
  # Should vertex labels be displayed?
  hbox <- gtkHBoxNew(F, 8)
  vbox$packStart(hbox, F, F, 0)
  label <- gtkLabelNewWithMnemonic('Display vertex labels?')
  hbox$packStart(label, F, F, 0)
  vertLabels <- gtkCheckButton()
  vertLabels$active <- FALSE
  hbox$packStart(vertLabels, F, F, 0)
  label$setMnemonicWidget(vertLabels)
  
  # Vertex colors based on community membership?
  #---------------------------------------
  hboxVcolor <- gtkHBoxNew(F, 8)
  vbox$packStart(hboxVcolor, F, F, 0)
  choices <- c('None (lightblue)', 'Communities', 'Lobes')
  comboVcolor <- gtkComboBoxNewText()
  comboVcolor$show()
  for (choice in choices) comboVcolor$appendText(choice)
  comboVcolor$setActive(1)

  label <- gtkLabelNewWithMnemonic('Vertex color')
  hboxVcolor$packStart(label, F, F, 0)
  hboxVcolor$add(comboVcolor)

  # Edge width?
  hbox <- gtkHBoxNew(F, 8)
  vbox$packStart(hbox, F, F, 0)
  label <- gtkLabelNewWithMnemonic('Edge width')
  hbox$packStart(label, F, F, 0)
  edgeWidth <- gtkEntryNew()
  edgeWidth$setWidthChars(3)
  edgeWidth$setText('1.5')
  hbox$packStart(edgeWidth, F, F, 0)
  label$setMnemonicWidget(edgeWidth)
   
  # Vertex size?
  #---------------------------------------
  hboxVsize <- gtkHBoxNew(F, 8)
  vbox$packStart(hboxVsize, F, F, 0)
  choices <- c('Constant', 'Degree', 'EV centrality', 'Bwtn centrality',
               'Subgraph centrality', 'Coreness',
               'Clustering coeff.', 'Part. coeff.', 'Loc. eff.',
               'Within-module degree z-score', 'Hub score')
  combo <- gtkComboBoxNewText()
  combo$show()
  for (choice in choices) combo$appendText(choice)
  combo$setActive(1)

  label <- gtkLabelNewWithMnemonic('Vertex size')
  hboxVsize$packStart(label, F, F, 0)
  hboxVsize$add(combo)

  vertSize.const <- gtkEntryNew()
  vertSize.const$setWidthChars(3)
  vertSize.const$setText('5')
  hboxVsize$packStart(vertSize.const, F, F, 0)
  vertSize.const$setSensitive(F)

  gSignalConnect(combo, 'changed', function(combo, ...) {
      if (combo$getActive() == 0) {
        vertSize.const$setSensitive(T)
      } else {
        vertSize.const$setSensitive(F)
      }
    })

  # Both, single hemisphere, or inter-hemispheric only?
  #---------------------------------------------------------
  hboxHemi <- gtkHBoxNew(F, 8)
  vbox$packStart(hboxHemi, F, F, 0)
  choices <- c('Both', 'Left only', 'Right only', 'Interhemispheric only')
  
  comboHemi <- gtkComboBoxNewText()
  comboHemi$show()
  for (choice in choices) comboHemi$appendText(choice)
  comboHemi$setActive(0)

  label <- gtkLabelNew('Hemisphere')
  hboxHemi$packStart(label, F, F, 0)
  hboxHemi$add(comboHemi)

  # Major lobe number, if applicable
  hboxLobe <- gtkHBoxNew(F, 8)
  vbox$packStart(hboxLobe, F, F, 0)
  choices <- c('All', 'Frontal', 'Parietal', 'Temporal', 'Occipital',
               'Frontal & Parietal', 'Frontal & Temporal', 'Frontal & Occipital',
               'Parietal & Temporal', 'Parietal & Occipital',
               'Temporal & Occipital')

  comboLobe <- gtkComboBoxNewText()
  comboLobe$show()
  for (choice in choices) comboLobe$appendText(choice)
  comboLobe$setActive(0)

  label <- gtkLabelNew('Lobe')
  hboxLobe$packStart(label, F, F, 0)
  hboxLobe$add(comboLobe)


  # Community number, if applicable
  hboxComm <- gtkHBoxNew(F, 8)
  vbox$packStart(hboxComm, F, F, 0)
  label <- gtkLabelNew('Which community? (ordered by size)')
  hboxComm$packStart(label, F, F, 0)
  commNum <- gtkEntryNew()
  commNum$setWidthChars(3)
  commNum$setText('1')
  hboxComm$packStart(commNum, F, F, 0)
  commNum$setSensitive(F)

  # Vertex and its neighborhood, if applicable
  hboxNeighb <- gtkHBoxNew(F, 8)
  vbox$packStart(hboxNeighb, F, F, 0)
  label <- gtkLabelNew('Which vertex?')
  hboxNeighb$packStart(label, F, F, 0)
  neighbVert <- gtkEntryNew()
  neighbVert$setWidthChars(4)
  neighbVert$setText('1')
  hboxNeighb$packStart(neighbVert, F, F, 0)
  neighbVert$setSensitive(F)
   
  #-----------------------------------------------------------------------------
  # Add two horizontal containers to check if the results will be saved to a file
  # and if so, to specify the file's name
  #-----------------------------------------------------------------------------
  vbox <- gtkVBoxNew(F, 8)
  vboxMain$add(vbox)
  hbox <- gtkHBoxNew(F, 8)
  vbox$packStart(hbox, F, F, 0)
  label <- gtkLabelNewWithMnemonic('Save figure 1?')
  hbox$packStart(label, F, F, 0)
  toSave1 <- gtkCheckButton()
  hbox$packStart(toSave1, F, F, 0)
  label$setMnemonicWidget(toSave1)
  label <- gtkLabelNewWithMnemonic('File name:')
  hbox$packStart(label, F, F, 0)
  exportFileName1 <- gtkEntryNew()
  exportFileName1$setWidthChars(20)
  exportFileName1$setText('output1')
  hbox$packStart(exportFileName1, F, F, 0)
  label$setMnemonicWidget(exportFileName1)
  fileExt <- gtkLabel('.png')
  hbox$packStart(fileExt, F, F, 0)
  
  # Save option for the second plotting area
  hbox <- gtkHBoxNew(F, 8)
  vbox$packStart(hbox, F, F, 0)
  label <- gtkLabelNewWithMnemonic('Save figure 2?')
  hbox$packStart(label, F, F, 0)
  toSave2 <- gtkCheckButton()
  hbox$packStart(toSave2, F, F, 0)
  label$setMnemonicWidget(toSave2)
  label <- gtkLabelNewWithMnemonic('File name:')
  hbox$packStart(label, F, F, 0)
  exportFileName2 <- gtkEntryNew()
  exportFileName2$setWidthChars(20)
  exportFileName2$setText('output2')
  hbox$packStart(exportFileName2, F, F, 0)
  label$setMnemonicWidget(exportFileName2)
  fileExt <- gtkLabel('.png')
  hbox$packStart(fileExt, F, F, 0)

  #-------------------------------------
  # Create 2 drawing areas for the plotting
  #-------------------------------------
  graphics <- vector('list', 2)
  vboxPlot <- vector('list', 2)
  groupplot <- vector('list', 2)

  for (i in 1:2) {
    graphics[[i]] <- gtkDrawingArea()
    graphics[[i]]$setSizeRequest(640, 640)

    vboxPlot[[i]] <- gtkVBox()
    vboxPlot[[i]]$packStart(graphics[[i]], expand=T, fill=T, padding=0)
    hboxMain$add(vboxPlot[[i]])

    asCairoDevice(graphics[[i]])
    par(pty='s', mar=rep(0, 4))
    groupplot[[i]] <- dev.cur()
  }

  #-----------------------------------------------------------------------------
  # Function to dynamically draw the graph
  #-----------------------------------------------------------------------------
  update.adj <- function(graphname1, graphname2, vertLabels, vertSize,
                         edgeWidth, vertColor, toSave1=FALSE, toSave2=FALSE,
                         exportFileName1, exportFileName2, fileExt, firstplot,
                         secondplot, vertSize.const=NULL, plotFunc, commNum,
                         neighbVert, hemi, lobe, orient1, orient2) {
    graph1 <- eval(parse(text=graphname1$getText()))
    g2 <- eval(parse(text=graphname2$getText()))


    #===========================================================================
    #===========================================================================
    # Function that does all the work
    #===========================================================================
    #===========================================================================
    make.plot <- function(dev, g, toSave, exportFileName, orient, ...) {
      dev.set(dev)

      # Lobe to plot
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
        if (g$atlas == 'aal90' || g$atlas == 'lpba40' || g$atlas == 'hoa112') {
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
        if (g$atlas == 'aal90' || g$atlas == 'lpba40' || g$atlas == 'hoa112') {
          g <- induced.subgraph(g, seq(1, Nv, 2))
          layout.g <- matrix(c(-eval(parse(text=g$atlas))$brainnet.coords[V(g)$name, 2],
                               eval(parse(text=g$atlas))$brainnet.coords[V(g)$name, 3]),
                             ncol=2, byrow=F)
        } else {
          g <- induced.subgraph(g, 1:(Nv/2))
          layout.g <- matrix(c(-eval(parse(text=g$atlas))$brainnet.coords[V(g)$name, 2],
                               eval(parse(text=g$atlas))$brainnet.coords[V(g)$name, 3]),
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
        if (g$atlas == 'aal90' || g$atlas == 'lpba40' || g$atlas == 'hoa112') {
          g <- induced.subgraph(g, seq(2, Nv, 2))
          layout.g <- matrix(c(eval(parse(text=g$atlas))$brainnet.coords[V(g)$name, 2],
                               eval(parse(text=g$atlas))$brainnet.coords[V(g)$name, 3]),
                             ncol=2, byrow=F)
        } else {
          g <- induced.subgraph(g, ((Nv/2 + 1):Nv))
          layout.g <- matrix(c(eval(parse(text=g$atlas))$brainnet.coords[V(g)$name, 2],
                               eval(parse(text=g$atlas))$brainnet.coords[V(g)$name, 3]),
                             ncol=2, byrow=F)
        }
        xlim.g <- c(-125, 85)
        ylim.g <- c(-85, 125)
        mult <- 100

      #---------------------------------
      # Plot in a circular layout
      #---------------------------------
      } else if (orient$getActive() == 3) {
        # Hemisphere to plot
        circ <- V(g)$circle.layout
        lobe <- V(g)$lobe
        lobe.color <- V(g)$lobe.color
        atlas <- g$atlas
        if (hemi$getActive()+1 == 1) {
          sg <- g
        } else if (hemi$getActive()+1 == 2) {
          sg <- g - E(g) + subgraph.edges(g, E(g)[1:(Nv/2) %--% 1:(Nv/2)])
        } else if (hemi$getActive()+1 == 3) {
          sg <- g - E(g) + subgraph.edges(g, E(g)[(Nv/2 + 1):Nv %--% (Nv/2 + 1):Nv])
        } else if (hemi$getActive()+1 == 4) {
          sg <- g - E(g) + subgraph.edges(g, E(g)[1:(Nv/2) %--% (Nv/2 + 1):Nv])
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
      # Vertex colors
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
      if (commNum$getSensitive()) {
        cNum <- eval(parse(text=commNum$getText()))
        plotFunc(g, n=cNum,
                 vertex.label=vertex.label, vertex.label.cex=vertex.label.cex,
                 vertex.size=vsize,
                 edge.width=ewidth,
                 vertex.color=vertex.color,
                 edge.color=edge.color)
      }

      # Vertex neighborhood, if applicable
      if (neighbVert$getSensitive()) {
        v <- V(g)$name[eval(parse(text=neighbVert$getText()))]
        plotFunc(g, v=v,
                 vertex.label=vertex.label,
                 vertex.label.cex=vertex.label.cex,
                 vertex.size=vsize,
                 edge.width=ewidth,
                 vertex.color=vertex.color,
                 edge.color=edge.color)
      }

      if ((!commNum$getSensitive()) && (!neighbVert$getSensitive())) {
        if (orient$getActive() == 3) {
          plotFunc <- plot
        }
        plotFunc(g, vertex.label=vertex.label, vertex.label.cex=vertex.label.cex,
               vertex.size=vsize,
               edge.width=ewidth,
               vertex.color=vertex.color,
               edge.color=edge.color,
               xlim=xlim.g,
               ylim=ylim.g,
               layout=layout.g)
        if (orient$getActive() == 3) {
          lobe.cols <- unique(V(g)$lobe.color[order(V(g)$lobe)])
          kNumLobes <- max(V(g)$lobe)
          legend('topleft', eval(parse(text=g$atlas))$lobe.names[1:kNumLobes],
                 fill=lobe.cols[1:kNumLobes])
        }
      }

      if (toSave$active == TRUE) {
        fname <- paste0(exportFileName$getText(), fileExt$getText())
        dev.copy(png, filename=fname)
        dev.off()
      }
    }

    make.plot(dev=firstplot, g=graph1, toSave=toSave1,
              exportFileName=exportFileName1, vertLabels, vertSize, edgeWidth,
              vertColor, fileExt, vertSize.const=NULL, commNum, hemi,
              orient=orient1)
    if (nchar(graphname2$getText()) > 0) {
      make.plot(dev=secondplot, g=g2, vertLabels, vertSize, edgeWidth, vertColor,
                toSave=toSave2, exportFileName=exportFileName2, fileExt,
                vertSize.const=NULL, commNum, hemi, orient=orient2)
    }
  }
  #-----------------------------------------------------------------------------

  # Add buttons
  the.buttons <- gtkHButtonBoxNew()
  the.buttons$setBorderWidth(5)
  vboxMain$add(the.buttons)
  the.buttons$setLayout("spread")
  the.buttons$setSpacing(40)
  buttonOK <- gtkButtonNewFromStock("gtk-ok")
  gSignalConnect(buttonOK, "clicked",
                 function(widget) update.adj(graphname[[1]], graphname[[2]], vertLabels,
                                   vertSize=combo, edgeWidth,
                                   vertColor=comboVcolor, toSave1,
                                   toSave2, exportFileName1, exportFileName2,
                                   fileExt, firstplot=groupplot[[1]],
                                   secondplot=groupplot[[2]],
                                   vertSize.const, plotFunc, commNum,
                                   neighbVert, hemi=comboHemi, lobe=comboLobe,
                                   orient1=comboOrient[[1]],
                                   orient2=comboOrient[[2]]))
  the.buttons$packStart(buttonOK, fill=F)
  buttonClose <- gtkButtonNewFromStock("gtk-close")
  gSignalConnect(buttonClose, "clicked", window$destroy)
  the.buttons$packStart(buttonClose, fill=F)

  hboxMain$showAll()

}
