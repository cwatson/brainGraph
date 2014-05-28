#' GUI for plotting graphs overlaid on MNI axial image.
#'
#' This function creates a GUI for plotting the graph over an axial image from
#' the MNI template. It gives the user control over several plotting parameters.
#'
#' @export

plot.adj.gui <- function() {
  window <- gtkWindow('toplevel')
  window['title'] <- 'brainGraph plotting'

  # Create main (horizontal) container
  hboxMain <- gtkHBoxNew(F, 8)
  window$add(hboxMain)

  # Create main (left side) vertical container
  vboxMain <- gtkVBoxNew(F, 8)
  hboxMain$add(vboxMain)

  # Create a frame and vertical container for plotting parameter specification
  frame <- gtkFrameNew('Specify plotting parameters')
  vboxMain$add(frame)
  
  vbox <- gtkVBoxNew(F, 8)
  vbox$setBorderWidth(24)
  frame$add(vbox)

  # Add horizontal container for every widget line
  hbox <- gtkHBoxNew(F, 8)
  vbox$packStart(hbox, expand=F, fill=F, padding=0)
  
  label <- gtkLabelNewWithMnemonic('Adjacency graph name')
  hbox$packStart(label, F, F, 0)
  # Add entry in the second column; named "graphname"
  graphname1 <- gtkEntryNew()
  graphname1$setWidthChars(50)
  graphname1$setText('adj.group1[[N]]')
  label$setMnemonicWidget(graphname1)
  hbox$packStart(graphname1, F, F, 0)

  hbox <- gtkHBoxNew(F, 8)
  vbox$packStart(hbox, expand=F, fill=F, padding=0)
  label <- gtkLabelNewWithMnemonic('Adjacency graph name')
  hbox$packStart(label, F, F, 0)
  # Add entry in the second column; named "graphname"
  graphname2 <- gtkEntryNew()
  graphname2$setWidthChars(50)
  graphname2$setText('adj.group2[[N]]')
  label$setMnemonicWidget(graphname2)
  hbox$packStart(graphname2, F, F, 0)

  # Add horizontal container to specify options
  # Should vertex labels be displayed?
  hbox <- gtkHBoxNew(F, 8)
  vbox$packStart(hbox, F, F, 0)
  label <- gtkLabelNewWithMnemonic('Display vertex labels?')
  hbox$packStart(label, F, F, 0)
  vertLabels <- gtkCheckButton()
  vertLabels$active <- TRUE
  hbox$packStart(vertLabels, F, F, 0)
  label$setMnemonicWidget(vertLabels)
  
  # Vertex colors based on community membership?
  hbox <- gtkHBoxNew(F, 8)
  vbox$packStart(hbox, F, F, 0)
  label <- gtkLabelNewWithMnemonic('Vertex community color?')
  hbox$packStart(label, F, F, 0)
  vertColor <- gtkCheckButton()
  hbox$packStart(vertColor, F, F, 0)
  label$setMnemonicWidget(vertColor)

  # Edge width?
  hbox <- gtkHBoxNew(F, 8)
  vbox$packStart(hbox, F, F, 0)
  label <- gtkLabelNewWithMnemonic('Edge width')
  hbox$packStart(label, F, F, 0)
  edgeWidth <- gtkEntryNew()
  edgeWidth$setWidthChars(1)
  edgeWidth$setText('2')
  hbox$packStart(edgeWidth, F, F, 0)
  label$setMnemonicWidget(edgeWidth)
   
  # Vertex size?
  #---------------------------------------
  hbox <- gtkHBoxNew(F, 8)
  vbox$packStart(hbox, F, F, 0)
  choices <- c('Constant', 'Degree', 'EV centrality', 'Bwtn centrality',
               'Clustering coeff.', 'Part. coeff.', 'Loc. eff.',
               'Within-module degree z-score')
  combo <- gtkComboBoxNewText()
  combo$show()
  for (choice in choices) combo$appendText(choice)
  combo$setActive(0)

  label <- gtkLabelNewWithMnemonic('Vertex size')
  hbox$packStart(label, F, F, 0)
  hbox$add(combo)
   
  #-----------------------------------------------------------------------------
  # Add two horizontal containers to check if the results will be saved to a file
  # and if so, to specify the file's name
  #-----------------------------------------------------------------------------
  vbox <- gtkVBoxNew(F, 8)
  vboxMain$add(vbox)
  hbox <- gtkHBoxNew(F, 8)
  vbox$packStart(hbox, F, F, 0)
  label <- gtkLabelNewWithMnemonic('Save figure 1?')
  hbox$packStart(label,F,F,0)
  toSave1 <- gtkCheckButton()
  hbox$packStart(toSave1,F,F,0)
  label$setMnemonicWidget(toSave1)
  label <- gtkLabelNewWithMnemonic('File name:')
  hbox$packStart(label,F,F,0)
  exportFileName1 <- gtkEntryNew()
  exportFileName1$setWidthChars(50)
  exportFileName1$setText('output1')
  hbox$packStart(exportFileName1,F,F,0)
  label$setMnemonicWidget(exportFileName1)
  fileExt <- gtkLabel('.png')
  hbox$packStart(fileExt,F,F,0)
  
  # Save option for the second plotting area
  vbox <- gtkVBoxNew(F, 8)
  vboxMain$add(vbox)
  hbox <- gtkHBoxNew(F,8)
  vbox$packStart(hbox, F, F, 0)
  label <- gtkLabelNewWithMnemonic('Save figure 2?')
  hbox$packStart(label,F,F,0)
  toSave2 <- gtkCheckButton()
  hbox$packStart(toSave2,F,F,0)
  label$setMnemonicWidget(toSave2)
  label <- gtkLabelNewWithMnemonic('File name:')
  hbox$packStart(label,F,F,0)
  exportFileName2 <- gtkEntryNew()
  exportFileName2$setWidthChars(50)
  exportFileName2$setText('output2')
  hbox$packStart(exportFileName2,F,F,0)
  label$setMnemonicWidget(exportFileName2)
  fileExt <- gtkLabel('.png')
  hbox$packStart(fileExt,F,F,0)

  # Create 2 drawing areas for the plotting
  graphics1 <- gtkDrawingArea()
  graphics1$setSizeRequest(560, 560)

  graphics2 <- gtkDrawingArea()
  graphics2$setSizeRequest(560, 560)

  #-----------------------------------------------------------------------------
  # Function to dynamically draw the graph
  #-----------------------------------------------------------------------------
  update.adj <- function(graphname1, graphname2, vertLabels, vertSize,
                         edgeWidth, vertColor, toSave1=FALSE, toSave2=FALSE,
                         exportFileName1, exportFileName2, fileExt, group1plot,
                         group2plot) {
    g1 <- eval(parse(text=graphname1$getText()))
    g2 <- eval(parse(text=graphname2$getText()))


    make.plot <- function(dev, g, toSave, exportFileName, ...) {
      dev.set(dev)
      plot.over.brain.axial(0)
      par(pty='s', mar=rep(0, 4))
      # Show vertex labels?
      if (vertLabels$active == FALSE) {
        vertex.label <- NA
        vertex.label.cex <- NA
      } else {
        vertex.label <- V(g)$name
        vertex.label.cex <- 0.75
      }
      # Show community membership colors?
      if (vertColor$active == TRUE) {
        vertex.color <- V(g)$color
        edge.color <- E(g)$color
      } else {
        vertex.color <- 'lightblue3'
        edge.color='red'
      }
      vsize <- switch(vertSize$getActive()+1,
        10,
        range.transform(V(g)$degree, 0, 15),
        25*V(g)$ev.cent,
        3*log1p(V(g)$btwn.cent)+.05,
        20*V(g)$transitivity,
        range.transform(V(g)$PC, 2.5, 15),
        range.transform(V(g)$l.eff, 0, 15),
        range.transform(abs(V(g)$z), 0, 15))
      ewidth <- eval(parse(text=edgeWidth$getText()))
      plot.adj(g, vertex.label=vertex.label, vertex.label.cex=vertex.label.cex,
               vertex.size=vsize,
               edge.width=ewidth,
               vertex.color=vertex.color,
               edge.color=edge.color)
      if (toSave$active == TRUE) {
        fname <- paste0(exportFileName$getText(), fileExt$getText())
        dev.copy(png, filename=fname)
        dev.off()
      }
    }

    make.plot(dev=group1plot, g=g1, vertLabels, vertSize, edgeWidth, vertColor,
              toSave=toSave1, exportFileName=exportFileName1, fileExt)
    if (nchar(graphname2$getText()) > 0) {
      make.plot(dev=group2plot, g=g2, vertLabels, vertSize, edgeWidth, vertColor,
                toSave=toSave2, exportFileName=exportFileName2, fileExt)
    } else {
    }

  }
  #-----------------------------------------------------------------------------

  vboxPlot1 <- gtkVBox()
  vboxPlot1$packStart(graphics1, expand=T, fill=T, padding=0)
  hboxMain$add(vboxPlot1)

  vboxPlot2 <- gtkVBox()
  vboxPlot2$packStart(graphics2, expand=T, fill=T, padding=0)
  hboxMain$add(vboxPlot2)

  asCairoDevice(graphics1)
  par(pty='s', mar=rep(0, 4))
  plot.over.brain.axial(0)
  group1plot <- dev.cur()

  asCairoDevice(graphics2)
  par(pty='s', mar=rep(0, 4))
  plot.over.brain.axial(0)
  group2plot <- dev.cur()


  # Add buttons
  the.buttons <- gtkHButtonBoxNew()
  the.buttons$setBorderWidth(5)
  vbox$add(the.buttons)
  the.buttons$setLayout("spread")
  the.buttons$setSpacing(40)
  buttonOK <- gtkButtonNewFromStock("gtk-ok")
  gSignalConnect(buttonOK, "clicked",
                 function(widget) update.adj(graphname1, graphname2, vertLabels,
                                   vertSize=combo, edgeWidth, vertColor, toSave1,
                                   toSave2, exportFileName1, exportFileName2,
                                   fileExt, group1plot, group2plot))
  the.buttons$packStart(buttonOK,fill=F)
  buttonClose <- gtkButtonNewFromStock("gtk-close")
  gSignalConnect(buttonClose, "clicked", window$destroy)
  the.buttons$packStart(buttonClose,fill=F)

  hboxMain$showAll()
}
