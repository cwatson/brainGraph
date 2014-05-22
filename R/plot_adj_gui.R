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
  graphname <- gtkEntryNew()
  graphname$setWidthChars(50)
  graphname$setText('adj.group1[[N]]')
  label$setMnemonicWidget(graphname)
  hbox$packStart(graphname, F, F, 0)

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
#  vertSize <- gtkEntryNew()
#  vertSize$setWidthChars(2)
#  vertSize$setText('10')
  hbox$add(combo)
#  label$setMnemonicWidget(vertSize)
   
  
   
  # Add two horizontal containers to check if the results will be saved to a file
  # and if so, to specify the file's name
  vbox <- gtkVBoxNew(F, 8)
  vboxMain$add(vbox)
  hbox <- gtkHBoxNew(F,8)
  vbox$packStart(hbox, F, F, 0)
  label <- gtkLabelNewWithMnemonic('Save figure?')
  hbox$packStart(label,F,F,0)
  toSave <- gtkCheckButton()
  hbox$packStart(toSave,F,F,0)
  label$setMnemonicWidget(toSave)
  label <- gtkLabelNewWithMnemonic('File name:')
  hbox$packStart(label,F,F,0)
  exportFileName <- gtkEntryNew()
  exportFileName$setWidthChars(50)
  exportFileName$setText('output')
  hbox$packStart(exportFileName,F,F,0)
  label$setMnemonicWidget(exportFileName)
  label <- gtkLabel('.png')
  hbox$packStart(label,F,F,0)
  
  # Create a drawing area for the plot
  graphics <- gtkDrawingArea()
  graphics$setSizeRequest(640, 640)

  #-----------------------------------------------------------------------------
  # Function to dynamically draw the graph
  #-----------------------------------------------------------------------------
  update.adj <- function(graphname, vertLabels, vertSize, edgeWidth, vertColor) {
    n <- nchar(graphname$getText())
    N <- eval(parse(text=substr(graphname$getText(), n-2, n-2)))
    g <- eval(parse(text=graphname$getText()))

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
  }
  #-----------------------------------------------------------------------------

  vboxPlot <- gtkVBox()
  vboxPlot$packStart(graphics, expand=T, fill=T, padding=0)
  hboxMain$add(vboxPlot)


  # Add buttons
  the.buttons <- gtkHButtonBoxNew()
  the.buttons$setBorderWidth(5)
  vbox$add(the.buttons)
  the.buttons$setLayout("spread")
  the.buttons$setSpacing(40)
  buttonOK <- gtkButtonNewFromStock("gtk-ok")
  gSignalConnect(buttonOK, "clicked",
                 function(widget) update.adj(graphname, vertLabels,
                                   vertSize=combo, edgeWidth, vertColor))
  the.buttons$packStart(buttonOK,fill=F)
  buttonClose <- gtkButtonNewFromStock("gtk-close")
  gSignalConnect(buttonClose, "clicked", window$destroy)
  the.buttons$packStart(buttonClose,fill=F)

  hboxMain$showAll()
  asCairoDevice(graphics)
  
  par(pty='s', mar=rep(0, 4))
  plot.over.brain.axial(0)
}
