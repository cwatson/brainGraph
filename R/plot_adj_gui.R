#' GUI for plotting graphs overlaid on an MNI152 image or in a circle.
#'
#' This function creates a GUI for plotting graphs over an image from the MNI152
#' template. It gives the user control over several plotting parameters. Also
#' possible is a circular plot (in addition to the axial and sagittal views).
#'
#' @export

plot.adj.gui <- function() {
  window <- gtkWindow('toplevel')
  window['title'] <- 'brainGraph'
  window['icon'] <- gdkPixbuf(filename=system.file('extdata',
                                                   'brainGraph_icon.png',
                                                   package='brainGraph'))

  #=============================================================================
  # Add a menubar
  #=============================================================================

  #---------------------------------------------------------
  # Callback functions for the "File" menu
  #---------------------------------------------------------
  save_cb1 <- function(widget, window, plot.dev=gui.params$plot1) {
    dialog <- gtkFileChooserDialog('Enter a name for the file', window, 'save',
                                   'gtk-cancel', GtkResponseType['cancel'],
                                   'gtk-save', GtkResponseType['accept'])
    gSignalConnect(dialog, 'response', function(dialog, response) {
                   if(response == GtkResponseType['accept']) {
                     fname <- dialog$getFilename()
                     dev.set(plot.dev)
                     dev.copy(png, fname)
                     dev.off()
                   }
                  dialog$destroy()
                  })
  }
  save_cb2 <- function(widget, window, plot.dev=gui.params$plot2) {
    dialog <- gtkFileChooserDialog('Enter a name for the file', window, 'save',
                                   'gtk-cancel', GtkResponseType['cancel'],
                                   'gtk-save', GtkResponseType['accept'])
    gSignalConnect(dialog, 'response', function(dialog, response) {
                   if(response == GtkResponseType['accept']) {
                     fname <- dialog$getFilename()
                     dev.set(plot.dev)
                     dev.copy(png, fname)
                     dev.off()
                   }
                 dialog$destroy()
                 })
  }
  quit_cb <- function(widget, window) window$destroy()

  # Create a "GtkAction" class
  file.actions <- list(
    list('FileMenu', NULL, '_File'),
    list('Save1', 'gtk-save', 'Save_1', '<control>1', 'Save plot 1?', save_cb1),
    list('Save2', 'gtk-save', 'Save_2', '<control>2', 'Save plot 2?', save_cb2),
    list('Quit', 'gtk-quit', '_Quit', '<control>Q', 'Quit the application',
         quit_cb)
  )
  action_group <- gtkActionGroup('brainGraphActions')
  action_group$addActions(file.actions, window)

  #---------------------------------------------------------
  # Callback functions for the "Plot" menu
  #---------------------------------------------------------
  plot_entire_cb <- function(widget, window) {
    plotFunc <<- plot.adj
    # Get number of children of vboxMainMenu and kill all but 'menubar'
    kNumChildren <- length(window[[1]]$getChildren())
    if (kNumChildren > 1) {
      for (i in 2:kNumChildren) window[[1]][[i]]$destroy()
    }
    hboxMain <- gtkHBoxNew(F, 8)
    window[[1]]$add(hboxMain)
    gui.params <<- build.gui(hboxMain)

    comboComm <<- NULL
    kNumComms <<- NULL
    comboNeighb <<- NULL
  }

  #-------------------------------------
  # Plot neighborhood of individual vertices
  #-------------------------------------
  plot_neighb_cb <- function(widget, window) {
    plotFunc <<- plot.neighborhood
    # Get number of children of vboxMainMenu and kill all but 'menubar'
    kNumChildren <- length(window[[1]]$getChildren())
    if (kNumChildren > 1) {
      for (i in 2:kNumChildren) window[[1]][[i]]$destroy()
    }
    hboxMain <- gtkHBoxNew(F, 8)
    window[[1]]$add(hboxMain)
    gui.params <<- build.gui(hboxMain)

    hboxNeighb <- gtkHBoxNew(F, 8)
    hboxMain[[1]][[1]][[1]]$packStart(hboxNeighb, F, F, 0)

    comboNeighb <<- gtkComboBoxNewText()
    comboNeighb$show()
    label <- gtkLabelNew('Which vertex?')
    hboxNeighb$packStart(label, F, F, 0)
    hboxNeighb$add(comboNeighb)

    graph1 <- eval(parse(text=gui.params$graphname[[1]]$getText()))
    choicesNeighb <- data.frame(seq_len(vcount(graph1)), V(graph1)$name)
    for (r in seq_len(nrow(choicesNeighb))) {
      comboNeighb$appendText(paste(sprintf('%02i', choicesNeighb[r, 1]),
                                   ': ', '\t',
                                   as.character(choicesNeighb[r, 2])))
    }
    comboNeighb$setActive(0)

    comboComm <<- NULL
    kNumComms <<- NULL
    gui.params$comboVcolor$setActive(0)
    gui.params$showDiameter$setSensitive(F)
    gui.params$edgeDiffs$setSensitive(F)
    gui.params$comboHemi$setSensitive(F)
    gui.params$comboLobe$setSensitive(F)
  }

  #-------------------------------------
  # Plot individual communities
  #-------------------------------------
  plot_comm_cb <- function(widget, window) {
    plotFunc <<- plot.community
    # Get number of children of vboxMainMenu and kill all but 'menubar'
    if (length(window[[1]]$getChildren()) > 1) {
      kNumChildren <- length(window[[1]]$getChildren())
      for (i in 2:kNumChildren) window[[1]][[i]]$destroy()
    }
    hboxMain <- gtkHBoxNew(F, 8)
    window[[1]]$add(hboxMain)
    gui.params <<- build.gui(hboxMain)

    hboxComm <- gtkHBoxNew(F, 8)
    hboxMain[[1]][[1]][[1]]$packStart(hboxComm, F, F, 0)
  
    comboComm <<- gtkComboBoxNewText()
    comboComm$show()
    label <- gtkLabelNew('Which community? (ordered by size)')
    hboxComm$packStart(label, F, F, 0)
    hboxComm$add(comboComm)

    # Check if the text in the boxes represent igraph objects
    if (!gui.params$graphname[[1]]$getText() == '') {
      graph1 <- eval(parse(text=gui.params$graphname[[1]]$getText()))
      c1 <- which(rev(sort(table(V(graph1)$comm))) > 2)
      if (!gui.params$graphname[[2]]$getText() == '') {
        graph2 <- eval(parse(text=gui.params$graphname[[2]]$getText()))
        c2 <- which(rev(sort(table(V(graph2)$comm))) > 2)
        all.comms <- union(c1, c2)
      } else {
        all.comms <- c1
      }
    }
    kNumComms <<- length(all.comms)
    combos <- sapply(seq_len(kNumComms), function(x)
                     combn(seq_len(kNumComms), x))
    comms <- sapply(2:(kNumComms), function(z)
                    apply(t(apply(t(combos[[z]]), 1, function(x)
                                  all.comms[x])), 1, paste, collapse=', '))
    comms <- do.call('c', comms)
    choices <- c(as.character(all.comms), comms)
    for (choice in choices) comboComm$appendText(choice)
    comboComm$setActive(0)
    comboComm$setSensitive(T)

    comboNeighb <<- NULL
    gui.params$comboVcolor$setActive(1)
    gui.params$showDiameter$setSensitive(F)
    gui.params$edgeDiffs$setSensitive(F)
    gui.params$comboHemi$setSensitive(F)
    gui.params$comboLobe$setSensitive(F)
  }

  # Action group for the "Plot" menu
  plot.actions <- list(
    list('PlotMenu', NULL, '_Plot'),
    list('PlotAll', NULL, '_Entire Graph', '<control>E', 'Plot entire graph',
         plot_entire_cb),
    list('PlotNeighb', NULL, '_Neighborhood Graph', '<control>N', 'Plot
         neighborhood of a single vertex', plot_neighb_cb),
    list('PlotComm', NULL, '_Community Graph', '<control>C', 'Plot a single
         community', plot_comm_cb)
  )
  action_group$addActions(plot.actions, window)

  #---------------------------------------------------------
  # Create a GtkUIManager instance
  #---------------------------------------------------------
  ui_manager <- gtkUIManager()
  ui_manager$insertActionGroup(action_group, 0)

  merge <- ui_manager$newMergeId()
  ui_manager$addUi(merge.id=merge, path='/', name='menubar', action=NULL,
                   type='menubar', top=F)
  ui_manager$addUi(merge, '/menubar', 'file', 'FileMenu', 'menu', FALSE)
  ui_manager$addUi(merge, '/menubar/file', 'save1', 'Save1', 'menuitem', FALSE)
  ui_manager$addUi(merge, '/menubar/file', 'save2', 'Save2', 'menuitem', FALSE)
  ui_manager$addUi(merge, '/menubar/file', 'sep', action=NULL, type=NULL, FALSE)
  ui_manager$addUi(merge, '/menubar/file', 'quit', 'Quit', 'menuitem', FALSE)

  ui_manager$addUi(merge, '/menubar', 'plot', 'PlotMenu', 'menu', FALSE)
  ui_manager$addUi(merge, '/menubar/plot', 'plotall', 'PlotAll', 'menuitem',
                   FALSE)
  ui_manager$addUi(merge, '/menubar/plot', 'plotneighb', 'PlotNeighb',
                   'menuitem', FALSE)
  ui_manager$addUi(merge, '/menubar/plot', 'plotcomm', 'PlotComm', 'menuitem',
                   FALSE)

  menubar <- ui_manager$getWidget('/menubar')
  window$addAccelGroup(ui_manager$getAccelGroup())

  vboxMainMenu <- gtkVBoxNew(F, 8)
  window$add(vboxMainMenu)
  vboxMainMenu$packStart(menubar, F, F)
  #=============================================================================
  #=============================================================================
  
  # Temporary window to get graph object names
  set.names <- function() {
    graphObj <<- vector('list', length=2)
    tempWindow <- gtkWindow()
    tempVbox <- gtkVBoxNew(F, 8)
    tempWindow$add(tempVbox)
    label <- vector('list', length=2)
    graphObjEntry <- vector('list', length=2)
    tempHbox <- vector('list', length=2)
    for (i in 1:2) {
      tempHbox[[i]] <- gtkHBoxNew(F, 8)
      tempVbox$packStart(tempHbox[[i]], F, F, 0)
      label[[i]] <- gtkLabelNew(paste0('Graph ', i, ' name:'))
      tempHbox[[i]]$packStart(label[[i]], F, F, 0)
      graphObjEntry[[i]] <- gtkEntryNew()
      graphObjEntry[[i]]$setWidthChars(20)
      tempHbox[[i]]$packStart(graphObjEntry[[i]], F, F, 0)
    }
    button <- gtkButtonNewFromStock('gtk-ok')
    gSignalConnect(button, 'clicked', function(widget) {
                   graphObj[[1]] <<- graphObjEntry[[1]]$getText()
                   graphObj[[2]] <<- graphObjEntry[[2]]$getText()
                   tempWindow$destroy()
                   })
    tempVbox$packStart(button, fill=F)
  }
  set.names()


  build.gui <- function(container) {
  #=============================================================================
  # Create main (left side) vertical container
  #=============================================================================
  vboxMain <- gtkVBoxNew(F, 8)
  container$add(vboxMain)

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
  frameG <- vector('list', length=2)
  vboxG <- vector('list', length=2)
  hboxOrient <- vector('list', length=2)
  comboOrient <- vector('list', length=2)
  graphname <- vector('list', length=2)
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
    graphname[[i]]$setText(graphObj[[i]])
    graphname[[i]]$setSensitive(F)
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

  # Should vertex labels be displayed?
  #---------------------------------------
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
  choices <- c('None (lightblue)', 'Communities', 'Lobes', 'Components')
  comboVcolor <- gtkComboBoxNewText()
  comboVcolor$show()
  for (choice in choices) comboVcolor$appendText(choice)
  comboVcolor$setActive(1)

  label <- gtkLabelNewWithMnemonic('Vertex color')
  hboxVcolor$packStart(label, F, F, 0)
  hboxVcolor$add(comboVcolor)
   
  #-----------------------------------------------------------------------------
  # Vertex size?
  hboxVsize <- gtkHBoxNew(F, 8)
  vbox$packStart(hboxVsize, F, F, 0)
  choices <- c('Constant', 'Degree', 'EV centrality', 'Bwtn centrality',
               'Subgraph centrality', 'Coreness',
               'Clustering coeff.', 'Part. coeff.', 'Loc. eff.',
               'Within-module degree z-score', 'Hub score', 'Other')
  comboVsize <- gtkComboBoxNewText()
  comboVsize$show()
  for (choice in choices) comboVsize$appendText(choice)
  comboVsize$setActive(1)

  label <- gtkLabelNewWithMnemonic('Vertex size')
  hboxVsize$packStart(label, F, F, 0)
  hboxVsize$add(comboVsize)

  vertSize.const <- gtkEntryNew()
  vertSize.const$setWidthChars(3)
  vertSize.const$setText('5')
  hboxVsize$packStart(vertSize.const, F, F, 0)
  vertSize.const$setSensitive(F)

  hboxVsizeOther <- gtkHBoxNew(F, 8)
  vbox$packStart(hboxVsizeOther, F, F, 0)
  label <- gtkLabelNewWithMnemonic(paste0('\t\t', 'Other:'))
  hboxVsizeOther$packStart(label, F, F, 0)
  vertSize.other <- gtkEntryNew()
  vertSize.other$setWidthChars(10)
  hboxVsizeOther$packStart(vertSize.other, F, F, 0)
  vertSize.other$setSensitive(F)

  gSignalConnect(comboVsize, 'changed', function(comboVsize, ...) {
      if (comboVsize$getActive() == 0) {
        vertSize.const$setSensitive(T)
        vertSize.other$setSensitive(F)
      } else if (comboVsize$getActive() == 11) {
        vertSize.const$setSensitive(F)
        vertSize.other$setSensitive(T)
      } else {
        vertSize.const$setSensitive(F)
        vertSize.other$setSensitive(F)
      }
    })
  #-----------------------------------------------------------------------------
  # Edge width?
  hboxEwidth <- gtkHBoxNew(F, 8)
  vbox$packStart(hboxEwidth, F, F, 0)
  choices <- c('Constant', 'Edge betweenness')
  comboEwidth <- gtkComboBoxNewText()
  comboEwidth$show()
  for (choice in choices) comboEwidth$appendText(choice)
  comboEwidth$setActive(1)

  label <- gtkLabelNewWithMnemonic('Edge width')
  hboxEwidth$packStart(label, F, F, 0)
  hboxEwidth$add(comboEwidth)

  edgeWidth.const <- gtkEntryNew()
  edgeWidth.const$setWidthChars(3)
  edgeWidth.const$setText('1')
  hboxEwidth$packStart(edgeWidth.const, F, F, 0)
  edgeWidth.const$setSensitive(F)

  gSignalConnect(comboEwidth, 'changed', function(widget, ...) {
      if (widget$getActive() == 0) {
        edgeWidth.const$setSensitive(T)
      } else {
        edgeWidth.const$setSensitive(F)
      }
    })
  #-----------------------------------------------------------------------------

  # Highlight the diameter of each graph?
  hbox <- gtkHBoxNew(F, 8)
  vbox$packStart(hbox, F, F, 0)
  label <- gtkLabelNewWithMnemonic('Show diameter?')
  hbox$packStart(label, F, F, 0)
  showDiameter <- gtkCheckButton()
  showDiameter$active <- FALSE
  hbox$packStart(showDiameter, F, F, 0)
  label$setMnemonicWidget(showDiameter)

  # Show edge set differences?
  label <- gtkLabelNewWithMnemonic('Show edge differences?')
  hbox$packStart(label, F, F, 0)
  edgeDiffs <- gtkCheckButton()
  edgeDiffs$active <- FALSE
  hbox$packStart(edgeDiffs, F, F, 0)
  label$setMnemonicWidget(edgeDiffs)

  #-----------------------------------------------------------------------------
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
  atlas <- eval(parse(text=graphname[[1]]$getText()))$atlas
  atlas.list <- eval(parse(text=atlas))
  hboxLobe <- gtkHBoxNew(F, 8)
  vbox$packStart(hboxLobe, F, F, 0)

  kNumLobes <- length(atlas.list$lobe.names)
  combos <- sapply(seq_len(kNumLobes - 1),
                   function(x) combn(seq_along(atlas.list$lobe.names), x))
  lobes <- sapply(2:(kNumLobes-1),
          function(z) apply(t(apply(t(combos[[z]]), 1,
                function(x) atlas.list$lobe.names[x])), 1, paste, collapse=', '))
  lobes <- do.call('c', lobes)
  choices <- c('All', atlas.list$lobe.names, lobes)

  comboLobe <- gtkComboBoxNewText()
  comboLobe$show()
  for (choice in choices) comboLobe$appendText(choice)
  comboLobe$setActive(0)

  label <- gtkLabelNew('Lobe')
  hboxLobe$packStart(label, F, F, 0)
  hboxLobe$add(comboLobe)

  #---------------------------------------------------------
  # Create 2 drawing areas for the plotting
  #---------------------------------------------------------
  graphics <- vector('list', 2)
  vboxPlot <- vector('list', 2)
  groupplot <- vector('list', 2)

  # Check screen resolution so the plotting window fits
  display.size <- system('xdpyinfo | grep dimensions', intern=T)
  m <- regexpr(' [0-9]+x', display.size)
  res <- regmatches(display.size, m)
  screen.x <- as.numeric(sub('x', '', sub(' ', '', res)))

  if (screen.x <= 1440) {
    graphics.res <- 460
  } else if (screen.x > 1440 && screen.x <= 1680) {
    graphics.res <- 560
  } else {
    graphics.res <- 640
  }

  for (i in 1:2) {
    graphics[[i]] <- gtkDrawingArea()
    graphics[[i]]$setSizeRequest(graphics.res, graphics.res)

    vboxPlot[[i]] <- gtkVBox()
    vboxPlot[[i]]$packStart(graphics[[i]], expand=T, fill=T, padding=0)
    container$add(vboxPlot[[i]])

    asCairoDevice(graphics[[i]])
    par(pty='s', mar=rep(0, 4))
    groupplot[[i]] <- dev.cur()
  }

  # Slider for curvature of edges in circle plots
  slider <<- vector('list', length=2)
  orient_cb1 <- function(widget) {
                if (widget$getActive() == 1) {
                  comboHemi$setActive(1)
                  comboHemi$setSensitive(F)
                  kNumChildren <- length(vboxPlot[[1]]$getChildren())
                  if (kNumChildren > 1) {
                    for (i in 2:kNumChildren) vboxPlot[[1]][[i]]$destroy()
                  }
                } else if (widget$getActive() == 2) {
                  comboHemi$setActive(2)
                  comboHemi$setSensitive(F)
                  kNumChildren <- length(vboxPlot[[1]]$getChildren())
                  if (kNumChildren > 1) {
                    for (i in 2:kNumChildren) vboxPlot[[1]][[i]]$destroy()
                  }
                } else if (widget$getActive() == 3) {
                  comboHemi$setSensitive(T)
                  slider[[1]] <<- gtkHScale(min=-1, max=1, step=0.05)
                  vboxPlot[[1]]$packStart(slider[[1]], F, F, 0)
                  slider[[1]]$setValue(0.25)
                } else {
                  comboHemi$setSensitive(T)
                  kNumChildren <- length(vboxPlot[[1]]$getChildren())
                  if (kNumChildren > 1) {
                    for (i in 2:kNumChildren) vboxPlot[[1]][[i]]$destroy()
                  }
                }
  }
  gSignalConnect(comboOrient[[1]], 'changed', orient_cb1)
  orient_cb2 <- function(widget) {
                if (widget$getActive() == 3) {
                  slider[[2]] <<- gtkHScale(min=-1, max=1, step=0.05)
                  vboxPlot[[2]]$packStart(slider[[2]], F, F, 0)
                  slider[[2]]$setValue(0.25)
                } else {
                  kNumChildren <- length(vboxPlot[[2]]$getChildren())
                  if (kNumChildren > 1) {
                    for (i in 2:kNumChildren) vboxPlot[[2]][[i]]$destroy()
                  }
                }
  }
  gSignalConnect(comboOrient[[2]], 'changed', orient_cb2)

  #-----------------------------------------------------------------------------
  # Add buttons
  the.buttons <- gtkHButtonBoxNew()
  the.buttons$setBorderWidth(5)
  vboxMain$add(the.buttons)
  buttonOK <- gtkButtonNewFromStock('gtk-ok')
  gSignalConnect(buttonOK, 'clicked',
                 function(widget) update.adj(graphname[[1]], graphname[[2]],
                                   vertLabels, vertSize=comboVsize,
                                   edgeWidth=comboEwidth, edgeDiffs,
                                   vertColor=comboVcolor,
                                   firstplot=groupplot[[1]],
                                   secondplot=groupplot[[2]],
                                   vertSize.const, edgeWidth.const, plotFunc,
                                   comboNeighb=comboNeighb, hemi=comboHemi,
                                   lobe=comboLobe, orient1=comboOrient[[1]],
                                   orient2=comboOrient[[2]], comm=comboComm,
                                   showDiameter=showDiameter,
                                   slider1=slider[[1]], slider2=slider[[2]],
                                   vertSize.other, kNumComms=kNumComms))
  the.buttons$packStart(buttonOK, fill=F)

  container$showAll()

  list(plot1=groupplot[[1]], plot2=groupplot[[2]], graphname=graphname,
       comboVcolor=comboVcolor, showDiameter=showDiameter, edgeDiffs=edgeDiffs,
       comboHemi=comboHemi, comboLobe=comboLobe)
  }
}
