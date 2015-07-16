#' GUI for plotting graphs overlaid on an MNI152 image or in a circle.
#'
#' This function creates a GUI for plotting graphs over an image from the MNI152
#' template. It gives the user control over several plotting parameters. Also
#' possible is a circular plot (in addition to the axial and sagittal views). It
#' is necessary for the graphs to have an \emph{atlas} attribute, and several
#' vertex- and edge-level attributes (set by
#' \code{\link{set.brainGraph.attributes}}).
#'
#' @export

plot_brainGraph_gui <- function() {
  window <- gtkWindow('toplevel')
  window['title'] <- 'brainGraph'
  window['icon'] <- gdkPixbuf(filename=system.file('extdata',
                                                   'brainGraph_icon.png',
                                                   package='brainGraph'))

  gui.params <- plotFunc <- comboComm <- kNumComms <- comboNeighb <- NULL
  graphObj <- slider <- lobe <- comboNeighbMult <- NULL

  # Function to add an entry
  add_entry <- function(container, label.text=NULL, char.width, entry.text=NULL) {
    if (!is.null(label.text)) {
      label <- gtkLabelNewWithMnemonic(label.text)
      container$packStart(label, F, F, 0)
    }
    entry <- gtkEntryNew()
    entry$setWidthChars(char.width)
    if (!is.null(entry.text)) entry$setText(entry.text)
    container$packStart(entry, F, F, 0)
    entry$setSensitive(F)
    return(entry)
  }

  # Function to add a check button & label
  add_check <- function(container, label.text) {
    label <- gtkLabelNewWithMnemonic(label.text)
    container$packStart(label, F, F, 0)
    chk.button <- gtkCheckButton()
    chk.button$active <- FALSE
    container$packStart(chk.button, F, F, 0)
    label$setMnemonicWidget(chk.button)
    return(chk.button)
  }

  # Function to add a combobox
  add_combo <- function(container, choices, label.text) {
    combo <- gtkComboBoxNewText()
    combo$show()
    for (choice in choices) combo$appendText(choice)
    combo$setActive(1)
    label <- gtkLabelNewWithMnemonic(label.text)
    container$packStart(label, F, F, 0)
    container$add(combo)
    return(combo)
  }
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
    list('Save1', 'gtk-save', 'Save plot _1', '<control>1', 'Save plot 1?', save_cb1),
    list('Save2', 'gtk-save', 'Save plot _2', '<control>2', 'Save plot 2?', save_cb2),
    list('Quit', 'gtk-quit', '_Quit', '<control>Q', 'Quit the application',
         quit_cb)
  )
  action_group <- gtkActionGroup('brainGraphActions')
  action_group$addActions(file.actions, window)

  #---------------------------------------------------------
  # Callback functions for the "Plot" menu
  #---------------------------------------------------------
  plot_entire_cb <- function(widget, window) {
    plotFunc <<- plot_brainGraph
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
    plotFunc <<- 'plot_neighborhood'
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

    graph1 <- eval(parse(text=gui.params$graphname[[1]]$getText()))
    choicesNeighb <- data.frame(seq_len(vcount(graph1)), V(graph1)$name)
    choices <- vector('character', length=(nrow(choicesNeighb) + 1))
    choices[1] <- '[Multiple]'
    for (r in seq_len(nrow(choicesNeighb))) {
      choices[r+1] <- paste(sprintf('%02i', choicesNeighb[r, 1]), ': ', '\t',
                          as.character(choicesNeighb[r, 2]))
    }
    comboNeighb <<- add_combo(hboxNeighb, choices, 'Which vertex?')
    comboNeighb$setActive(1)

    comboComm <<- NULL
    kNumComms <<- NULL
    gui.params$comboVcolor$setActive(0)
    gui.params$showDiameter$setSensitive(F)
    gui.params$edgeDiffs$setSensitive(F)
    gui.params$comboHemi$setSensitive(F)
    gui.params$comboLobe$setSensitive(F)

    gSignalConnect(comboNeighb, 'changed', function(comboNeighb, ...) {
      if (comboNeighb$getActive() == 0) {
        tempComboWindow <- gtkWindow()
        tempComboVbox <- gtkVBoxNew(F, 8)
        tempComboWindow$add(tempComboVbox)
        tempComboHbox <- gtkHBoxNew(F, 8)
        tempComboVbox$packStart(tempComboHbox, F, F, 0)
        helpText <- sprintf('%s\n\n%s\n%s',
                            'Choose vertices (comma-separated list)',
                            'Example A: 24, 58',
                            'Example B: lPCUN, rPCUN')
        tempLabel <- gtkLabelNew(helpText)
        tempComboVbox$packStart(tempLabel, F, F, 0)
        tempComboEntry <- add_entry(tempComboVbox, char.width=25)
        tempComboEntry$setSensitive(T)

        button <- gtkButtonNewFromStock('gtk-ok')
        gSignalConnect(button, 'clicked', function(widget) {
                       comboNeighbMult <<- tempComboEntry$getText()
                       tempComboWindow$destroy()
                       })
        tempComboVbox$packStart(button, fill=F)
      }
    })
  }

  #-------------------------------------
  # Plot individual communities
  #-------------------------------------
  plot_comm_cb <- function(widget, window) {
    plotFunc <<- 'plot_community'
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
    comboComm <<- add_combo(hboxComm, choices, 'Which community? (ordered by size)')
    comboComm$setActive(0)

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
    dialog <- gtkDialogNewWithButtons(title='Choose graph objects',
                                      parent=window, flags='destroy-with-parent',
                                      'gtk-ok', GtkResponseType['ok'],
                                      'gtk-cancel', GtkResponseType['cancel'],
                                      show=FALSE)
    tempVbox <- dialog$getContentArea()
    label <- graphObjEntry <- tempHbox <- vector('list', length=2)
    for (i in 1:2) {
      tempHbox[[i]] <- gtkHBoxNew(F, 8)
      tempVbox$packStart(tempHbox[[i]], F, F, 0)
      graphObjEntry[[i]] <- add_entry(tempHbox[[i]],
                                      label.text=paste0('Graph ', i, ' name:'),
                                      char.width=20)
      graphObjEntry[[i]]$setSensitive(T)
    }
    gSignalConnect(dialog, 'response', function(dialog, response, user.data) {
        if (response == GtkResponseType['ok']) {
          graphObj[[1]] <<- graphObjEntry[[1]]$getText()
          graphObj[[2]] <<- graphObjEntry[[2]]$getText()

          if (!is.igraph(eval(parse(text=graphObj[[1]])))) {
            warnDialog <- gtkMessageDialog(parent=dialog, flags='destroy-with-parent',
                                           type='error', buttons='close',
                                           'Error: Not an igraph object!')
            response <- warnDialog$run()
            if (response == GtkResponseType['close']) warnDialog$destroy()
          } else {
            if (nchar(graphObj[[2]]) > 0 & !is.igraph(eval(parse(text=graphObj[[2]])))) {
              warnDialog <- gtkMessageDialog(parent=dialog, flags='destroy-with-parent',
                                             type='error', buttons='close',
                                             'Error: Not an igraph object!')
              response <- warnDialog$run()
              if (response == GtkResponseType['close']) warnDialog$destroy()
            } else {
              dialog$Destroy()
            }
          }
        } else {
          window$destroy()
        }
    })
    dialog$showAll()
    dialog$setModal(TRUE)
  }
  set.names()

  #-------------------------------------------------------------------------------
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

  vbox <- gtkVBoxNew(F, 6)
  vbox$setBorderWidth(5)
  frame$add(vbox)

  # Create a frame for each plot area
  #---------------------------------------------------------
  frameG <- vboxG <- hboxOrient <- comboOrient <- graphname <- vector('list', 2)
  for (i in 1:2) {
    frameG[[i]] <- gtkFrameNew(sprintf('Graph %i', i))
    vbox$add(frameG[[i]])

    vboxG[[i]] <- gtkVBoxNew(F, 8)
    vboxG[[i]]$setBorderWidth(10)
    frameG[[i]]$add(vboxG[[i]])

    hbox <- gtkHBoxNew(F, 8)
    vboxG[[i]]$packStart(hbox, expand=F, fill=F, padding=0)

    graphname[[i]] <- add_entry(hbox, label.text='Name:', char.width=20,
                                entry.text=graphObj[[i]])

    hboxOrient[[i]] <- gtkHBoxNew(F, 8)
    vboxG[[i]]$add(hboxOrient[[i]])
    choices <- c('Axial', 'Sagittal (left)', 'Sagittal (right)', 'Circular')
    comboOrient[[i]] <- add_combo(hboxOrient[[i]], choices, 'Orientation')
    comboOrient[[i]]$setActive(0)
  }
  #---------------------------------------------------------

  # Should vertex labels be displayed?
  #---------------------------------------
  hbox <- gtkHBoxNew(F, 8)
  vbox$packStart(hbox, F, F, 0)
  vertLabels <- add_check(hbox, 'Display vertex _labels?')

  # Vertex colors based on community membership?
  #---------------------------------------
  atlas <- eval(parse(text=graphname[[1]]$getText()))$atlas
  hboxVcolor <- gtkHBoxNew(F, 8)
  vbox$packStart(hboxVcolor, F, F, 0)
  choices <- c('None (lightblue)', 'Communities', 'Lobes', 'Components')
  if (atlas == 'destrieux') {
    choices <- c(choices, 'Class')
  }

  comboVcolor <- add_combo(hboxVcolor, choices, 'Vertex _color')

  #-----------------------------------------------------------------------------
  # Vertex size?
  #-----------------------------------------------------------------------------
  frameVsize <- gtkFrameNew('Vertex size')
  vbox$add(frameVsize)

  vboxVsize <- gtkVBoxNew(F, 8)
  frameVsize$add(vboxVsize)

  hboxVsize <- gtkHBoxNew(F, 8)
  vboxVsize$packStart(hboxVsize, F, F, 0)
  choices <- c('Constant', 'Degree', 'EV centrality', 'Btwn centrality',
               'Subgraph centrality', 'Coreness', 'Clustering coeff.', 'PC',
               'E.local', 'E.nodal', 'Within-module degree z-score', 'Hub score',
               'Vulnerability', 'NN degree', 'Other', 'Equation')
  comboVsize <- add_combo(hboxVsize, choices, 'Attribute')
  vertSize.const <- add_entry(hboxVsize, char.width=3, entry.text='5')

  hboxVsizeOther <- gtkHBoxNew(F, 8)
  vboxVsize$packStart(hboxVsizeOther, F, F, 0)
  vertSize.other <- add_entry(hboxVsizeOther, label.text=paste0('\t', 'Other:'),
                              char.width=10)
  # Have 2 boxes to allow for different minimums for 2 groups
  hboxVsizeMin <- vertSize.min <- vector('list', 2)
  for (i in 1:2) {
    hboxVsizeMin[[i]] <- gtkHBoxNew(F, 8)
    hboxVsizeOther$packStart(hboxVsizeMin[[i]], F, F, 0)
    vertSize.min[[i]] <- add_entry(hboxVsizeMin[[i]],
                                   label.text=sprintf('Min. %i:', i),
                                   char.width=6, entry.text='0')
    vertSize.min[[i]]$setSensitive(T)
  }

  # Create 2 entries for entering a more complicated eqn
  hboxVsizeEqnMain <- gtkHBoxNew(F, 8)
  vboxVsize$packStart(hboxVsizeEqnMain, F, F, 0)
  vboxVsizeEqnMain <- gtkVBoxNew(F, 8)
  hboxVsizeEqnMain$packStart(vboxVsizeEqnMain, F, F, 0)
  hboxVsizeEqn <- vertSizeEqn <- vector('list', 2)
  for (i in 1:2) {
    hboxVsizeEqn[[i]] <- gtkHBoxNew(F, 8)
    vboxVsizeEqnMain$packStart(hboxVsizeEqn[[i]], F, F, 0)
    vertSizeEqn[[i]] <- add_entry(hboxVsizeEqn[[i]],
                                  label.text=sprintf('Eqn. %i', i),
                                  char.width=40)
  }

  vertSizeHelp <- gtkButtonNewWithLabel('?')
  hboxVsizeEqnMain$packStart(vertSizeHelp, F, F, 0)
  gSignalConnect(vertSizeHelp, 'clicked', function(widget, ...) {
    helpWin <- gtkWindow()
    helpVbox <- gtkVBoxNew(F, 8)
    helpWin$add(helpVbox)
    helpHbox <- gtkHBoxNew(F, 8)
    helpVbox$packStart(helpHbox, F, F, 0)
    helpText <- sprintf('%s\n\n\n%s',
                        'Instructions: enter simple logical expressions
                        separated by a single ampersand (&)
                        The vertex attributes must already exist.',
                        'Example: degree > 23 & btwn.cent > 50')
    helpLabel <- gtkLabelNew(helpText)
    helpVbox$packStart(helpLabel, F, F, 0)
    helpButton <- gtkButtonNewFromStock('gtk-ok')
    gSignalConnect(helpButton, 'clicked', function(widget) helpWin$destroy())
    helpVbox$packStart(helpButton, fill=F)
                                  }
  )

  gSignalConnect(comboVsize, 'changed', function(widget, ...) {
      if (widget$getActive() == 0) {  # 'Constant'
        vertSize.const$setSensitive(T)
        vertSize.other$setSensitive(F)
        vertSize.min[[1]]$setSensitive(F)
        vertSize.min[[2]]$setSensitive(F)
        vertSizeEqn[[1]]$setSensitive(F)
        vertSizeEqn[[2]]$setSensitive(F)
      } else if (widget$getActive() == 14) {  # 'Other'
        vertSize.const$setSensitive(F)
        vertSize.other$setSensitive(T)
        vertSizeEqn[[1]]$setSensitive(F)
        vertSizeEqn[[2]]$setSensitive(F)
        vertSize.min[[1]]$setSensitive(T)
        vertSize.min[[2]]$setSensitive(T)
      } else if (widget$getActive() == 15) {  # equation
        vertSize.const$setSensitive(F)
        vertSizeEqn[[1]]$setSensitive(T)
        vertSizeEqn[[2]]$setSensitive(T)
        vertSize.min[[1]]$setSensitive(F)
        vertSize.min[[2]]$setSensitive(F)
      } else {
        vertSize.const$setSensitive(F)
        vertSize.other$setSensitive(F)
        vertSizeEqn[[1]]$setSensitive(F)
        vertSizeEqn[[2]]$setSensitive(F)
        vertSize.min[[1]]$setSensitive(T)
        vertSize.min[[2]]$setSensitive(T)
      }
    })
  #-----------------------------------------------------------------------------
  # Edge width?
  hboxEwidth <- gtkHBoxNew(F, 8)
  vbox$packStart(hboxEwidth, F, F, 0)
  choices <- c('Constant', 'Edge betweenness', 'Distance')
  comboEwidth <- add_combo(hboxEwidth, choices, 'Edge width')

  edgeWidth.const <- add_entry(hboxEwidth, char.width=3, entry.text='1')

  # Have 2 entries to allow for min's of 2 groups
  hboxEwidthOther <- gtkHBoxNew(F, 8)
  vbox$packStart(hboxEwidthOther, F, F, 0)
  hboxEwidthMin <- edgeWidth.min <- vector('list', 2)
  for (i in 1:2) {
    hboxEwidthMin[[i]] <- gtkHBoxNew(F, 8)
    hboxEwidthOther$packStart(hboxEwidthMin[[i]], F, F, 0)
    edgeWidth.min[[i]] <- add_entry(hboxEwidthMin[[i]],
                                    label.text=sprintf('\t\tMin. %i:', i),
                                    char.width=6, entry.text='0')
    edgeWidth.min[[i]]$setSensitive(T)
  }

  gSignalConnect(comboEwidth, 'changed', function(widget, ...) {
      if (widget$getActive() == 0) {
        edgeWidth.const$setSensitive(T)
        edgeWidth.min[[1]]$setSensitive(F)
        edgeWidth.min[[2]]$setSensitive(F)
      } else {
        edgeWidth.const$setSensitive(F)
        edgeWidth.min[[1]]$setSensitive(T)
        edgeWidth.min[[2]]$setSensitive(T)
      }
    })
  #-----------------------------------------------------------------------------

  # Highlight the diameter of each graph?
  hbox <- gtkHBoxNew(F, 8)
  vbox$packStart(hbox, F, F, 0)
  showDiameter <- add_check(hbox, 'Show _diameter?')

  # Show edge set differences?
  edgeDiffs <- add_check(hbox, 'Show _edge differences?')

  #-----------------------------------------------------------------------------
  # Both, single hemisphere, inter-hemispheric, or homologous only?
  #---------------------------------------------------------
  hboxHemi <- gtkHBoxNew(F, 8)
  vbox$packStart(hboxHemi, F, F, 0)
  choices <- c('Both', 'Left only', 'Right only', 'Interhemispheric only',
               'Homologous only')
  comboHemi <- add_combo(hboxHemi, choices, 'Hemisphere')
  comboHemi$setActive(0)
  gSignalConnect(comboHemi, 'changed', function(widget, ...) {
      if (widget$getActive() == 0) {
        showDiameter$setSensitive(T)
        edgeDiffs$setSensitive(T)
      } else {
        showDiameter$setSensitive(F)
        edgeDiffs$setSensitive(F)
      }
  })

  # Major lobe number, if applicable
  atlas.dt <- eval(parse(text=data(list=atlas)))
  hboxLobe <- gtkHBoxNew(F, 8)
  vbox$packStart(hboxLobe, F, F, 0)

  kNumLobes <- nlevels(atlas.dt[, lobe])
  combos <- sapply(seq_len(kNumLobes - 1),
                   function(x) combn(seq_along(levels(atlas.dt[, lobe])), x))
  lobes <- sapply(2:(kNumLobes-1),
          function(z) apply(t(apply(t(combos[[z]]), 1,
                function(x) levels(atlas.dt[, lobe])[x])), 1, paste, collapse=', '))
  lobes <- do.call('c', lobes)
  choices <- c('All', levels(atlas.dt[, lobe]), lobes)

  comboLobe <- add_combo(hboxLobe, choices, 'Lobe')
  comboLobe$setActive(0)
  gSignalConnect(comboLobe, 'changed', function(widget, ...) {
      if (widget$getActive() == 0) {
        showDiameter$setSensitive(T)
      } else {
        showDiameter$setSensitive(F)
      }
  })

  #---------------------------------------------------------
  # Create 2 drawing areas for the plotting
  #---------------------------------------------------------
  graphics <- vboxPlot <- groupplot <- vector('list', 2)

  # Check screen resolution so the plotting window fits
  OS <- .Platform$OS.type
  if (OS == 'windows') {
    screen.x <- 1600
  } else {
  display.size <- system('xdpyinfo | grep dimensions', intern=T)
  m <- regexpr(' [0-9]+x', display.size)
  res <- regmatches(display.size, m)
  screen.x <- as.numeric(sub('x', '', sub(' ', '', res)))
  }

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
  orient_cb <- function(widget, ind) {
    if (widget$getActive() == 3) {
      comboHemi$setSensitive(T)
      slider[[ind]] <<- gtkHScale(min=-1, max=1, step=0.05)
      vboxPlot[[ind]]$packStart(slider[[ind]], F, F, 0)
      slider[[ind]]$setValue(0.25)
    } else {
      kNumChildren <- length(vboxPlot[[ind]]$getChildren())
      if (kNumChildren > 1) {
        for (i in 2:kNumChildren) vboxPlot[[ind]][[i]]$destroy()
      }
      if (widget$getActive() == 1) {
        comboHemi$setActive(1)
        comboHemi$setSensitive(F)
      } else if (widget$getActive() == 2) {
        comboHemi$setActive(2)
        comboHemi$setSensitive(F)
      } else {
        comboHemi$setSensitive(T)
      }
    }
  }
  gSignalConnect(comboOrient[[1]], 'changed', orient_cb, 1)
  gSignalConnect(comboOrient[[2]], 'changed', orient_cb, 2)

  #-----------------------------------------------------------------------------
  # Add buttons
  the.buttons <- gtkHButtonBoxNew()
  the.buttons$setBorderWidth(5)
  vboxMain$add(the.buttons)
  buttonOK <- gtkButtonNewFromStock('gtk-ok')
  gSignalConnect(buttonOK, 'clicked',
                 function(widget)
                     update_brainGraph_gui(graphname[[1]], graphname[[2]],
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
                                   vertSize.other, vertSize.min1=vertSize.min[[1]],
                                   vertSize.min2=vertSize.min[[2]],
                                   edgeWidth.min1=edgeWidth.min[[1]],
                                   edgeWidth.min2=edgeWidth.min[[2]],
                                   kNumComms=kNumComms, comboNeighbMult,
                                   vertSize.eqn1=vertSizeEqn[[1]],
                                   vertSize.eqn2=vertSizeEqn[[2]]))
  buttonRename <- gtkButtonNewWithMnemonic('Pick _new graphs')
  gSignalConnect(buttonRename, 'clicked',
                  function(widget) set.names())
  the.buttons$packStart(buttonOK, expand=T, fill=F)
  the.buttons$packStart(buttonRename, expand=T, fill=F)

  container$showAll()

  list(plot1=groupplot[[1]], plot2=groupplot[[2]], graphname=graphname,
       comboVcolor=comboVcolor, showDiameter=showDiameter, edgeDiffs=edgeDiffs,
       comboHemi=comboHemi, comboLobe=comboLobe)
  }
}
