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

  gui.params <- plotFunc <- comboComm <- kNumComms <- comboNeighb <-
  graphObj <- slider <- lobe <- comboNeighbMult <- graph1 <- kNumGroups <-
  vsize.opts <- ewidth.opts <- NULL

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
    combo$setActive(0)
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
                   if (response == GtkResponseType['accept']) {
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
                   if (response == GtkResponseType['accept']) {
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
    comboNeighbMult <<- NULL
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
    lapply(gui.params$comboHemi, function(x) x$setSensitive(F))
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

        button <- gtkButtonNewFromStock('gtk-ok')
        gSignalConnect(button, 'clicked', function(widget) {
                       comboNeighbMult <<- tempComboEntry$getText()
                       tempComboWindow$destroy()
                       })
        tempComboVbox$packStart(button, fill=F)
      } else {
        comboNeighbMult <<- NULL
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
      graph1 <<- eval(parse(text=gui.params$graphname[[1]]$getText()))
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
    combos <- lapply(seq_len(kNumComms), function(x)
                     combn(seq_len(kNumComms), x))
    comms <- lapply(2:kNumComms, function(z)
                    apply(t(apply(t(combos[[z]]), 1, function(x)
                                  all.comms[x])), 1, paste, collapse=', '))
    comms <- do.call('c', comms)
    choices <- c(as.character(all.comms), comms)
    comboComm <<- add_combo(hboxComm, choices, 'Which community? (ordered by size)')

    comboNeighb <<- NULL
    comboNeighbMult <<- NULL
    gui.params$comboVcolor$setActive(1)
    gui.params$showDiameter$setSensitive(F)
    gui.params$edgeDiffs$setSensitive(F)
    lapply(gui.params$comboHemi, function(x) x$setSensitive(F))
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

  vboxMainMenu <- gtkVBoxNew(F, 2)
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
    }
    gSignalConnect(dialog, 'response', function(dialog, response) {
        if (response == GtkResponseType['ok']) {
          graphObj[[1]] <<- graphObjEntry[[1]]$getText()
          graphObj[[2]] <<- graphObjEntry[[2]]$getText()

          if (!is_igraph(eval(parse(text=graphObj[[1]])))) {
            warnDialog <- gtkMessageDialog(parent=dialog, flags='destroy-with-parent',
                                           type='error', buttons='close',
                                           'Error: Not an igraph object!')
            response <- warnDialog$run()
            if (response == GtkResponseType['close']) warnDialog$destroy()
          } else {
            if (nchar(graphObj[[2]]) > 0 & !is_igraph(eval(parse(text=graphObj[[2]])))) {
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
  vboxMain <- gtkVBoxNew(F, 2)
  container$add(vboxMain)

  #---------------------------------------------------------
  # Create frame + vert container for plot parameters
  #---------------------------------------------------------
  frame <- gtkFrameNew('Specify plotting parameters')
  vboxMain$add(frame)

  vbox <- gtkVBoxNew(F, 2)
  vbox$setBorderWidth(3)
  frame$add(vbox)

  # Create a frame for each plot area
  #---------------------------------------------------------
  frameG <- vboxG <- hboxOrient <- comboOrient <- graphname <- graphs <-
  hboxVsizeMin <- vertSize.adj <- vertSize.spin <- hboxVsizeEqn <-
  vertSizeEqn <- hboxEwidthMin <- edgeWidth.adj <- edgeWidth.spin <-
  graphics <- vboxPlot <- groupplot <- comboHemi <- vector('list', 2)
  for (i in 1:2) {
    frameG[[i]] <- gtkFrameNew(sprintf('Graph %i', i))
    vbox$add(frameG[[i]])

    vboxG[[i]] <- gtkVBoxNew(F, 2)
    vboxG[[i]]$setBorderWidth(5)
    frameG[[i]]$add(vboxG[[i]])

    hbox <- gtkHBoxNew(F, 6)
    vboxG[[i]]$packStart(hbox, expand=F, fill=F, padding=0)
    graphname[[i]] <- add_entry(hbox, label.text='Name:', char.width=20,
                                entry.text=graphObj[[i]])
    graphname[[i]]$setSensitive(F)

    hboxOrient[[i]] <- gtkHBoxNew(F, 6)
    vboxG[[i]]$add(hboxOrient[[i]])
    choices <- c('Axial', 'Sagittal (left)', 'Sagittal (right)', 'Circular')
    comboOrient[[i]] <- add_combo(hboxOrient[[i]], choices, 'Orientation')
    if (nchar(graphname[[i]]$getText()) > 0) {
      graphs[[i]] <- eval(parse(text=graphname[[i]]$getText()))
    }
    #-----------------------------------------------------------------------------
    # Both, single hemi, inter-hemi, homologous, inter/intra community/lobe?
    #-----------------------------------------------------------------------------
    hboxHemi <- gtkHBoxNew(F, 6)
    vboxG[[i]]$packStart(hboxHemi, F, F, 0)
    choices <- c('Both hemispheres', 'Left only', 'Right only',
                 'Interhemispheric only', 'Homologous only', 'Inter-community only',
                 'Intra-community only', 'Inter-lobar only', 'Intra-lobar only')
    comboHemi[[i]] <- add_combo(hboxHemi, choices, 'Hemi/edges')
  }
  kNumGroups <<- sum(vapply(graphname, function(x) nchar(x$getText()) > 0, logical(1)))
  #---------------------------------------------------------

  # Should vertex labels be displayed?
  #---------------------------------------
  hbox <- gtkHBoxNew(F, 24)
  vbox$packStart(hbox, F, F, 0)

  hboxLabels <- gtkHBoxNew(F, 6)
  hbox$packStart(hboxLabels, F, F, 0)
  vertLabels <- add_check(hboxLabels, 'Display vertex _labels?')

  # Show a legend (for lobe colors)?
  #---------------------------------------
  hboxLegend <- gtkHBoxNew(F, 6)
  hbox$packStart(hboxLegend, F, F, 0)
  showLegend <- add_check(hboxLegend, 'Show le_gend?')
  showLegend$setSensitive(F)

  # Vertex colors based on community membership?
  #---------------------------------------
  atlas <- graphs[[1]]$atlas
  hboxVcolor <- gtkHBoxNew(F, 6)
  vbox$packStart(hboxVcolor, F, F, 0)
  choices <- c('None (lightblue)', 'Communities', 'Lobes', 'Components',
               'Communities (weighted)')
  if (atlas %in% c('destrieux', 'destrieux.scgm')) choices <- c(choices, 'Class')

  comboVcolor <- add_combo(hboxVcolor, choices, 'Vertex _color')
  gSignalConnect(comboVcolor, 'changed', function(widget, ...) {
      i <- widget$getActive()
      if (i %in% c(2, 5)) {
        showLegend$setSensitive(T)
      } else {
        showLegend$setSensitive(F)
      }
  })

  #-----------------------------------------------------------------------------
  # Vertex size?
  #-----------------------------------------------------------------------------
  frameVsize <- gtkFrameNew('Vertex size')
  vbox$add(frameVsize)

  vboxVsize <- gtkVBoxNew(F, 2)
  frameVsize$add(vboxVsize)

  hboxVsize <- gtkHBoxNew(F, 4)
  vboxVsize$packStart(hboxVsize, F, F, 0)
  choices <- c('Constant', 'Degree', 'EV centrality', 'Btwn centrality',
               'Coreness', 'Clustering coeff.', 'PC',
               'E.local', 'E.nodal', 'Within-module degree z-score', 'Hub score',
               'Vulnerability', 'NN degree', 'Asymmetry', 'Eccentricity',
               'Distance', 'Distance strength', 'Lp', 'Other', 'Equation')
  if (is_weighted(graphs[[1]])) {
    choices <- c(choices, 'strength', 'knn.wt', 'E.local.wt', 'E.nodal.wt')
  }
  comboVsize <- add_combo(hboxVsize, choices, 'Attribute')
  vertSize.const <- add_entry(hboxVsize, char.width=3, entry.text='5')

  hboxVsizeOther <- gtkHBoxNew(F, 6)
  vboxVsize$packStart(hboxVsizeOther, F, F, 0)
  vertSize.other <- add_entry(hboxVsizeOther, label.text=paste0('\t', 'Other:'),
                              char.width=10)
  vertSize.other$setSensitive(F)

  # Have 2 boxes to allow for different minimums for 2 groups
  for (i in 1:2) {
    hboxVsizeMin[[i]] <- gtkHBoxNew(F, 6)
    hboxVsizeOther$packStart(hboxVsizeMin[[i]], F, F, 0)
    labelMin <- gtkLabelNew(sprintf('Min. %i:', i))
    hboxVsizeMin[[i]]$packStart(labelMin, F, F, 0)
    #vertSize.adj[[i]] <- gtkAdjustmentNew(value=0, lower=0, upper=10, step.incr=1)
    #vertSize.spin[[i]] <- gtkSpinButtonNew(vertSize.adj[[i]], climb.rate=1, digits=0)
    vertSize.spin[[i]] <- gtkSpinButtonNewWithRange(min=0, max=10, step=1)
    hboxVsizeMin[[i]]$packStart(vertSize.spin[[i]], F, F, 0)
    vertSize.spin[[i]]$setSensitive(F)
  }

  # Create 2 entries for entering a more complicated eqn
  hboxVsizeEqnMain <- gtkHBoxNew(F, 6)
  vboxVsize$packStart(hboxVsizeEqnMain, F, F, 0)
  vboxVsizeEqnMain <- gtkVBoxNew(F, 2)
  hboxVsizeEqnMain$packStart(vboxVsizeEqnMain, F, F, 0)
  for (i in 1:2) {
    hboxVsizeEqn[[i]] <- gtkHBoxNew(F, 6)
    vboxVsizeEqnMain$packStart(hboxVsizeEqn[[i]], F, F, 0)
    vertSizeEqn[[i]] <- add_entry(hboxVsizeEqn[[i]],
                                  label.text=sprintf('Eqn. %i', i),
                                  char.width=40)
    vertSizeEqn[[i]]$setSensitive(F)
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
  })

  vsize.opts <<- c('const', 'degree', 'ev.cent', 'btwn.cent',
                  'coreness', 'transitivity', 'PC', 'E.local', 'E.nodal',
                  'z.score', 'hub.score', 'vulnerability', 'knn', 'asymm',
                  'eccentricity', 'dist', 'dist.strength', 'Lp', 'other', 'eqn',
                  'strength', 'knn.wt', 'E.local.wt', 'E.nodal.wt')
  gSignalConnect(comboVsize, 'changed', function(widget, ...) {
      i <- widget$getActive()
      if (i == 0) {  # 'Constant'
        vertSize.const$setSensitive(T)
        vertSize.other$setSensitive(F)
        lapply(vertSize.spin, function(x) x$setSensitive(F))
        lapply(vertSizeEqn, function(x) x$setSensitive(F))
      } else {
        vertSize.const$setSensitive(F)
        vertSize.other$setSensitive(F)
        lapply(vertSize.spin, function(x) x$setSensitive(T))
        lapply(vertSizeEqn, function(x) x$setSensitive(F))
        if (i < 18 | i > 19) {
          if (i == 10 && !is_directed(graphs[[1]])) i <- 2
          for (j in seq_len(kNumGroups)) {
            rangeX <- range(vertex_attr(graphs[[j]], vsize.opts[i + 1]), na.rm=T)
            newMin <- rangeX[1]
            newMax <- rangeX[2]
            newStep <- ifelse(diff(rangeX) > 10 & newMin >= 0, 1,
                              ifelse(diff(rangeX) < 1, 0.01, 0.1))
            newDigits <- ifelse(diff(rangeX) > 10 & newMin >= 0, 0,
                                ifelse(diff(rangeX) < 1, 2, 1))
            gtkSpinButtonSetDigits(vertSize.spin[[j]], newDigits)
            gtkSpinButtonSetIncrements(vertSize.spin[[j]], step=newStep, page=0)
            gtkSpinButtonSetRange(vertSize.spin[[j]], min=newMin, max=newMax)
            gtkSpinButtonSetValue(vertSize.spin[[j]], newMin)
          }
        } else if (i == 18) {  # 'Other'
          vertSize.other$setSensitive(T)
          lapply(vertSizeEqn, function(x) x$setSensitive(F))
          lapply(vertSize.spin, function(x) x$setSensitive(T))
          for (j in seq_len(kNumGroups)) {
            gtkSpinButtonSetDigits(vertSize.spin[[j]], 2)
            gtkSpinButtonSetIncrements(vertSize.spin[[j]], step=0.01, page=0)
            gtkSpinButtonSetRange(vertSize.spin[[j]], min=-100, max=100)
            gtkSpinButtonSetValue(vertSize.spin[[j]], 0)
          }
        } else if (i == 19) {  # equation
          lapply(vertSizeEqn, function(x) x$setSensitive(T))
          lapply(vertSize.spin, function(x) x$setSensitive(F))
        }
      }
  })
  # Commented out b/c clicks were causing spin value to change multiple times
  #gSignalConnect(vertSize.spin[[1]], 'value-changed', function(x) updateFun(1))
  #gSignalConnect(vertSize.spin[[2]], 'value-changed', function(x) updateFun(2))
  #-----------------------------------------------------------------------------
  # Edge width?
  #-----------------------------------------------------------------------------
  frameEwidth <- gtkFrameNew('Edge width')
  vbox$add(frameEwidth)

  vboxEwidth <- gtkVBoxNew(F, 2)
  frameEwidth$add(vboxEwidth)

  hboxEwidth <- gtkHBoxNew(F, 6)
  vboxEwidth$packStart(hboxEwidth, F, F, 0)
  choices <- c('Constant', 'Edge betweenness', 'Distance', 'Other')
  if ('weight' %in% edge_attr_names(graphs[[1]])) choices <- c(choices, 'Weight')
  comboEwidth <- add_combo(hboxEwidth, choices, 'Attribute')
  edgeWidth.const <- add_entry(hboxEwidth, char.width=3, entry.text='1')

  hboxEwidthOther <- gtkHBoxNew(F, 6)
  vboxEwidth$packStart(hboxEwidthOther, F, F, 0)
  edgeWidth.other <- add_entry(hboxEwidthOther, label.text=paste0('\t', 'Other:'),
                              char.width=10)
  edgeWidth.other$setSensitive(F)

  # Have 2 entries to allow for min's of 2 groups
  for (i in 1:2) {
    hboxEwidthMin[[i]] <- gtkHBoxNew(F, 6)
    hboxEwidthOther$packStart(hboxEwidthMin[[i]], T, F, 0)
    labelMin <- gtkLabelNew(sprintf('Min. %i:', i))
    hboxEwidthMin[[i]]$packStart(labelMin, F, F, 0)
    edgeWidth.adj[[i]] <- gtkAdjustmentNew(value=0, lower=0, upper=100, step.incr=1)
    edgeWidth.spin[[i]] <- gtkSpinButtonNew(edgeWidth.adj[[i]], climb.rate=1, digits=0)
    hboxEwidthMin[[i]]$packStart(edgeWidth.spin[[i]], F, F, 0)
    edgeWidth.spin[[i]]$setSensitive(F)
  }

  ewidth.opts <<- c('const', 'btwn', 'dist', 'other', 'weight')
  gSignalConnect(comboEwidth, 'changed', function(widget, ...) {
    i <- widget$getActive()
    if (i == 0) {  # 'Constant'
      edgeWidth.const$setSensitive(T)
      edgeWidth.other$setSensitive(F)
      lapply(edgeWidth.spin, function(x) x$setSensitive(F))
    } else if (i == 3) {  # 'Other'
      edgeWidth.const$setSensitive(F)
      edgeWidth.other$setSensitive(T)
      lapply(edgeWidth.spin, function(x) x$setSensitive(T))
      for (j in seq_len(kNumGroups)) {
        gtkSpinButtonSetDigits(edgeWidth.spin[[j]], 2)
        gtkSpinButtonSetIncrements(edgeWidth.spin[[j]], step=0.01, page=0)
        gtkSpinButtonSetRange(edgeWidth.spin[[j]], min=-100, max=100)
        gtkSpinButtonSetValue(edgeWidth.spin[[j]], 0)
      }
    } else {
      edgeWidth.const$setSensitive(F)
      edgeWidth.other$setSensitive(F)
      lapply(edgeWidth.spin, function(x) x$setSensitive(T))
      for (j in seq_len(kNumGroups)) {
        rangeX <- range(edge_attr(graphs[[j]], ewidth.opts[i + 1]), na.rm=T)
        newMin <- rangeX[1]
        newMax <- rangeX[2]
        newStep <- ifelse(diff(rangeX) > 1 & newMin >= 0, 1,
                          ifelse(diff(rangeX) < 1, 0.01, 0.1))
        newDigits <- ifelse(diff(rangeX) > 1 & newMin >= 0, 0, 2)
        gtkSpinButtonSetDigits(edgeWidth.spin[[j]], newDigits)
        gtkSpinButtonSetIncrements(edgeWidth.spin[[j]], step=newStep, page=0)
        gtkSpinButtonSetRange(edgeWidth.spin[[j]], min=newMin, max=newMax)
        gtkSpinButtonSetValue(edgeWidth.spin[[j]], newMin)
      }
    }
  })
  #-----------------------------------------------------------------------------

  # Highlight the diameter of each graph?
  hbox <- gtkHBoxNew(F, 6)
  vbox$packStart(hbox, F, F, 0)
  showDiameter <- add_check(hbox, 'Show _diameter?')

  # Show edge set differences?
  edgeDiffs <- add_check(hbox, 'Show _edge differences?')
  if (kNumGroups <= 1) edgeDiffs$setSensitive(F)

  # Major lobe number, if applicable
  atlas.dt <- eval(parse(text=atlas))
  hboxLobe <- gtkHBoxNew(F, 6)
  vbox$packStart(hboxLobe, F, F, 0)

  kNumLobes <- nlevels(atlas.dt[, lobe])
  combos <- lapply(seq_len(kNumLobes - 1),
                   function(x) combn(seq_along(levels(atlas.dt[, lobe])), x))
  lobes <- lapply(2:(kNumLobes-1),
          function(z) apply(t(apply(t(combos[[z]]), 1,
                function(x) levels(atlas.dt[, lobe])[x])), 1, paste, collapse=', '))
  lobes <- do.call('c', lobes)
  choices <- c('All', levels(atlas.dt[, lobe]), lobes)

  comboLobe <- add_combo(hboxLobe, choices, 'Lobe')
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
  if (.Platform$OS.type == 'windows') {
    scr_width <- system("wmic desktopmonitor get screenwidth", intern=TRUE)
    scr_height <- system("wmic desktopmonitor get screenheight", intern=TRUE)
    res <- as.numeric(c(scr_width[-c(1, length(scr_width))],
                        scr_height[-c(1, length(scr_height))]))
  } else {
    display.size <- system("xrandr | grep '*' | awk '{print $1}'", intern=T)
    res <- as.numeric(strsplit(display.size, 'x')[[1]])
  }

  ydim <- ifelse(res[2] < 800, 0.8 * res[2], 700)

  for (i in 1:2) {
    graphics[[i]] <- gtkDrawingArea()
    graphics[[i]]$setSizeRequest(res[1] / 3 - 10, ydim)#0.8 * res[2])#graphics.res, graphics.res)

    vboxPlot[[i]] <- gtkVBox()
    vboxPlot[[i]]$packStart(graphics[[i]], expand=T, fill=T, padding=0)
    container$add(vboxPlot[[i]])

    asCairoDevice(graphics[[i]])
    par(pty='s', mar=rep(0, 4))
    groupplot[[i]] <- dev.cur()
  }

  #-----------------------------------------------------------------------------
  # Add buttons
  the.buttons <- gtkHButtonBoxNew()
  the.buttons$setBorderWidth(5)
  vboxMain$add(the.buttons)
  btnOK <- gtkButtonNewFromStock('gtk-ok')
  buttonRename <- gtkButtonNewWithMnemonic('Pick _new graphs')
  the.buttons$packStart(btnOK, expand=T, fill=F)
  the.buttons$packStart(buttonRename, expand=T, fill=F)
  container$showAll()

  # Callback functions
  updateFun <- function(j) {
    other.ind <- setdiff(seq_len(kNumGroups), j)
    update_brainGraph_gui(plotDev=groupplot[[j]], graph1=graphs[[j]],
                          graph2=graphs[[other.ind]], plotFunc=plotFunc,
                          vertSize=comboVsize, edgeWidth=comboEwidth,
                          vertColor=comboVcolor, hemi=comboHemi[[j]], lobe=comboLobe,
                          orient=comboOrient[[j]], vertSize.min=vertSize.spin[[j]],
                          edgeWidth.min=edgeWidth.spin[[j]],
                          vertSize.const=vertSize.const, edgeWidth.const=edgeWidth.const,
                          vertLabels=vertLabels, showLegend=showLegend,
                          comm=comboComm, kNumComms=kNumComms,
                          neighb=comboNeighb, neighbMult=comboNeighbMult,
                          slider=slider[[j]], vertSize.other, edgeWidth.other,
                          vertSize.eqn=vertSizeEqn[[j]], showDiameter=showDiameter, edgeDiffs)
  }
  gSignalConnect(btnOK, 'clicked', function(widget) {
                 updateFun(1)
                 if (kNumGroups > 1) updateFun(2)
  })
  gSignalConnect(buttonRename, 'clicked', function(widget) set.names())

  # Horizontal slider for edge curvature in circle plots
  slider <<- vector('list', length=2)
  orient_cb <- function(widget, ind) {
    if (widget$getActive() == 3) {
      comboHemi[[ind]]$setSensitive(T)
      slider[[ind]] <<- gtkHScaleNewWithRange(min=-1, max=1, step=0.05)
      vboxPlot[[ind]]$packStart(slider[[ind]], F, F, 0)
      slider[[ind]]$setValue(0.25)
      gSignalConnect(slider[[ind]], 'value-changed',
                     function(widget) updateFun(ind))
    } else {
      kNumChildren <- length(vboxPlot[[ind]]$getChildren())
      if (kNumChildren > 1) {
        for (i in 2:kNumChildren) vboxPlot[[ind]][[i]]$destroy()
      }
      if (widget$getActive() == 1) {
        comboHemi[[ind]]$setActive(1)
        comboHemi[[ind]]$setSensitive(F)
      } else if (widget$getActive() == 2) {
        comboHemi[[ind]]$setActive(2)
        comboHemi[[ind]]$setSensitive(F)
      } else {
        comboHemi[[ind]]$setSensitive(T)
      }
    }
    updateFun(ind)
  }
  gSignalConnect(comboOrient[[1]], 'changed', orient_cb, 1)
  gSignalConnect(comboOrient[[2]], 'changed', orient_cb, 2)

  list(plot1=groupplot[[1]], plot2=groupplot[[2]], graphname=graphname,
       comboVcolor=comboVcolor, showDiameter=showDiameter, edgeDiffs=edgeDiffs,
       comboHemi=comboHemi, comboLobe=comboLobe)
  }
}
