#' GUI for plotting graphs overlaid on an MNI152 image or in a circle.
#'
#' This function creates a GUI for plotting graphs over an image from the MNI152
#' template. It gives the user control over several plotting parameters. Also
#' possible is a circular plot (in addition to the axial and sagittal views). It
#' is necessary for the graphs to have an \emph{atlas} attribute, and several
#' vertex- and edge-level attributes (set by
#' \code{\link{set_brainGraph_attr}}).
#'
#' @export
#' @family Plotting functions

plot_brainGraph_gui <- function() {
  if (!(requireNamespace("RGtk2", quietly = TRUE) && (requireNamespace("cairoDevice", quietly = TRUE)))) {
    warning(paste(c("You need to install RGtk2 and cairoDevice to use the brainGraph GUI.")))
    return(NULL)
  } else {
    require(cairoDevice)
    require(RGtk2)
  }
  
  window <- gtkWindow('toplevel')
  window['title'] <- 'brainGraph'
  window['icon'] <- gdkPixbuf(filename=system.file('extdata',
                                                   'brainGraph_icon.png',
                                                   package='brainGraph'))

  gui.params <- hboxMain <- plotFunc <- neighbs <- view <- comms <- myComm <-
    myNeighb <- myNeighbInd <- graphObj <- slider <- graph1 <- kNumGroups <-
    vsize.opts <- vsize.measure <- ewidth.opts <- ewidth.measure <- showDiameter <-
    edgeDiffs <- myLobe <- lobe <- Lobe <- NULL

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
  save_cb <- function(widget, window, plot.dev) {
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
                   dialog$destroy()})
  }
  save_cb1 <- function(widget, window) save_cb(widget, window, plot.dev=gui.params$plot1)
  save_cb2 <- function(widget, window) save_cb(widget, window, plot.dev=gui.params$plot2)
  quit_cb <- function(widget, window) window$destroy()

  # Create a "GtkAction" class
  file.actions <- list(
    list('FileMenu', NULL, '_File'),
    list('Save1', 'gtk-save', 'Save plot _1', '<control>1', 'Save plot 1?', save_cb1),
    list('Save2', 'gtk-save', 'Save plot _2', '<control>2', 'Save plot 2?', save_cb2),
    list('Quit', 'gtk-quit', '_Quit', '<control>Q', 'Quit the application', quit_cb)
  )
  action_group <- gtkActionGroup('brainGraphActions')
  action_group$addActions(file.actions, window)

  #---------------------------------------------------------
  # Callback functions for the "Plot" menu
  #---------------------------------------------------------
  kill_others <- function(window) {
    # Get number of children of vboxMainMenu and kill all but 'menubar'
    kNumChildren <- length(window[[1]]$getChildren())
    if (kNumChildren > 1) {
      for (i in 2:kNumChildren) window[[1]][[i]]$destroy()
    }
    hboxMain <<- gtkHBoxNew(F, 8)
    window[[1]]$add(hboxMain)
    gui.params <<- build.gui(hboxMain)
  }

  # Plot whole graph
  #-------------------------------------
  plot_entire_cb <- function(widget, window) {
    plotFunc <<- plot.brainGraph
    kill_others(window)
  }

  # Plot neighborhoods
  #-------------------------------------
  plot_neighb_cb <- function(widget, window) {
    plotFunc <<- 'plot_neighborhood'
    kill_others(window)
    gui.params$comboVcolor$appendText('Neighborhoods')
    gui.params$comboVcolor$setActive(0)
    lapply(gui.params$comboHemi, function(x) x$setSensitive(F))

    graph1 <- eval(parse(text=gui.params$graphname[[1]]$getText()))
    neighbs <<- data.frame(index=seq_along(V(graph1)), name=V(graph1)$name)
    model <- rGtkDataFrame(neighbs)
    view <<- gtkTreeView(model)
    view$getSelection()$setMode('multiple')
    column1 <- gtkTreeViewColumn('Index', gtkCellRendererText(), text=0)
    view$appendColumn(column1)
    column2 <- gtkTreeViewColumn('Name', gtkCellRendererText(), text=1)
    view$appendColumn(column2)
    scrolled_window <- gtkScrolledWindow()
    scrolled_window$setSizeRequest(-1, 124)
    scrolled_window$add(view)
    gui.params$vbox$add(scrolled_window)
  }

  # Plot individual communities
  #-------------------------------------
  plot_comm_cb <- function(widget, window) {
    plotFunc <<- 'plot_community'
    kill_others(window)
    gui.params$comboVcolor$setActive(1)
    lapply(gui.params$comboHemi, function(x) x$setSensitive(F))

    graph1 <- eval(parse(text=gui.params$graphname[[1]]$getText()))
    all.comms <- which(table(V(graph1)$comm) > 2)
    if (!gui.params$graphname[[2]]$getText() == '') {
      graph2 <- eval(parse(text=gui.params$graphname[[2]]$getText()))
      c2 <- which(table(V(graph2)$comm) > 2)
      all.comms <- union(all.comms, c2)
    }
    comms <<- data.frame(Community=all.comms)
    model <- rGtkDataFrame(comms)
    view <<- gtkTreeView(model)
    view$getSelection()$setMode('multiple')
    column <- gtkTreeViewColumn('Community', gtkCellRendererText(), text=0)
    view$appendColumn(column)
    scrolled_window <- gtkScrolledWindow()
    scrolled_window$setSizeRequest(-1, 124)
    scrolled_window$add(view)
    gui.params$vbox$add(scrolled_window)
  }

  # Action group for the "Plot" menu
  plot.actions <- list(
    list('PlotMenu', NULL, '_Plot'),
    list('PlotAll', NULL, '_Entire Graph', '<control>E', 'Plot entire graph', plot_entire_cb),
    list('PlotNeighb', NULL, '_Neighborhood Graph', '<control>N', 'Plot
         neighborhood of a single vertex', plot_neighb_cb),
    list('PlotComm', NULL, '_Community Graph', '<control>C', 'Plot a single
         community', plot_comm_cb))
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
  ui_manager$addUi(merge, '/menubar/plot', 'plotall', 'PlotAll', 'menuitem', FALSE)
  ui_manager$addUi(merge, '/menubar/plot', 'plotneighb', 'PlotNeighb', 'menuitem', FALSE)
  ui_manager$addUi(merge, '/menubar/plot', 'plotcomm', 'PlotComm', 'menuitem', FALSE)

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

  build.gui <- function(container) {  # "container" is "hboxMain"
  #=============================================================================
  # Create main (left side) vertical container
  #=============================================================================
  vboxMain <- gtkVBoxNew(F, 2)
  container$add(vboxMain)

  #---------------------------------------------------------
  # Create frame + vert container for *all* plot parameters
  #---------------------------------------------------------
  frame <- gtkFrameNew('Specify plotting parameters')
  vboxMain$add(frame)

  vbox <- gtkVBoxNew(F, 2)
  vbox$setBorderWidth(3)
  frame$add(vbox)

  # Create a frame for each plot area
  #---------------------------------------------------------
  frameG <- vboxG <- hboxOrient <- comboOrient <- graphname <- graphs <- hboxMin <-
    spinButtons <- hboxVsizeEqn <- vertSizeEqn <- slider <- hboxEwidthMax <-
    edgeWidthMax.spin <- graphics <- vboxPlot <- groupplot <- comboHemi <- vector('list', 2)
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
    # Both, single hemi, inter-hemi, homologous, inter/intra community/lobe?
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
               'Communities (weighted)', 'Class', 'Network')
  comboVcolor <- add_combo(hboxVcolor, choices, 'Vertex _color')
  gSignalConnect(comboVcolor, 'changed', function(widget, ...) {
      if (widget$getActive() %in% c(2, 5, 6)) {
        showLegend$setSensitive(T)
      } else {
        showLegend$setSensitive(F)
      }
  })

  #-----------------------------------------------------------------------------
  add_frame <- function(container, name, choices, const, max) {
    newFrame <- gtkFrameNew(name)
    container$add(newFrame)
    vbox <- gtkVBoxNew(F, 2)
    newFrame$add(vbox)

    hbox <- gtkHBoxNew(F, 4)
    vbox$packStart(hbox, F, F, 0)
    comboBox <- add_combo(hbox, choices, 'Scale by:')
    constEntry <- add_entry(hbox, char.width=3, entry.text=const)
    hboxOther <- gtkHBoxNew(F, 6)
    vbox$packStart(hboxOther, F, F, 0)
    otherEntry <- add_entry(hboxOther, label.text=paste0('\t', 'Other:'), char.width=10)
    otherEntry$setSensitive(F)
    for (i in 1:2) {
      hboxMin[[i]] <- gtkHBoxNew(F, 6)
      hboxOther$packStart(hboxMin[[i]], T, F, 0)
      labelMin <- gtkLabelNew(sprintf('Min. %i:', i))
      hboxMin[[i]]$packStart(labelMin, F, F, 0)
      spinButtons[[i]] <- gtkSpinButtonNewWithRange(min=0, max=10, step=1)
      hboxMin[[i]]$packStart(spinButtons[[i]], F, F, 0)
      spinButtons[[i]]$setSensitive(F)
    }
    return(list(vbox, comboBox, constEntry, otherEntry, spinButtons))
  }
  # Vertex size?
  #-----------------------------------------------------------------------------
  choices <- c('Constant', 'Degree', 'EV centrality', 'Btwn centrality',
               'K-core', 'Clustering coeff.', 'PC', 'GC', 'E.local', 'E.nodal',
               'Within-module degree z-score', 'Hub score', 'Vulnerability',
               'NN degree', 'Asymmetry', 'Eccentricity', 'Distance',
               'Distance strength', 'Lp', 'Other', 'Equation')
  if (is_weighted(graphs[[1]])) {
    choices <- c(choices, 'strength', 'knn.wt', 'E.local.wt', 'E.nodal.wt', 'S-core')
  }
  Vsize <- add_frame(vbox, 'Vertex size', choices, '5', max=10)

  # Create 2 entries for entering a more complicated eqn
  hboxVsizeEqnMain <- gtkHBoxNew(F, 6)
  Vsize[[1]]$packStart(hboxVsizeEqnMain, F, F, 0)
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
                  'k.core', 'transitivity', 'PC', 'GC', 'E.local', 'E.nodal',
                  'z.score', 'hub.score', 'vulnerability', 'knn', 'asymm',
                  'eccentricity', 'dist', 'dist.strength', 'Lp', 'other', 'eqn',
                  'strength', 'knn.wt', 'E.local.wt', 'E.nodal.wt', 's.core')
  vsize.measure <<- 'const'
  gSignalConnect(Vsize[[2]], 'changed', function(widget, ...) {
      vsize.measure <<- vsize.opts[widget$getActive() + 1]
      if (vsize.measure == 'const') {
        Vsize[[3]]$setSensitive(T)
        Vsize[[4]]$setSensitive(F)
        lapply(Vsize[[5]], function(x) x$setSensitive(F))
        lapply(vertSizeEqn, function(x) x$setSensitive(F))
      } else {
        Vsize[[3]]$setSensitive(F)
        Vsize[[4]]$setSensitive(F)
        lapply(Vsize[[5]], function(x) x$setSensitive(T))
        lapply(vertSizeEqn, function(x) x$setSensitive(F))
        if (vsize.measure == 'other') {
          Vsize[[4]]$setSensitive(T)
          lapply(vertSizeEqn, function(x) x$setSensitive(F))
          lapply(Vsize[[5]], function(x) x$setSensitive(T))
          for (j in seq_len(kNumGroups)) {
            gtkSpinButtonSetDigits(Vsize[[5]][[j]], 2)
            gtkSpinButtonSetIncrements(Vsize[[5]][[j]], step=0.01, page=0)
            gtkSpinButtonSetRange(Vsize[[5]][[j]], min=-100, max=100)
            gtkSpinButtonSetValue(Vsize[[5]][[j]], 0)
          }
        } else if (vsize.measure == 'eqn') {
          lapply(vertSizeEqn, function(x) x$setSensitive(T))
          lapply(Vsize[[5]], function(x) x$setSensitive(F))
        } else {
          if (vsize.measure == 'hub.score' && !is_directed(graphs[[1]])) vsize.measure <<- 'ev.cent'
          for (j in seq_len(kNumGroups)) {
            rangeX <- range(vertex_attr(graphs[[j]], vsize.measure), na.rm=TRUE)
            newMin <- rangeX[1]
            newMax <- rangeX[2]
            newStep <- ifelse(diff(rangeX) > 10 & newMin >= 0, 1,
                              ifelse(diff(rangeX) < 1, 0.01, 0.1))
            newDigits <- ifelse(diff(rangeX) > 10 & newMin >= 0, 0,
                                ifelse(diff(rangeX) < 1, 2, 1))
            gtkSpinButtonSetDigits(Vsize[[5]][[j]], newDigits)
            gtkSpinButtonSetIncrements(Vsize[[5]][[j]], step=newStep, page=0)
            gtkSpinButtonSetRange(Vsize[[5]][[j]], min=newMin, max=newMax)
            gtkSpinButtonSetValue(Vsize[[5]][[j]], newMin)
          }
        }
      }
  })
  # Commented out b/c clicks were causing spin value to change multiple times
  #gSignalConnect(vertSize.spin[[1]], 'value-changed', function(x) updateFun(1))
  #gSignalConnect(vertSize.spin[[2]], 'value-changed', function(x) updateFun(2))

  # Edge width?
  #-----------------------------------------------------------------------------
  choices <- c('Constant', 'Edge betweenness', 'Distance', 'Other')
  if ('weight' %in% edge_attr_names(graphs[[1]])) choices <- c(choices, 'Weight')
  Ewidth <- add_frame(vbox, 'Edge width', choices, '1', max=100)

  hboxEwidthOtherMax <- gtkHBoxNew(F, 6)
  Ewidth[[1]]$packStart(hboxEwidthOtherMax, F, F, 0)
  for (i in 1:2) {
    hboxEwidthMax[[i]] <- gtkHBoxNew(F, 6)
    if (i == 1) {
      hboxEwidthOtherMax$packStart(hboxEwidthMax[[i]], T, F, 0)
    } else {
      hboxEwidthOtherMax$packEnd(hboxEwidthMax[[i]], F, F, 0)
    }
    labelMax <- gtkLabelNew(sprintf('Max. %i:', i))
    hboxEwidthMax[[i]]$packStart(labelMax, F, F, 0)
    edgeWidthMax.spin[[i]] <- gtkSpinButtonNewWithRange(min=0, max=100, step=1)
    hboxEwidthMax[[i]]$packStart(edgeWidthMax.spin[[i]], F, F, 0)
    edgeWidthMax.spin[[i]]$setSensitive(F)
  }

  ewidth.opts <<- c('const', 'btwn', 'dist', 'other', 'weight')
  ewidth.measure <<- 'const'
  gSignalConnect(Ewidth[[2]], 'changed', function(widget, ...) {
    ewidth.measure <<- ewidth.opts[widget$getActive() + 1]
    if (ewidth.measure == 'const') {
      Ewidth[[3]]$setSensitive(T)
      Ewidth[[4]]$setSensitive(F)
      lapply(Ewidth[[5]], function(x) x$setSensitive(F))
      lapply(edgeWidthMax.spin, function(x) x$setSensitive(F))
    } else {
      Ewidth[[3]]$setSensitive(F)
      Ewidth[[4]]$setSensitive(F)
      lapply(Ewidth[[5]], function(x) x$setSensitive(T))
      lapply(edgeWidthMax.spin, function(x) x$setSensitive(T))
      if (ewidth.measure == 'other') {
        Ewidth[[4]]$setSensitive(T)
        lapply(Ewidth[[5]], function(x) x$setSensitive(T))
        lapply(edgeWidthMax.spin, function(x) x$setSensitive(T))
        for (j in seq_len(kNumGroups)) {
          gtkSpinButtonSetDigits(Ewidth[[5]][[j]], 2)
          gtkSpinButtonSetIncrements(Ewidth[[5]][[j]], step=0.01, page=0)
          gtkSpinButtonSetRange(Ewidth[[5]][[j]], min=-100, max=100)
          gtkSpinButtonSetValue(Ewidth[[5]][[j]], 0)
          gtkSpinButtonSetDigits(edgeWidthMax.spin[[j]], 2)
          gtkSpinButtonSetIncrements(edgeWidthMax.spin[[j]], step=0.01, page=0)
          gtkSpinButtonSetRange(edgeWidthMax.spin[[j]], min=-100, max=100)
          gtkSpinButtonSetValue(edgeWidthMax.spin[[j]], 100)
        }
      } else {
        for (j in seq_len(kNumGroups)) {
          rangeX <- range(edge_attr(graphs[[j]], ewidth.measure), na.rm=TRUE)
          newMin <- rangeX[1]
          newMax <- rangeX[2]
          newStep <- ifelse(diff(rangeX) > 1 & newMin >= 0, 1,
                            ifelse(diff(rangeX) < 1, 0.01, 0.1))
          newDigits <- ifelse(diff(rangeX) > 1 & newMin >= 0, 0, 2)
          gtkSpinButtonSetDigits(Ewidth[[5]][[j]], newDigits)
          gtkSpinButtonSetIncrements(Ewidth[[5]][[j]], step=newStep, page=0)
          gtkSpinButtonSetRange(Ewidth[[5]][[j]], min=newMin, max=newMax)
          gtkSpinButtonSetValue(Ewidth[[5]][[j]], newMin)
          gtkSpinButtonSetDigits(edgeWidthMax.spin[[j]], newDigits)
          gtkSpinButtonSetIncrements(edgeWidthMax.spin[[j]], step=newStep, page=0)
          gtkSpinButtonSetRange(edgeWidthMax.spin[[j]], min=newMin, max=newMax)
          gtkSpinButtonSetValue(edgeWidthMax.spin[[j]], newMax)
        }
      }
    }
  })
  #-----------------------------------------------------------------------------

  myLobe <<- 'All'
  if (identical(plotFunc, plot.brainGraph)) {
    # Highlight the diameter of each graph? Show edge set differences?
    hbox <- gtkHBoxNew(F, 6)
    vbox$packStart(hbox, F, F, 0)
    showDiameter <- add_check(hbox, 'Show _diameter?')
    edgeDiffs <- add_check(hbox, 'Show _edge differences?')
    if (kNumGroups <= 1) edgeDiffs$setSensitive(F)

    # Major lobe
    atlas.dt <- get(atlas)
    lobes <- data.frame(Lobe=c('All', as.character(atlas.dt[, levels(lobe)])))
    model <- rGtkDataFrame(lobes)
    setDT(lobes)
    lobes[, Lobe := as.character(Lobe)]
    view <<- gtkTreeView(model)
    view$getSelection()$setMode('multiple')
    column <- gtkTreeViewColumn('Lobe', gtkCellRendererText(), text=0)
    view$appendColumn(column)
    scrolled_window <- gtkScrolledWindow()
    scrolled_window$setSizeRequest(-1, 124)
    scrolled_window$add(view)
    vbox$add(scrolled_window)
  }

  #---------------------------------------------------------
  # Create 2 drawing areas for the plots
  #---------------------------------------------------------
  if (.Platform$OS.type == 'windows') {
    scr_width <- system("wmic desktopmonitor get screenwidth", intern=TRUE)
    scr_height <- system("wmic desktopmonitor get screenheight", intern=TRUE)
    res <- as.numeric(c(scr_width[-c(1, length(scr_width))],
                        scr_height[-c(1, length(scr_height))]))
  } else {
    display.size <- system("xrandr | grep '*' | awk '{print $1}'", intern=TRUE)
    res <- as.numeric(strsplit(display.size, 'x')[[1]])
  }

  ydim <- ifelse(res[2] < 800, 0.8 * res[2], 700)
  for (i in 1:2) {
    graphics[[i]] <- gtkDrawingArea()
    graphics[[i]]$setSizeRequest(res[1] / 3 - 10, ydim)
    vboxPlot[[i]] <- gtkVBox()
    vboxPlot[[i]]$packStart(graphics[[i]], expand=TRUE, fill=TRUE, padding=0)
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
  btnRename <- gtkButtonNewWithMnemonic('Pick _new graphs')
  the.buttons$packStart(btnOK, expand=TRUE, fill=F)
  the.buttons$packStart(btnRename, expand=TRUE, fill=F)
  container$showAll()

  # Callback functions
  updateFun <- function(j) {
    selection <- view$getSelection()
    sel_paths <- selection$getSelectedRows()$retval
    if (length(sel_paths) > 0) {
      sel_rows <- sapply(sel_paths, function(x) x$getIndices())
      if (identical(plotFunc, plot.brainGraph)) {
        myLobe <- paste(lobes[sel_rows+1, Lobe], collapse=', ')
      } else if (plotFunc == 'plot_community') {
        myComm <- comms[sel_rows+1, 'Community']
      } else if (plotFunc == 'plot_neighborhood') {
        myNeighb <- neighbs[sel_rows+1, 'name']
        myNeighbInd <- neighbs[sel_rows+1, 'index']
      }
    }
    k <- setdiff(seq_len(kNumGroups), j)
    update_brainGraph_gui(plotDev=groupplot[[j]], g=graphs[[j]], g2=graphs[[k]],
      plotFunc=plotFunc, vsize.measure=vsize.measure, ewidth.measure=ewidth.measure,
      vertColor=comboVcolor, hemi=comboHemi[[j]], lobe=myLobe, orient=comboOrient[[j]],
      vertSize.min=Vsize[[5]][[j]], edgeWidth.min=Ewidth[[5]][[j]],
      edgeWidth.max=edgeWidthMax.spin[[j]], vertSize.const=Vsize[[3]],
      edgeWidth.const=Ewidth[[3]], vertLabels=vertLabels, showLegend=showLegend,
      comms=myComm, neighb=myNeighb, neighbInd=myNeighbInd, slider=slider[[j]],
      Vsize[[4]], Ewidth[[4]], vertSize.eqn=vertSizeEqn[[j]],
      showDiameter=showDiameter, edgeDiffs=edgeDiffs)
  }
  gSignalConnect(btnOK, 'clicked', function(widget) {
                 for (i in seq_len(kNumGroups)) updateFun(i)})
  gSignalConnect(btnRename, 'clicked', function(widget) set.names())

  # Horizontal slider for edge curvature in circle plots
  orient_cb <- function(widget, ind) {
    if (widget$getActive() == 3) {
      comboHemi[[ind]]$setSensitive(T)
      slider[[ind]] <<- gtkHScaleNewWithRange(min=-1, max=1, step=0.05)
      vboxPlot[[ind]]$packStart(slider[[ind]], F, F, 0)
      slider[[ind]]$setValue(0.25)
      gSignalConnect(slider[[ind]], 'value-changed', function(widget) updateFun(ind))
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
  for (i in 1:2) gSignalConnect(comboOrient[[i]], 'changed', orient_cb, i)

  list(plot1=groupplot[[1]], plot2=groupplot[[2]], graphname=graphname,
       comboVcolor=comboVcolor, comboHemi=comboHemi, vbox=vbox)
  }
}
