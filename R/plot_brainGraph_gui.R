#' GUI for plotting graphs overlaid on an MNI152 image or in a circle
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
  if (!all(vapply(c('RGtk2', 'cairoDevice'), requireNamespace, logical(1L), quietly=TRUE))) {
    warning('You must install RGtk2 and cairoDevice to use the brainGraph GUI.')
    return(NULL)
  }

  window <- RGtk2::gtkWindow('toplevel')
  window$setTitle('brainGraph')
  window$setIconFromFile(system.file('extdata', 'brainGraph_icon.png', package='brainGraph'))
  action_group <- RGtk2::gtkActionGroup('brainGraphActions')

  #=============================================================================
  # Add a menubar
  #=============================================================================
  # Action group for the "File" menu
  #-----------------------------------------------
  file.actions <- list(
    list('FileMenu', NULL, '_File'),
    list('SelectGraphs', 'gtk-open', 'Select graphs', '<control>O', 'Select graphs', select_graphs),
    list('Save1', 'gtk-save', 'Save plot _1', '<control>1', 'Save plot 1?', save_cb1),
    list('Save2', 'gtk-save', 'Save plot _2', '<control>2', 'Save plot 2?', save_cb2),
    list('Quit', 'gtk-quit', '_Quit', '<control>Q', 'Quit the application', quit_cb)
  )
  action_group$addActions(file.actions, window)

  # Action group for the "Plot" menu
  #-----------------------------------------------
  plot.actions <- list(  # Name, ID, label, accelerator, tooltip, callback func
    list('PlotMenu', NULL, '_Plot'),
    list('PlotAll', NULL, '_Entire Graph', '<control>E', 'Plot entire graph',
         function(widget, win) build_gui(widget, win, 'lobe')),
    list('PlotNeighb', NULL, '_Neighborhood Graph', '<control>N',
         'Plot neighborhood of a single vertex',
         function(widget, win) build_gui(widget, win, 'neighborhood')),
    list('PlotComm', NULL, '_Community Graph', '<control>C', 'Plot by communities',
         function(widget, win) build_gui(widget, win, 'comm')),
    list('PlotNetwork', NULL, '_Functional networks (Dosenbach, Power, Gordon)', '<control>F',
         'Plot functional networks', function(widget, win) build_gui(widget, win, 'network')),
    list('PlotArea', NULL, 'Cortical _areas (HCP)', '<control>A',
         'Plot HCP cortical areas', function(widget, win) build_gui(widget, win, 'area')),
    list('PlotYeo7', NULL, '_Yeo 7 networks (brainnetome)', '<control>Y',
         'Plot Yeo 7 networks', function(widget, win) build_gui(widget, win, 'Yeo_7network')),
    list('PlotYeo17', NULL, 'Yeo 1_7 networks (brainnetome)', '<control>7',
         'Plot Yeo 17 networks', function(widget, win) build_gui(widget, win, 'Yeo_17network')),
    list('PlotGyri', NULL, '_Gyri (brainnetome)', '<control>G',
         'Plot Gyri', function(widget, win) build_gui(widget, win, 'gyrus')))
  action_group$addActions(plot.actions, window)

  #---------------------------------------------------------
  # Create a GtkUIManager instance
  #---------------------------------------------------------
  ui_manager <- RGtk2::gtkUIManager()
  ui_manager$insertActionGroup(action_group, pos=0)

  merge <- ui_manager$newMergeId()
  ui_manager$addUi(merge.id=merge, path='/', name='menubar', action=NULL,
                   type='menubar', top=FALSE)
  ui_manager$addUi(merge, '/menubar', 'file', 'FileMenu', 'menu', FALSE)
  ui_manager$addUi(merge, '/menubar/file', 'sel', 'SelectGraphs', 'menuitem', FALSE)
  ui_manager$addUi(merge, '/menubar/file', 'save1', 'Save1', 'menuitem', FALSE)
  ui_manager$addUi(merge, '/menubar/file', 'save2', 'Save2', 'menuitem', FALSE)
  ui_manager$addUi(merge, '/menubar/file', 'sep', action=NULL, type=NULL, FALSE)
  ui_manager$addUi(merge, '/menubar/file', 'quit', 'Quit', 'menuitem', FALSE)

  ui_manager$addUi(merge, '/menubar', 'plot', 'PlotMenu', 'menu', FALSE)
  ui_manager$addUi(merge, '/menubar/plot', 'plotall', 'PlotAll', 'menuitem', FALSE)
  ui_manager$addUi(merge, '/menubar/plot', 'plotneighb', 'PlotNeighb', 'menuitem', FALSE)
  ui_manager$addUi(merge, '/menubar/plot', 'plotcomm', 'PlotComm', 'menuitem', FALSE)
  ui_manager$addUi(merge, '/menubar/plot', 'plotnetwork', 'PlotNetwork', 'menuitem', FALSE)
  ui_manager$addUi(merge, '/menubar/plot', 'plotarea', 'PlotArea', 'menuitem', FALSE)
  ui_manager$addUi(merge, '/menubar/plot', 'plotyeo7', 'PlotYeo7', 'menuitem', FALSE)
  ui_manager$addUi(merge, '/menubar/plot', 'plotyeo17', 'PlotYeo17', 'menuitem', FALSE)
  ui_manager$addUi(merge, '/menubar/plot', 'plotgyri', 'PlotGyri', 'menuitem', FALSE)

  menubar <- ui_manager$getWidget('/menubar')
  window$addAccelGroup(ui_manager$getAccelGroup())
  vboxMainMenu <- RGtk2::gtkVBoxNew(homogeneous=FALSE, spacing=2)
  window$add(vboxMainMenu)
  vboxMainMenu$packStart(menubar, expand=FALSE)

  # Select graph objects to plot
  select_graphs(NULL, window)

  #=============================================================================
  # Main GUI-building function
  #=============================================================================
  build_gui <- function(widget, window, vGroupAttr) {  # "widget" is unused
    kill_others(window[[1L]])
    container <- RGtk2::gtkHBoxNew(homogeneous=FALSE, spacing=8)
    window[[1L]]$add(container)

    # Defaults / global variables
    vsizeName <- ewidthName <- 'const'
    vcols <- c('None (lightblue)', 'Communities', 'Lobes', 'Components',
               'Communities (weighted)', 'Class', 'Network', 'Cortical Area (HCP)',
               'Yeo_7network (brainnetome)', 'Yeo_17network (brainnetome)')
    if (vGroupAttr == 'neighborhood') vcols <- c(vcols, 'Neighborhoods')
    vsizes <- c('Constant', 'Degree', 'EV Centrality', 'Btwn. Centrality',
                'K-core', 'Clustering Coeff.', 'Participation Coeff.', 'Gateway Coeff.',
                'Local Eff.', 'Nodal Eff.', 'Within-module degree z-score',
                'Vulnerability', 'K-Nearest-Neighb. Degree', 'Asymmetry', 'Eccentricity',
                'Distance', 'Distance Strength', 'Shortest Path Length', 'Other', 'Equation')
    vsizeAttr <- c('const', 'degree', 'ev.cent', 'btwn.cent',
                   'k.core', 'transitivity', 'PC', 'GC', 'E.local', 'E.nodal',
                   'z.score', 'vulnerability', 'knn', 'asymm',
                   'eccentricity', 'dist', 'dist.strength', 'Lp', 'other', 'eqn',
                   'strength', 'knn.wt', 'E.local.wt', 'E.nodal.wt', 's.core')
    ewidths <- c('Constant', 'Edge Betweenness', 'Distance', 'Other')
    ewidthAttr <- c('const', 'btwn', 'dist', 'other', 'weight')
    hboxVsizeEqn <- vsizeEqn <- slider <- hboxEwidthMax <- ewidthMaxSpin <- graphics <-
      vboxPlot <- vector('list', 2L)
    plotdev <<- vector('list', 2L)

    # Create the main (L side) vertical container
    # And create a frame + vert box for *all* plot params
    #=================================================================
    vboxMain <- RGtk2::gtkVBoxNew(homogeneous=FALSE, spacing=2)
    container$packStart(vboxMain)
    frameMain <- RGtk2::gtkFrameNew(label='Specify plotting parameters')
    vboxMain$packStart(frameMain)
    vbox <- RGtk2::gtkVBoxNew(homogeneous=FALSE, spacing=2)
    vbox$setBorderWidth(5)
    frameMain$add(vbox)

    # Create a frame for each graph selected
    #-------------------------------------------------------
    kNumGroups <- sum(vapply(graphObj, function(x) nchar(x) > 0L, logical(1L)))
    groupSeq <- seq_len(kNumGroups)
    res <- add_orient_hemi(vbox, graphObj, kNumGroups, vGroupAttr)
    graphs <- res$graphs; comboOrient <- res$orient; comboHemi <- res$hemi

    # Frame for the annotation-related options
    #---------------------------------------------------------------------------
    frameAnn <- RGtk2::gtkFrameNew(label='Annotations')
    vbox$packStart(frameAnn)
    vboxAnn <- RGtk2::gtkVBoxNew(homogeneous=FALSE, spacing=2)
    vboxAnn$setBorderWidth(5)
    frameAnn$add(vboxAnn)
    hboxChecks <- RGtk2::gtkHBoxNew(homogeneous=TRUE, spacing=6)
    vboxAnn$packStart(hboxChecks)

    # Show a title and/or subtitle?
    vboxTitle <- RGtk2::gtkVBoxNew(homogeneous=FALSE, spacing=2)
    mainTitle <- add_check(vboxTitle, 'Display _main title?')
    subTitle <- add_check(vboxTitle, 'Display _subtitle?')
    hboxChecks$packStart(vboxTitle, expand=FALSE)

    # Should vertex labels and/or a legend displayed?
    vboxLabels <- RGtk2::gtkVBoxNew(homogeneous=FALSE, spacing=2)
    vertLabels <- add_check(vboxLabels, 'Display vertex _labels?')
    showLegend <- add_check(vboxLabels, 'Show le_gend?')
    showLegend$setSensitive(FALSE)
    hboxChecks$packStart(vboxLabels, expand=FALSE)

    # Vertex colors?
    comboVcol <- add_combo(vboxAnn, vcols, 'Vertex _color', spacing=4)

    # Vertex size?
    #---------------------------------------------------------------------------
    if (is_weighted(graphs[[1L]])) {
      vsizes <- c(vsizes, 'Strength', 'K-Nearest-Neighb. Degree (wt)', 'Local Eff. (wt)',
                  'Nodal Eff. (wt)', 'S-core')
    }
    Vsize <- add_frame(vbox, 'Vertex size', vsizes, '5', kNumGroups=kNumGroups, max=10)

    # Create 1-2 entries for entering a more complicated eqn
    hboxVsizeEqnMain <- RGtk2::gtkHBoxNew(homogeneous=FALSE, spacing=6)
    Vsize$box$packStart(hboxVsizeEqnMain, expand=FALSE)
    vboxVsizeEqnMain <- RGtk2::gtkVBoxNew(homogeneous=FALSE, spacing=2)
    hboxVsizeEqnMain$packStart(vboxVsizeEqnMain, expand=FALSE)
    for (i in groupSeq) {
      hboxVsizeEqn[[i]] <- RGtk2::gtkHBoxNew(homogeneous=FALSE, spacing=6)
      vboxVsizeEqnMain$packStart(hboxVsizeEqn[[i]], expand=FALSE)
      vsizeEqn[[i]] <- add_entry(hboxVsizeEqn[[i]], label=sprintf('Eqn. %i', i), nchars=40)
      vsizeEqn[[i]]$setSensitive(FALSE)
    }

    vertSizeHelp <- RGtk2::gtkButtonNewWithLabel('?')
    hboxVsizeEqnMain$packStart(vertSizeHelp, expand=FALSE)
    RGtk2::gSignalConnect(vertSizeHelp, signal='clicked', f=help_win)

    RGtk2::gSignalConnect(Vsize$combo, signal='changed', function(widget, ...) {
        vsizeName <<- vsizeAttr[widget$getActive() + 1]
        if (vsizeName == 'const') {
          Vsize$const$setSensitive(TRUE)
          Vsize$other$setSensitive(FALSE)
          lapply(Vsize$spin[groupSeq], function(x) x$setSensitive(FALSE))
          lapply(vsizeEqn[groupSeq], function(x) x$setSensitive(FALSE))
        } else {
          Vsize$const$setSensitive(FALSE)
          Vsize$other$setSensitive(FALSE)
          lapply(Vsize$spin[groupSeq], function(x) x$setSensitive(TRUE))
          lapply(vsizeEqn[groupSeq], function(x) x$setSensitive(FALSE))

          if (vsizeName == 'other') {
            Vsize$other$setSensitive(TRUE)
            for (j in groupSeq) {
              Vsize$spin[[j]]$setDigits(2)
              Vsize$spin[[j]]$setIncrements(step=0.01, page=0)
              Vsize$spin[[j]]$setRange(min=-100, max=100)
              Vsize$spin[[j]]$setValue(0)
            }
          } else if (vsizeName == 'eqn') {
            lapply(vsizeEqn[groupSeq], function(x) x$setSensitive(TRUE))
            lapply(Vsize$spin[groupSeq], function(x) x$setSensitive(FALSE))
          } else {  # Everything but "const", "other", and "eqn"
            for (j in groupSeq) {
              rangeX <- range(vertex_attr(graphs[[j]], vsizeName), na.rm=TRUE)
              vrange <- diff(rangeX)
              newMin <- rangeX[1L]
              newStep <- ifelse(vrange > 10 && newMin >= 0, 1,
                                ifelse(vrange < 1, 0.01, 0.1))
              newDigits <- ifelse(vrange > 10 && newMin >= 0, 0,
                                  ifelse(vrange < 1, 2, 1))
              Vsize$spin[[j]]$setDigits(newDigits)
              Vsize$spin[[j]]$setIncrements(step=newStep, page=0)
              Vsize$spin[[j]]$setRange(min=newMin, max=rangeX[2L])
              Vsize$spin[[j]]$setValue(newMin)
            }
          }
        }
    })

    # Edge width?
    #-----------------------------------------------------------------------------
    if (is_weighted(graphs[[1L]])) ewidths <- c(ewidths, 'Weight')
    Ewidth <- add_frame(vbox, 'Edge width', ewidths, '1', kNumGroups=kNumGroups, max=100)

    hboxEwidthOtherMax <- RGtk2::gtkHBoxNew(homogeneous=FALSE, spacing=6)
    Ewidth$box$packStart(hboxEwidthOtherMax, expand=FALSE)
    for (i in groupSeq) {
      hboxEwidthMax[[i]] <- RGtk2::gtkHBoxNew(homogeneous=FALSE, spacing=4)
      hboxEwidthOtherMax$packStart(hboxEwidthMax[[i]], expand=TRUE, fill=FALSE)
      labelMax <- RGtk2::gtkLabelNew(sprintf('Max. %i:', i))
      hboxEwidthMax[[i]]$packStart(labelMax, expand=FALSE)
      ewidthMaxSpin[[i]] <- RGtk2::gtkSpinButtonNewWithRange(min=0, max=100, step=1)
      hboxEwidthMax[[i]]$packStart(ewidthMaxSpin[[i]], expand=FALSE)
      ewidthMaxSpin[[i]]$setSensitive(FALSE)
    }

    RGtk2::gSignalConnect(Ewidth$combo, signal='changed', function(widget, ...) {
      ewidthName <<- ewidthAttr[widget$getActive() + 1]
      if (ewidthName == 'const') {
        Ewidth$const$setSensitive(TRUE)
        Ewidth$other$setSensitive(FALSE)
        lapply(Ewidth$spin[groupSeq], function(x) x$setSensitive(FALSE))
        lapply(ewidthMaxSpin[groupSeq], function(x) x$setSensitive(FALSE))
      } else {
        Ewidth$const$setSensitive(FALSE)
        lapply(Ewidth$spin[groupSeq], function(x) x$setSensitive(TRUE))
        lapply(ewidthMaxSpin[groupSeq], function(x) x$setSensitive(TRUE))
        if (ewidthName == 'other') {
          Ewidth$other$setSensitive(TRUE)
          newDigits <- rep.int(2, kNumGroups)
          newStep <- rep_len(0.01, kNumGroups)
          newMin <- rep_len(-100, kNumGroups); newMax <- -newMin
        } else {
          Ewidth$other$setSensitive(FALSE)
          newDigits <- newStep <- newMin <- newMax <- rep(0, kNumGroups)
          for (j in groupSeq) {
            rangeX <- range(edge_attr(graphs[[j]], ewidthName), na.rm=TRUE)
            erange <- diff(rangeX)
            newMin[j] <- rangeX[1L]; newMax[j] <- rangeX[2L]
            newStep[j] <- ifelse(erange > 10 && newMin >= 0, 1,
                                 ifelse(erange < 1, 0.01, 0.1))
            newDigits[j] <- ifelse(erange > 1 && newMin >= 0, 0, 2)
          }
        }
        for (j in groupSeq) {
            Ewidth$spin[[j]]$setDigits(newDigits[j])
            Ewidth$spin[[j]]$setIncrements(step=newStep[j], page=0)
            Ewidth$spin[[j]]$setRange(min=newMin[j], max=newMax[j])
            Ewidth$spin[[j]]$setValue(newMin[j])
            ewidthMaxSpin[[j]]$setDigits(newDigits[j])
            ewidthMaxSpin[[j]]$setIncrements(step=newStep[j], page=0)
            ewidthMaxSpin[[j]]$setRange(min=newMin[j], max=newMax[j])
            ewidthMaxSpin[[j]]$setValue(newMax[j])
        }
      }
    })

    # For "lobe", "show diameter" and "edge diffs"
    #-------------------------------------------------------
    showDiameter <- edgeDiffs <- NULL
    if (vGroupAttr == 'lobe') {
      hbox <- RGtk2::gtkHBoxNew(homogeneous=TRUE, spacing=6)
      vbox$packStart(hbox, expand=FALSE)
      showDiameter <- add_check(hbox, 'Show _diameter?')
      edgeDiffs <- add_check(hbox, 'Show _edge differences?')
      if (kNumGroups < 2L) edgeDiffs$setSensitive(FALSE)
    }

    # Add a scrolling window and the buttons
    view <- add_scroll_window(graphs, vGroupAttr, kNumGroups, vbox)
    buttons <- add_buttons(vboxMain)

    #-------------------------------------------------------
    # Create 2 drawing areas for the plots
    #-------------------------------------------------------
    scr_res <- get_screen_size()
    for (i in groupSeq) {
      graphics[[i]] <- RGtk2::gtkDrawingArea()
      graphics[[i]]$setSizeRequest(scr_res[1L], scr_res[2L])
      vboxPlot[[i]] <- RGtk2::gtkVBox()
      vboxPlot[[i]]$packStart(graphics[[i]], expand=TRUE, fill=TRUE)
      container$packStart(vboxPlot[[i]])
      cairoDevice::asCairoDevice(graphics[[i]])
      par(pty='s', mar=rep(0, 4))
      plotdev[[i]] <<- dev.cur()
    }

    #-------------------------------------------------------
    # Main callback function to update the plot areas
    #-------------------------------------------------------
    updateFun <- function(index) {
      sel_paths <- view$getSelection()$getSelectedRows()$retval
      sel_rows <- vapply(sel_paths, function(x) x$getIndices(), integer(1L)) + 1L
      mod <- view$getModel()
      vgroups <- mod[sel_rows, dim(mod)[2L]]
      if (!vGroupAttr %in% c('neighborhood', 'comm')) vgroups <- paste(vgroups, collapse=', ')

      k <- setdiff(groupSeq, index)
      dev.set(plotdev[[index]])
      if (vsizeName == 'other') vsizeName <- Vsize$other$getText()
      if (ewidthName == 'other') ewidthName <- Ewidth$other$getText()
      update_brainGraph_gui(graphs[[index]], comboOrient[[index]], comboHemi[[index]],
                            vertLabels, showLegend, mainTitle, subTitle, comboVcol,
                            vsizeName, Vsize$const, Vsize$spin[[index]], vsizeEqn[[index]],
                            ewidthName, Ewidth$const, Ewidth$spin[[index]], ewidthMaxSpin[[index]],
                            vGroupAttr, vgroups,
                            slider[[index]], showDiameter, edgeDiffs, g2=graphs[[k]])
    }
    RGtk2::gSignalConnect(buttons[[1L]], signal='clicked',
                          function(widget) for (i in groupSeq) updateFun(i))
    RGtk2::gSignalConnect(buttons[[2L]], signal='clicked', select_graphs, window)

    # Horizontal slider for edge curvature in circle plots
    orient_cb <- function(widget, i) {
      if (widget$getActive() == 3) {
        comboHemi[[i]]$setSensitive(TRUE)
        slider[[i]] <<- RGtk2::gtkHScaleNewWithRange(min=-1, max=1, step=0.05)
        vboxPlot[[i]]$packStart(slider[[i]], expand=FALSE)
        slider[[i]]$setValue(0.25)
        updateFun(i)
        RGtk2::gSignalConnect(slider[[i]], signal='value-changed', function(widget) updateFun(i))
      } else {
        kill_others(vboxPlot[[i]])
        if (widget$getActive() %in% c(1, 2)) {
          comboHemi[[i]]$setActive(widget$getActive())
          comboHemi[[i]]$setSensitive(FALSE)
        } else {
          comboHemi[[i]]$setSensitive(TRUE)
          updateFun(i)
        }
      }
    }
    for (i in groupSeq) {
      RGtk2::gSignalConnect(comboOrient[[i]], signal='changed', orient_cb, i)
      RGtk2::gSignalConnect(comboHemi[[i]], signal='changed', function(widget, i) updateFun(i), i)
    }
    RGtk2::gSignalConnect(comboVcol, signal='changed', function(widget) {
        showLegend$setSensitive(TRUE)
        if (widget$getActive() %in% c(0, 1, 3, 4)) showLegend$setSensitive(FALSE)
        for (i in groupSeq) updateFun(i)
    })
    comboNum <- switch(vGroupAttr, neighborhood=0, comm=1, lobe=, gyrus=2, network=6,
                       area=7, Yeo_7network=8, Yeo_17network=9)
    comboVcol$setActive(comboNum)

    # Signals for the check buttons
    RGtk2::gSignalConnect(mainTitle, signal='toggled',
                          function(w) for (i in groupSeq) updateFun(i))
    RGtk2::gSignalConnect(subTitle, signal='toggled',
                          function(w) for (i in groupSeq) updateFun(i))
    RGtk2::gSignalConnect(vertLabels, signal='toggled',
                          function(w) for (i in groupSeq) updateFun(i))
    RGtk2::gSignalConnect(showLegend, signal='toggled',
                          function(w) for (i in groupSeq) updateFun(i))

    container$showAll()
  }
}

#===============================================================================
# Helper functions for adding/updating the GUI elements
#===============================================================================
# The graph selector window
select_graphs <- function(widget, window) {
  dialog <- RGtk2::gtkDialogNewWithButtons(title='Choose graph objects',
                                           parent=window, flags='destroy-with-parent',
                                           'gtk-ok', RGtk2::GtkResponseType['ok'],
                                           'gtk-cancel', RGtk2::GtkResponseType['cancel'],
                                           show=FALSE)
  tempVbox <- dialog$getContentArea()
  graphObjEntry <- tempHbox <- vector('list', length=2)
  for (i in 1:2) {
    tempHbox[[i]] <- RGtk2::gtkHBoxNew(homogeneous=FALSE, spacing=8)
    tempVbox$packStart(tempHbox[[i]], expand=FALSE)
    graphObjEntry[[i]] <- add_entry(tempHbox[[i]], label=paste0('Graph ', i, ' name:'), nchars=20)
  }
  RGtk2::gSignalConnect(dialog, signal='response', function(dialog, response) {
      if (response == RGtk2::GtkResponseType['ok']) {
        graphObj[[1L]] <<- graphObjEntry[[1L]]$getText()
        graphObj[[2L]] <<- graphObjEntry[[2L]]$getText()

        if (nchar(graphObj[[1L]]) > 0L) {
          if (!is_igraph(eval(parse(text=graphObj[[1L]])))) {
            warnDialog <- RGtk2::gtkMessageDialog(parent=dialog, flags='destroy-with-parent',
                                                  type='error', buttons='close',
                                                  'Error: Not an igraph object!')
            response <- warnDialog$run()
          }
          if (response == RGtk2::GtkResponseType['close']) warnDialog$destroy()
        } else {
          if (nchar(graphObj[[2L]]) > 0L) {
            if (!is_igraph(eval(parse(text=graphObj[[2L]])))) {
              warnDialog <- RGtk2::gtkMessageDialog(parent=dialog, flags='destroy-with-parent',
                                                    type='error', buttons='close',
                                                    'Error: Not an igraph object!')
              response <- warnDialog$run()
              if (response == RGtk2::GtkResponseType['close']) warnDialog$destroy()
            }
          }
        }
        dialog$destroy()
      }
  })
  dialog$showAll()
  dialog$setModal(TRUE)
}

# Add frames for orientation + hemi, and extract graph objects
#-----------------------------------------------------------
add_orient_hemi <- function(container, graphObj, kNumGroups, vGroupAttr) {
  orient <- c('Axial', 'Sagittal (left)', 'Sagittal (right)', 'Circular')
  hemis <- c('Both hemispheres',
             paste(c('Left', 'Right', 'Interhemispheric', 'Homologous', 'Inter-community',
                     'Intra-community', 'Inter-lobar', 'Intra-lobar', 'Inter-network',
                     'Intra-network', 'Inter-area', 'Intra-area',
                     'Inter-Yeo7', 'Intra-Yeo7', 'Inter-Yeo17', 'Intra-Yeo17'), 'only'))
  frame <- vbox <- comboOrient <- myGs <- comboHemi <- vector('list', kNumGroups)
  for (i in seq_len(kNumGroups)) {
    frame[[i]] <- RGtk2::gtkFrameNew(label=sprintf('Graph object: "%s"', graphObj[[i]]))
    container$packStart(frame[[i]])

    vbox[[i]] <- RGtk2::gtkVBoxNew(homogeneous=FALSE, spacing=2)
    vbox[[i]]$setBorderWidth(5)
    frame[[i]]$add(vbox[[i]])

    comboOrient[[i]] <- add_combo(vbox[[i]], orient, 'Orientation', spacing=6)
    if (nchar(graphObj[[i]]) > 0L) myGs[[i]] <- eval(parse(text=graphObj[[i]]))

    comboHemi[[i]] <- add_combo(vbox[[i]], hemis, 'Hemi/edges', spacing=6)
    if (vGroupAttr == 'neighborhood') comboHemi[[i]]$setSensitive(FALSE)
  }
  list(graphs=myGs, orient=comboOrient, hemi=comboHemi)
}

# Function to add an entry
#-----------------------------------------------------------
add_entry <- function(container, label=NULL, nchars, entry.text=NULL) {
  if (!is.null(label)) {
    label <- RGtk2::gtkLabelNewWithMnemonic(label)
    container$packStart(label, expand=FALSE)
  }
  entry <- RGtk2::gtkEntryNew()
  entry$setWidthChars(nchars)
  if (!is.null(entry.text)) entry$setText(entry.text)
  container$packStart(entry, expand=FALSE)
  return(entry)
}

# Function to add a check button & label
#-----------------------------------------------------------
add_check <- function(container, label.text) {
  chk.button <- RGtk2::gtkCheckButtonNewWithMnemonic(label=label.text)
  container$packStart(chk.button, expand=TRUE)
  return(chk.button)
}

# Create a "gtkHBox" and add a combobox to it
# Set "hbox=TRUE" to also return the "gtkHBox" (e.g., to add more to it)
#-----------------------------------------------------------
add_combo <- function(container, choices, label.text, spacing, retHbox=FALSE) {
  hbox <- RGtk2::gtkHBoxNew(homogeneous=FALSE, spacing=spacing)
  container$packStart(hbox, expand=FALSE)
  combo <- RGtk2::gtkComboBoxNewText()
  combo$show()
  for (choice in choices) combo$appendText(choice)
  combo$setActive(0)
  label <- RGtk2::gtkLabelNewWithMnemonic(label.text)
  hbox$packStart(label, expand=FALSE)
  hbox$packStart(combo)
  res <- if (isTRUE(retHbox)) list(combo=combo, hbox=hbox) else combo
  res
}

#---------------------------------------------------------
# Callback functions for the "File" menu
#---------------------------------------------------------
save_cb <- function(widget, window, plot.dev) {
  dialog <- RGtk2::gtkFileChooserDialog('Enter a filename', window, 'save',
                                        'gtk-cancel', RGtk2::GtkResponseType['cancel'],
                                        'gtk-save', RGtk2::GtkResponseType['accept'])
  RGtk2::gSignalConnect(dialog, signal='response', function(dialog, response) {
                          if (response == RGtk2::GtkResponseType['accept']) {
                            fname <- dialog$getFilename()
                            dev.set(plot.dev)
                            dev.copy(png, fname)
                            dev.off()
                          }
                          dialog$destroy()})
}

save_cb1 <- function(widget, win) save_cb(widget, win, plot.dev=plotdev[[1L]])
save_cb2 <- function(widget, win) save_cb(widget, win, plot.dev=plotdev[[2L]])
quit_cb <- function(widget, win) win$destroy()

# Add a frame for "vertex size" and "edge width"
#-----------------------------------------------------------
add_frame <- function(container, name, choices, const, kNumGroups, max) {
  newFrame <- RGtk2::gtkFrameNew(label=name)
  container$packStart(newFrame)
  vbox <- RGtk2::gtkVBoxNew(homogeneous=FALSE, spacing=2)
  vbox$setBorderWidth(5)
  newFrame$add(vbox)

  comboBox <- add_combo(vbox, choices, 'Scale by:', spacing=4, retHbox=TRUE)
  constEntry <- add_entry(comboBox$hbox, nchars=3, entry.text=const)
  hboxOther <- RGtk2::gtkHBoxNew(homogeneous=FALSE, spacing=6)
  vbox$packStart(hboxOther, expand=FALSE)
  otherEntry <- add_entry(hboxOther, label=paste0('\t', 'Other:'), nchars=10)
  otherEntry$setSensitive(FALSE)
  hboxMin <- spinButtons <- vector('list', 2L)
  for (i in seq_len(kNumGroups)) {
    hboxMin[[i]] <- RGtk2::gtkHBoxNew(homogeneous=FALSE, spacing=6)
    hboxOther$packStart(hboxMin[[i]], TRUE, fill=FALSE)
    labelMin <- RGtk2::gtkLabelNew(sprintf('Min. %i:', i))
    hboxMin[[i]]$packStart(labelMin, expand=FALSE)
    spinButtons[[i]] <- RGtk2::gtkSpinButtonNewWithRange(min=0, max=max, step=1)
    hboxMin[[i]]$packStart(spinButtons[[i]], expand=FALSE)
    spinButtons[[i]]$setSensitive(FALSE)
  }
  return(list(box=vbox, combo=comboBox$combo, const=constEntry, other=otherEntry, spin=spinButtons))
}

# Add a scroll window at bottom, for "lobes", "neighbors", "comms", etc.
#-----------------------------------------------------------
add_scroll_window <- function(graphs, vGroupAttr, kNumGroups, vbox) {
  g1 <- graphs[[1L]]
  if (vGroupAttr == 'neighborhood') {
    columns <- list(RGtk2::gtkTreeViewColumn('Index', RGtk2::gtkCellRendererText(), text=0),
                    RGtk2::gtkTreeViewColumn('Vertex Name', RGtk2::gtkCellRendererText(), text=1))
    DF <- data.frame(index=seq_along(V(g1)), name=V(g1)$name, stringsAsFactors=FALSE)
  } else {
    colname <- switch(vGroupAttr, lobe='Lobe', comm='Community', area='Cortical Area',
                      network='Functional network', Yeo_7network='Network (Yeo 7)',
                      Yeo_17network='Network (Yeo 17)', gyrus='Gyrus')
    columns <- list(RGtk2::gtkTreeViewColumn(colname, RGtk2::gtkCellRendererText(), text=0))
    if (vGroupAttr == 'comm') {
      all.comms <- which(table(V(g1)$comm) > 2L)
      if (kNumGroups > 1L) {
        c2 <- which(table(V(graphs[[2L]])$comm) > 2L)
        all.comms <- union(all.comms, c2)
      }
      DF <- data.frame(Data=all.comms)
    } else {
      DF <- data.frame(Data=sort(unique(vertex_attr(g1, vGroupAttr))))
    }
    if (vGroupAttr == 'lobe') DF <- rbind(data.frame(Data='All'), DF)
  }
  model <- RGtk2::rGtkDataFrame(DF)
  treeview <- RGtk2::gtkTreeView(model)
  treeview$getSelection()$setMode('multiple')
  lapply(columns, function(x) treeview$appendColumn(x))
  scrolled_window <- RGtk2::gtkScrolledWindow()
  scrolled_window$setSizeRequest(width=-1, height=150)
  scrolled_window$add(treeview)
  vbox$packStart(scrolled_window)
  treeview$setCursor(RGtk2::gtkTreePathNewFromIndices(0))  # Default value
  treeview
}

# Add the "OK" and "New graphs" buttons
#-----------------------------------------------------------
add_buttons <- function(container) {
  buttons <- RGtk2::gtkHButtonBoxNew()
  buttons$setBorderWidth(5)
  container$packStart(buttons)
  btnOK <- RGtk2::gtkButtonNewFromStock('gtk-ok')
  btnRename <- RGtk2::gtkButtonNewWithMnemonic('Pick _new graphs')
  buttons$packStart(btnOK, expand=TRUE, fill=FALSE)
  buttons$packStart(btnRename, expand=TRUE, fill=FALSE)
  list(btnOK, btnRename)
}

# Create a "help" window for vertex size equations
#-----------------------------------------------------------
help_win <- function(widget) {
    helpWin <- RGtk2::gtkWindow()
    helpVbox <- RGtk2::gtkVBoxNew(homogeneous=FALSE, spacing=8)
    helpWin$add(helpVbox)
    helpHbox <- RGtk2::gtkHBoxNew(homogeneous=FALSE, spacing=8)
    helpVbox$packStart(helpHbox, expand=FALSE)
    helpText <- sprintf('%s\n\n\n%s',
                        'Instructions: Enter simple logical expressions
                     separated by a single ampersand (&) or pipe (|)
                     The vertex attributes must already exist.',
                     'Example: degree > 23 & btwn.cent > 50')
    helpLabel <- RGtk2::gtkLabelNew(helpText)
    helpVbox$packStart(helpLabel, expand=FALSE)
    helpButton <- RGtk2::gtkButtonNewFromStock('gtk-ok')
    RGtk2::gSignalConnect(helpButton, signal='clicked', function(widget) helpWin$destroy())
    helpVbox$packStart(helpButton, fill=FALSE)
}

# Get the screen size
#-----------------------------------------------------------
get_screen_size <- function() {
  if (.Platform$OS.type == 'windows') {
    scr_width <- system("wmic desktopmonitor get screenwidth", intern=TRUE)
    scr_height <- system("wmic desktopmonitor get screenheight", intern=TRUE)
    res <- as.numeric(c(scr_width[-c(1L, length(scr_width))],
                        scr_height[-c(1L, length(scr_height))]))
  } else {
    display.size <- system("xrandr | grep '*' | awk '{print $1}'", intern=TRUE)
    res <- as.numeric(strsplit(display.size, 'x')[[1L]])
  }
  res[1L] <- res[1L] / 3 - 10  # To fit the left panel + 2 plot areas
  res[2L] <- if (res[2L] < 800) 0.8 * res[2L] else 700
  res
}

# Get number of children of the container and kill all but the 1st
# Used to kill all but "menubar" (if "container=window") and to kill "slider"
#-----------------------------------------------------------
kill_others <- function(container) {
  kNumChildren <- length(container$getChildren())
  if (kNumChildren > 1L) {
    for (i in 2:kNumChildren) container[[i]]$destroy()
  }
}
