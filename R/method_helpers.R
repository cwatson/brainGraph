# Get info from brainGraph and brainGraphList objects
#-------------------------------------------------------------------------------
print_bg_summary <- function(object) {
  name <- name_str <- dens.pct <- Group <- modality <- weighting <- thresh <- clustmethod <- 'N/A'

  ver <- vapply(object$version, as.character, character(1L))
  date_created <- sub('T', ' ', object$date)
  type <- simpleCap(object$type)
  atlasfull <-
    switch(object$atlas,
           aal116='AAL-116', aal2.120=, aal2.94='AAL2', aal90='AAL-90',
           brainsuite='Brainsuite', craddock200='Craddock-200',
           destrieux='Destrieux', destrieux.scgm='Destrieux + SCGM',
           dk='Desikan-Killiany', dk.scgm='Desikan-Killiany + SCGM',
           dkt='Desikan-Killiany-Tourville', dkt.scgm='Desikan-Killiany-Tourville + SCGM',
           dosenbach160='Dosenbach-160', hoa112='Harvard-Oxford cortical and subcortical',
           lpba40='LONI probabilistic brain atlas', hcp_mmp1.0='HCP MMP',
           brainnetome='Brainnetome', power264='Power-264', gordon333='Gordon-333', object$atlas)
  if (!is.null(object$modality)) {
    modality <-
      switch(object$modality, dti='DTI', fmri='fMRI', thickness='Cortical thickness',
             area='Cortical surface area', volume='Cortical/subcortical volume', object$modality)
  }
  if (!is.null(object$weighting)) {
    weighting <-
        switch(object$weighting, fa='FA (fractional anisotropy)',
               sld='Streamline density', pearson='Pearson correlation',
               spearman='Spearman\'s rank correlation',
               kendall='Kendall\'s rank correlation', partial='Partial correlation',
               object$weighting)
  }
  thresh_str <- 'Threshold:'
  if (!is.null(object$threshold)) {
    thresh <- prettyNum(object$threshold, ',', digits=4L)
    if (length(thresh) > 1L) {
      thresh_str <- sub(':', 's:', thresh_str)
      if (length(thresh) > 6L) thresh <- thresh[1L:6L]
      thresh <- paste(thresh, collapse='; ')
    }
  }
  if (!is.null(object$clust.method)) {
    clustmethod <-
      switch(object$clust.method, edge_betweenness='Edge betweenness',
             fast_greedy='Greedy optimization (hierarchical agglomeration)',
             infomap='Infomap', label_prop='Label propagation',
             leading_eigen='Leading eigenvector',
             louvain='Louvain (multi-level modularity optimization)', optimal='Optimal',
             spinglass='Potts spin glass model', walktrap='Walktrap algorithm',
             object$clust.method)
  }
  if (is.brainGraph(object)) {
    if (!is_weighted(object)) weighting <- 'Unweighted'
    dens.pct <- sprintf('%1.2f%s', 100 * graph.density(object), '%')
    if (!is.null(object$name)) name <- object$name
    if (!is.null(object$Group)) Group <- object$Group
    name_str <- switch(object$level,
                       group='Group:',
                       contrast='Contrast:',
                       'Subject ID:')
  } else if (is.brainGraphList(object)) {
    if (!is_weighted(object$graphs[[1L]])) weighting <- 'Unweighted'
  }

  df <- data.frame(
          A=c('Software versions',
              '       R release:', '      brainGraph:', '          igraph:',
              'Date created:', 'Observed or random?', 'Brain atlas used:',
              'Imaging modality:', 'Edge weighting:', 'Clustering method:',
              'Graph density:', thresh_str, name_str, 'Group:'),
          B=c('', ver, date_created, type, atlasfull, modality, weighting,
              clustmethod, dens.pct, thresh, name, Group))
  dimnames(df)[[2L]] <- rep.int('', 2L)
  return(df)
}

# Print a character vector "x" as a data.frame with "numCols" columns
#-------------------------------------------------------------------------------
print_text_vector <- function(x, numCols) {
  div <- length(x) %/% numCols
  splits <- split(x, ceiling(seq_along(x) / div))
  lens <- lengths(splits)
  nsplits <- length(splits)
  splits[[nsplits]] <- c(splits[[nsplits]], rep.int('', div - lens[nsplits]))
  attrs.df <- as.data.frame(splits, stringsAsFactors=FALSE)
  dimnames(attrs.df)[[2L]] <- rep.int('', nsplits)
  return(attrs.df)
}

# Print a simple title
#-------------------------------------------------------------------------------
print_title_summary <- function(...) {
  title <- paste0(..., collapse='')
  width <- max(getOption('width') / 2L, nchar(title))
  message('\n', rep.int('=', width))
  message(title)
  message(rep.int('=', width))
}

# Print info about the model
#-------------------------------------------------------------------------------
print_model_summary <- function(x) {
  cat('GLM formula:\n\t', formula(x), '\n')
  cat('based on', nobs(x), 'observations, with', df.residual(x), 'degrees of freedom.')
  cat('\n\n')
  invisible(x)
}

# Print details about the contrast matrix
#-------------------------------------------------------------------------------
print_contrast_type_summary <- function(x) {
  cat('Contrast type: ', paste(toupper(x$con.type), 'contrast'), '\n')
  symb <- switch(x$alt, two.sided='!=', greater='>', less='<')
  alt <- sprintf('C %s 0', symb)
  cat('Alternative hypothesis: ', alt, '\n')
  cat('Contrasts: ', '\n')

  con <- x$contrasts
  if (is.matrix(con)) con <- list(con)
  for (i in seq_along(con)) {
    con[[i]] <- as.character(MASS::fractions(con[[i]]))
    con[[i]] <- sub('^0$', '.', con[[i]])
    if (is.list(x$contrasts)) message(x$con.name[i])
    print(format(con[[i]], justify='right'), quote=FALSE)
    cat('\n')
  }

  invisible(x)
}

# Print the subjects removed due to incomplete data
#-------------------------------------------------------------------------------
print_subs_summary <- function(x) {
  n <- length(x$removed.subs)
  if (n > 0L) {
    message(n, ' subjects removed due to incomplete data:')
    cat('  ', paste(names(x$removed.subs), collapse=', '), '\n')
  }
  invisible(x)
}

# Print per-contrast statistics
#-------------------------------------------------------------------------------
print_contrast_stats_summary <- function(x, ...) {
  Contrast <- csize <- p.perm <- `p-value` <- NULL
  printCon <- if (is.null(x$printCon)) x$DT[, unique(Contrast)] else x$printCon

  cls <- class(x)
  if ('summary.NBS' %in% cls) {
    DT <- x$DT.sum[csize > 1]
    setnames(DT, c('csize', 'ecount'), c('# vertices', '# edges'))
    DT[, `p-value` := signif(p.perm)]
    DT[, c('alt', 'N', 'p.perm') := NULL]
  } else {
    DT <- x$DT.sum[, !c(intersect(names(x$DT.sum), 'Outcome')), with=FALSE]
    oldnames <- grep('p-value', names(x$DT.sum), value=TRUE)
    newname <- if (x$con.type == 'f') 'Pr(>F)' else 'Pr(>|t|)'
    newnames <- sub('p-value', newname, oldnames)
    setnames(DT, oldnames, newnames)
  }

  for (i in printCon) {
    message(i, ': ')
    if (dim(DT[Contrast == i])[1L] == 0L) {
      message('\tNo signficant results!\n')
    } else {
      if ('summary.NBS' %in% cls) {
        printCoefmat(DT[Contrast == i, !'Contrast'], tst.ind=2L, P.values=TRUE, has.Pvalue=TRUE, digits=x$digits, ...)
      } else {
        if (isTRUE(x$print.head)) {
          print(DT[Contrast == i, !'Contrast'], topn=5L, nrows=10L, digits=x$digits)
        } else {
          print(DT[Contrast == i, !'Contrast'], digits=x$digits)
        }
      }
      cat('\n')
    }
  }
  invisible(x)
}

# Print details regarding permutation analyses
#-------------------------------------------------------------------------------
print_permutation_summary <- function(x) {
  message('\n', 'Permutation analysis', '\n', rep.int('-', getOption('width') / 4L))
  perm.method <- switch(x$perm.method,
                        freedmanLane='Freedman-Lane',
                        terBraak='ter Braak',
                        smith='Smith')
  cat('Permutation method: ', perm.method, '\n')
  cat('Partition method: ', simpleCap(x$part.method), '\n')
  cat('# of permutations: ', prettyNum(x$N, ','), '\n')

  invisible(x)
}

# Print an abbreviated graph metric in full
#-------------------------------------------------------------------------------
print_full_metric <- function(x) {
  full <- switch(x,
      coreness=, degree=, eccentricity=, diameter=, strength=, vulnerability=simpleCap(x),
      knn='K-nearest neighbor degree', mod='Modularity', mod.wt='Modularity (weighted)',
      E.global='Global efficiency', E.global.wt='Global efficiency (weighted)',
      E.local='Local efficiency', E.local.wt='Local efficiency (weighted)',
      E.nodal='Nodal efficiency', E.nodal.wt='Nodal efficiency (weighted)',
      Cp='Clustering coefficient', Lp='Characteristic path length',
      assort=, assortativity='Degree assortativity', assort.lobe='Lobe assortativity',
      btwn.cent='Betweenness centrality', ev.cent='Eigenvector centrality',
      clo.cent='Closeness centrality', lev.cent='Leverage centrality', pagerank='PageRank',
      subg.cent='Subgraph centrality', communicability='Communicability betweenness centrality',
      asymm='Edge asymmetry', max.comp='Largest connected component', num.hubs='# of hubs',
      num.tri='# of triangles', spatial.dist='Mean Euclidean distance', transitivity='Local transitivity')
  return(full)
}
