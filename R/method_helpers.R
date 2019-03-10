# Print a simple title
print_title_summary <- function(title) {
  title <- paste(title)
  message('\n', title, '\n', rep('-', getOption('width') / 2))
}

# Print info about the outcome and/or graph measure of interest
print_measure_summary <- function(x) {
  if (x$outcome == x$measure) {
    cat('Graph metric of interest (outcome): ', x$outcome, '\n\n')
  } else {
    cat('Graph metric of interest (covariate): ', x$measure, '\n')
    cat('Outcome measure: ', x$outcome, '\n\n')
  }
  invisible(x)
}

# Print details about the contrast matrix
print_contrast_type_summary <- function(x) {
  cat('Contrast type: ', paste(toupper(x$con.type), 'contrast'), '\n')
  alt <- switch(x$alt,
                two.sided='C != 0',
                greater='C > 0',
                less='C < 0')
  cat('Alternative hypothesis: ', alt, '\n')
  cat('Contrast matrix: ', '\n')
  print(x$con.mat)

  cat('\n')
  invisible(x)
}

# Print the subjects removed due to incomplete data
print_subs_summary <- function(x) {
  n <- length(x$removed.subs)
  if (n > 0) {
    cat(n, 'subjects removed due to incomplete data:\n ', paste(x$removed.subs, collapse=', '), '\n')
  }
  invisible(x)
}

# Print per-contrast statistics
print_contrast_stats_summary <- function(x) {
  if (is.null(x$contrast)) {
    contrast <- x$DT[, unique(contrast)]
  } else {
    contrast <- x$contrast
  }
  for (i in contrast) {
    message('Contrast ', i, ': ', x$con.name[i])
    if (nrow(x$DT.sum[contrast == i]) == 0) {
      message('\tNo signficant results!\n')
    } else {
      if (isTRUE(x$print.head)) {
        print(x$DT.sum[contrast == i, !c('Contrast', 'contrast')], topn=5, nrows=10, digits=x$digits)
      } else {
        print(x$DT.sum[contrast == i, !c('Contrast', 'contrast')], digits=x$digits)
      }
      cat('\n')
    }
  }
  invisible(x)
}

# Print details regarding permutation analyses
print_permutation_summary <- function(x) {
  message('\n', 'Permutation analysis', '\n', rep('-', getOption('width') / 4))
  perm.method <- switch(x$perm.method,
                        freedmanLane='Freedman-Lane',
                        terBraak='ter Braak',
                        smith='Smith')
  part.method <- switch(x$part.method,
                        beckmann='Beckmann',
                        guttman='Guttman',
                        ridgway='Ridgway')
  cat('Permutation method: ', perm.method, '\n')
  cat('Partition method: ', part.method, '\n')
  cat('# of permutations: ', prettyNum(x$N, ','), '\n')

  invisible(x)
}
