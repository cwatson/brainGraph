#' Import data for structural connectivity analysis
#'
#' Given a directory, atlas name, and modality, this function imports data for
#' structural connectivity analysis. It expects files containing a table of
#' region-wise structural MRI measures (e.g., mean cortical thickness), one for
#' each hemisphere. The first column of all files should contain the
#' \emph{subject ID}; the column name will be changed to \code{Study.ID}.
#'
#' The files should have specific names; the second in the following list is
#' only required for atlases/parcellations that include \emph{subcortical gray
#' matter} (e.g., \code{dk.scgm}).
#' \itemize{
#'   \item \code{${parcellation}_${hemi}_${modality}.csv} for cortical volume,
#'     thickness, surface area, or local gyrification index (LGI). Here,
#'     \code{${parcellation}} can be \code{aparc}, \code{aparc.DKTatlas40},
#'     or \code{aparc.a2009s}. For example, for cortical thickness with the
#'     \emph{Desikan-Killiany} atlas, the filename should be
#'     \code{aparc_lh_thickness.csv}. If you are using a custom atlas, see the
#'     \emph{Note} below. The \code{${hemi}} variable is either \code{lh} or
#'     \code{rh}. Finally, \code{${modality}} should be either \code{volume},
#'     \code{thickness}, \code{area}, or \code{lgi}.
#'   \item \code{asegstats.csv} for SCGM volume
#' }
#'
#' @note When using a custom atlas, the name of the atlas's data.table should
#'   match the \code{${parcellation}} portion of the filename (specification
#'   shown above). Furthermore, it must conform to the output of Freesurfer's
#'   \code{aparcstats2table} (and \code{asegstats2table}, if applicable).
#'   Otherwise, please contact me for inclusion of a different data type.
#'
#' @param datadir The path to the directory containing the data files
#' @param atlas Character string specifying the atlas in use. For a custom
#'   atlas, please specify \code{'custom'}, and provide the name to the
#'   \code{custom.atlas} argument
#' @param modality The structural imaging measure (default: \code{'thickness'})
#' @param exclude.subs Vector indicating the subjects to exclude, if any
#'   (default: \code{NULL})
#' @param custom.atlas Character string specifying the name of the R object for
#'   the atlas in use, if \code{atlas='custom'} was also supplied (default:
#'   \code{NULL})
#' @export
#'
#' @return A list containing:
#'   \item{atlas}{Character string}
#'   \item{modality}{Character string}
#'   \item{lhrh}{A \code{data.table} of structural MRI measures for both
#'     hemispheres}
#'   \item{aseg}{A \code{data.table} of structural MRI measures for subcortical
#'     gray matter, if applicable}
#'   \item{excluded}{Vector of subject ID's that were excluded}
#'   \item{missing}{Vector of subject ID's that are not present in \emph{both}
#'     the cortical and subcortical tables (if applicable)}
#'
#' @family Structural covariance network functions
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#'   raw_data <- import_scn('/home/cwatson/data', atlas='dkt',
#'                          exclude.subs=c('con07', 'con23', 'pat15'))
#' }

import_scn <- function(datadir, atlas, modality='thickness', exclude.subs=NULL, custom.atlas=NULL) {
  missing <- scgm <- NULL

  atlas <- match.arg(atlas, choices=c(data(package='brainGraph')$results[, 3], 'custom'))
  if (atlas == 'custom') {
    stopifnot(!is.null(custom.atlas), exists(custom.atlas))
    atlas <- parc <- custom.atlas
  } else {
    parc <- switch(atlas,
                   dk=,dk.scgm='aparc',
                   dkt=,dkt.scgm='aparc.DKTatlas40',
                   destrieux=,destrieux.scgm='aparc.a2009s',
                   stop(paste0('Invalid atlas ', atlas, '. Choose from the following:\n',
                               paste(data(package='brainGraph')$results[, 3], collapse='\n'))))
  }

  # Cortical measures, both hemispheres
  files <- file.path(datadir, paste0(parc, c('_lh_', '_rh_'), modality, '.csv'))
  lh <- update_fs_names(files[1], modality, parc, 'lh', exclude.subs)
  rh <- update_fs_names(files[2], modality, parc, 'rh', exclude.subs)
  lhrh <- merge(lh, rh)

  # Subcortical measures
  if (grepl('scgm', atlas)) {
    asegfile <- file.path(datadir, 'asegstats.csv')
    scgm <- update_fs_names(asegfile, 'aseg', exclude.subs)
    sID <- getOption('bg.subject_id')
    missing <- union(setdiff(scgm[[sID]], lhrh[[sID]]), setdiff(lhrh[[sID]], scgm[[sID]]))
  }

  return(list(atlas=atlas, modality=modality, lhrh=lhrh, aseg=scgm,
              excluded=exclude.subs, missing=missing))
}

#' Update column names in a Freesurfer table
#'
#' Updates column names by abbreviating to match with the names in the
#' \code{brainGraph} atlases.
#' @keywords internal

update_fs_names <- function(filename, modality, parcellation, hemi, exclude.subs) {
  sID <- getOption('bg.subject_id')
  stopifnot(file.exists(filename))
  DT <- fread(filename, colClasses=list(character=1))
  names(DT)[1] <- sID

  if (modality == 'aseg') {
    DT <- DT[, 1:15, with=FALSE]
    names(DT) <- gsub('Left-', 'l', names(DT))
    names(DT) <- gsub('Right-', 'r', names(DT))
    names(DT) <- gsub('Thalamus.Proper', 'THAL', names(DT))
    names(DT) <- gsub('Putamen', 'PUT', names(DT))
    names(DT) <- gsub('Pallidum', 'PALL', names(DT))
    names(DT) <- gsub('Caudate', 'CAUD', names(DT))
    names(DT) <- gsub('Hippocampus', 'HIPP', names(DT))
    names(DT) <- gsub('Amygdala', 'AMYG', names(DT))
    names(DT) <- gsub('Accumbens.area', 'ACCU', names(DT))
  } else {
    # Remove the "mean" column, if present
    if (any(grepl('mean', ignore.case=TRUE, names(DT)))) {
      DT <- DT[, -grep('mean', ignore.case=TRUE, names(DT)), with=FALSE]
    }
    names(DT) <- gsub(paste0(hemi, '[_|\\.]'), substr(hemi, 1, 1), names(DT))
    names(DT) <- gsub(paste0('_', modality), '', names(DT))

    if (parcellation == 'aparc.a2009s') {
      names(DT) <- gsub('-', '.', names(DT))
    } else {
      full <- c('bankssts', 'caudalanteriorcingulate', 'caudalmiddlefrontal', 'precuneus', 'cuneus',
                'entorhinal', 'fusiform', 'inferiorparietal', 'inferiortemporal', 'isthmuscingulate',
                'lateraloccipital', 'lateralorbitofrontal', 'lingual', 'medialorbitofrontal',
                'middletemporal', 'parahippocampal', 'paracentral', 'parsopercularis',
                'parsorbitalis', 'parstriangularis', 'pericalcarine', 'postcentral',
                'posteriorcingulate', 'precentral', 'rostralanteriorcingulate',
                'rostralmiddlefrontal', 'superiorfrontal', 'superiorparietal', 'superiortemporal',
                'supramarginal', 'frontalpole', 'temporalpole', 'transversetemporal', 'insula')
      abbrev <- c('BSTS', 'cACC', 'cMFG', 'PCUN', 'CUN', 'ENT', 'FUS', 'IPL', 'ITG', 'iCC', 'LOG',
                  'LOF', 'LING', 'MOF', 'MTG', 'PARH', 'paraC', 'pOPER', 'pORB', 'pTRI', 'periCAL',
                  'postC', 'PCC', 'preC', 'rACC', 'rMFG', 'SFG', 'SPL', 'STG', 'SMAR', 'FP', 'TP',
                  'TT', 'INS')
      for (i in seq_along(full)) names(DT) <- gsub(full[i], abbrev[i], names(DT))
    }
  }

  setkeyv(DT, sID)
  if (!is.null(exclude.subs)) DT <- DT[!get(sID) %in% exclude.subs]
  return(DT)
}
