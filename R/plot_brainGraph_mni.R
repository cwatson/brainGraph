#' Draw an axial or sagittal slice of the MNI152 T1 image
#'
#' This function draws an axial or sagittal slice from the MNI152 T1 image, to
#' plot the vertices of a graph over it. It will optionally write to a filename
#' for output.
#'
#' @param plane Character string, either 'axial' or 'sagittal'
#' @param slice The x or z-coordinate of the slice to use
#' @param hemi Character string, either 'left' or 'right'
#' @param save Logical indicating whether or not a png file should be saved
#' (default: FALSE)
#' @param fname The name of the file to be saved
#' @export
#' @seealso \code{\link[oro.nifti]{image.nifti}}

plot_brainGraph_mni <- function(plane=c('axial', 'sagittal'), slice,
                                hemi=c('left', 'right'),
                                save=FALSE, fname=NULL) {
    if (isTRUE(save)) {
      png(filename=fname)
    } else {
      if (length(dev.list()) == 0) {
        dev.new()
      }
    }

    plane <- match.arg(plane)
    if (plane == 'axial') {
      slice <- 46
      X <- mni152
      slicemax <- max(X[, , slice])

    } else {
      hemi <- match.arg(hemi)
      if (hemi == 'right') {
        X <- mni152
      } else if (hemi == 'left') {
        tmp <- mni152@.Data
        X <- nifti(tmp[rev(seq_len(nrow(tmp))), rev(seq_len(ncol(tmp))), ])
      } else {
        stop('Argument "hemi" must be "left" or "right"!')
      }

      #slice <- 50
      slicemax <- max(X[slice, , ])
    }

    image(X, plot.type='single', plane=plane, z=slice, zlim=c(3500, slicemax))
    par(new=T, mai=c(0, 0, 0, 0), mar=c(0, 0, 0, 0))
}
