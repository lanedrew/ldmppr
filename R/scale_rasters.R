#' Scale a set of rasters
#'
#' @param raster_list a list of raster objects.
#' @param reference_resolution the resolution to resample the rasters to.
#'
#' @return a list of scaled raster objects.
#' @export
#'
#' @examples
#' # Create two example rasters
#' rast_a <- terra::rast(
#'   ncol = 10, nrow = 10,
#'   xmin = 0, xmax = 10,
#'   ymin = 0, ymax = 10,
#'   vals = runif(100)
#' )
#'
#' rast_b <- terra::rast(
#'   ncol = 10, nrow = 10,
#'   xmin = 0, xmax = 10,
#'   ymin = 0, ymax = 10,
#'   vals = runif(100)
#' )
#'
#' # Scale example rasters in a list
#' rast_list <- list(rast_a, rast_b)
#' scale_rasters(rast_list)
#'
scale_rasters <- function(raster_list, reference_resolution = NULL) {

  # Rescale each raster
  scaled_rasters <- lapply(raster_list, terra::scale)

  # Get the current extent
  current_extents <- lapply(scaled_rasters, terra::ext)

  # Compute the new extent
  new_extents <- lapply(current_extents, function(ext) {
    x_range <- ext[2] - ext[1]
    y_range <- ext[4] - ext[3]
    new_ext <- c(0, x_range, 0, y_range)
    return(new_ext)
  })

  # Align raster to the new extent
  scaled_rasters <- mapply(function(x, y) {
    terra::ext(x) <- y
    return(x)
  }, scaled_rasters, new_extents, SIMPLIFY = FALSE)

  # Resample if reference resolution is provided
  if (!is.null(reference_resolution)) {
    scaled_rasters <- lapply(scaled_rasters, function(x) terra::resample(x, reference_resolution))
  }

  return(scaled_rasters)
}
