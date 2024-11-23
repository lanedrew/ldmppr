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
  scaled_rasters <- lapply(raster_list, function(r) {
    # Scale the raster values
    r_scaled <- terra::scale(r)

    # Get the current extent
    current_extent <- terra::ext(r_scaled)

    # Compute the new extent
    x_range <- current_extent[2] - current_extent[1]
    y_range <- current_extent[4] - current_extent[3]
    new_extent <- c(0, x_range, 0, y_range)

    # Align raster to the new extent
    terra::ext(r_scaled) <- new_extent

    # Resample if reference resolution is provided
    if (!is.null(reference_resolution)) {
      r_scaled <- terra::resample(r_scaled, reference_resolution)
    }

    return(r_scaled)
  })

  return(scaled_rasters)
}
