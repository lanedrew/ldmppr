#' Scale a set of rasters
#'
#' @param raster_list a list of raster objects.
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
scale_rasters <- function(raster_list) {
  # Scale the rasters
  scaled_rasters <- base::lapply(raster_list, terra::scale)

  # Obtain the raster extents
  raster_extents <- base::lapply(scaled_rasters, terra::ext)

  # Replace the old extents with adjust extents
  new_raster_extents <- base::lapply(raster_extents, function(x) c(0, x[2] - x[1], 0, x[4] - x[3]))
  for (i in 1:base::length(scaled_rasters)) {
    terra::ext(scaled_rasters[[i]]) <- new_raster_extents[[i]]
  }

  return(scaled_rasters)
}
