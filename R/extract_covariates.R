#' Extract covariate values from a set of rasters
#'
#' @param locations a data frame of (x,y) locations with names "x" and "y".
#' @param raster_list a list of raster objects.
#'
#' @return a matrix of covariates drawn from the provided rasters.
#' @export
#'
#' @examples
#' # Load example raster data
#' file_path <- system.file("extdata", "Snodgrass_aspect_southness_1m.tif", package = "ldmppr")
#' south <- terra::rast(file_path)
#' file_path <- system.file("extdata", "Snodgrass_wetness_index_1m.tif", package = "ldmppr")
#' wet <- terra::rast(file_path)
#'
#' # Scale the rasters
#' raster_list <- list(south, wet)
#' scaled_raster_list <- scale_rasters(raster_list)
#'
#' # Load example locations
#' locations <- small_example_data %>%
#'   dplyr::select(x, y) %>%
#'   as.matrix()
#'
#' # Extract covariates
#' extract_covars(locations, scaled_raster_list)
#'
extract_covars <- function(locations, raster_list) {
  # Extract covariate values from the raster list and collate into a data frame
  base::do.call(base::cbind, base::lapply(raster_list, function(q) terra::extract(q, y = locations, method = "bilinear")))
}
