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
#' raster_paths <- list.files(system.file("extdata", package = "ldmppr"),
#'   pattern = "\\.tif$", full.names = TRUE
#' )
#' rasters <- lapply(raster_paths, terra::rast)
#'
#' # Scale the rasters
#' scaled_raster_list <- scale_rasters(rasters)
#'
#' # Load example locations
#' locations <- small_example_data %>%
#'   dplyr::select(x, y) %>%
#'   as.matrix()
#'
#' # Extract covariates
#' example_covars <- extract_covars(locations, scaled_raster_list)
#' head(example_covars)
#'
extract_covars <- function(locations, raster_list) {
  # Extract covariate values from the raster list and collate into a data frame
  base::do.call(base::cbind, base::lapply(raster_list, function(q) terra::extract(q, y = locations, method = "bilinear")))
}
