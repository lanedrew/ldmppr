#' Extract covariate values from a set of rasters
#'
#' @param locations a matrix/data.frame of (x,y) locations.
#' @param raster_list a list of SpatRaster objects.
#'
#' @return a data.frame of covariates (no ID column; unique names).
#' @export
#'
#' @examples
#' # Load example raster data
#' raster_paths <- list.files(system.file("extdata", package = "ldmppr"),
#'   pattern = "\\.tif$", full.names = TRUE
#' )
#' raster_paths <- raster_paths[!grepl("_med\\.tif$", raster_paths)]
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
  if (is.data.frame(locations)) locations <- as.matrix(locations)
  if (!is.matrix(locations) || ncol(locations) != 2) {
    stop("`locations` must be a 2-column matrix/data.frame of (x,y).", call. = FALSE)
  }
  if (is.null(raster_list) || !is.list(raster_list) || length(raster_list) == 0) {
    stop("`raster_list` must be a non-empty list of SpatRaster objects.", call. = FALSE)
  }

  # Stack into a multi-layer SpatRaster
  r_stack <- terra::rast(raster_list)

  # Ensure layer names exist & are unique
  nm <- names(r_stack)
  if (is.null(nm) || any(!nzchar(nm))) nm <- rep("cov", terra::nlyr(r_stack))
  names(r_stack) <- make.unique(nm, sep = "__")

  X <- terra::extract(r_stack, locations, method = "bilinear")

  # Drop terra's ID column
  if ("ID" %in% names(X)) X[["ID"]] <- NULL

  # Final safety: enforce unique names
  names(X) <- make.unique(names(X), sep = "__")

  X
}
