#' Extract covariate values from a set of rasters
#'
#' @param locations a data frame of (x,y) locations with names "x" and "y".
#' @param raster_list a list of raster objects.
#'
#' @return a matrix of covariates drawn from the provided rasters.
#' @export
#'
extract_covars <- function(locations, raster_list) {
  # Extract covariate values from the raster list and collate into a data frame
  base::do.call(base::cbind, base::lapply(raster_list,  function(q) terra::extract(q, y = locations, method = "bilinear")))
}
