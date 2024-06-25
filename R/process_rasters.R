#' Scale a set of rasters
#'
#' @param raster_list the list of raster objects
#'
#' @return a list of scaled raster objects
#' @export
#'
#' @examples
#' # Create two example rasters
#' rast_a <- terra::rast(ncol = 10, nrow = 10,
#'  xmin = 0, xmax = 10,
#'  ymin = 0, ymax = 10,
#'   vals = runif(100))
#'
#' rast_b <- terra::rast(ncol = 10, nrow = 10,
#'  xmin = 0, xmax = 10,
#'   ymin = 0, ymax = 10,
#'    vals = runif(100))
#'
#' # Scale example rasters in a list
#' rast_list <- list(rast_a, rast_b)
#' scale_rasters(rast_list)
scale_rasters <- function(raster_list) {
  scaled_rasters <- base::lapply(raster_list, terra::scale)
  raster_extents <- base::lapply(scaled_rasters, terra::ext)
  new_raster_extents <- base::lapply(raster_extents, function(x) c(0, x[2] - x[1], 0, x[4] - x[3]))

  for(i in 1:base::length(scaled_rasters)){
    terra::ext(scaled_rasters[[i]]) <- new_raster_extents[[i]]
  }

  return(scaled_rasters)
}


#' Extract covariate values from a set of rasters
#'
#' @param x a 2-column matrix or data.frame(x, y) or (lon, lat)
#' @param raster_list the list of raster objects
#'
#' @return a matrix of covariates drawn from the provided rasters
#' @export
#'
extract_covars <- function(x, raster_list) {
  base::do.call(base::cbind, base::lapply(raster_list,  function(q) terra::extract(q, y = x, method = "bilinear")))
}

