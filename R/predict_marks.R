#' Predict values from the mark distribution
#'
#' @param sim_realization a realization of simulated values
#' @param raster_list list of rasters
#' @param size_model a predictive model
#'
#' @return a vector of predictions
#' @export
#'
predict_marks <- function(sim_realization, raster_list, size_model){

  s <- sim_realization[,c("x", "y")]
  raster_trans <- scale_rasters(raster_list)
  X <- extract_covars(x = s, raster_list = raster_trans)
  X$x <- s[,1]
  X$y <- s[,2]
  X$time <- sim_realization[,1]

  return(stats::predict(size_model, X))
}
