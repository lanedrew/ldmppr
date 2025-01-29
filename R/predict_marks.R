#' Predict values from the mark distribution
#'
#' @param sim_realization a data frame containing a thinned or unthinned realization from \code{simulate_sc}.
#' @param raster_list a list of raster objects.
#' @param scaled_rasters `TRUE` or `FALSE` indicating whether the rasters have been scaled.
#' @param mark_model a model object (typically from \code{train_mark_model}).
#' @param xy_bounds  a vector of domain bounds (2 for x, 2 for y).
#' @param include_comp_inds `TRUE` or `FALSE` indicating whether to generate and use competition indices as covariates.
#' @param competition_radius distance for competition radius if \code{include_comp_inds} is `TRUE`.
#' @param correction type of correction to apply ("none" or "toroidal").
#'
#' @return a vector of predicted mark values.
#' @export
#'
#' @examples
#' # Simulate a realization
#' generating_parameters <- c(2, 8, .02, 2.5, 3, 1, 2.5, .2)
#' M_n <- matrix(c(10, 14), ncol = 1)
#' generated_locs <- simulate_sc(
#'   t_min = 0,
#'   t_max = 1,
#'   sc_params = generating_parameters,
#'   anchor_point = M_n,
#'   xy_bounds = c(0, 25, 0, 25)
#' )
#'
#' # Load the raster files
#' raster_paths <- list.files(system.file("extdata", package = "ldmppr"),
#'   pattern = "\\.tif$", full.names = TRUE
#' )
#' rasters <- lapply(raster_paths, terra::rast)
#'
#' # Scale the rasters
#' scaled_raster_list <- scale_rasters(rasters)
#'
#' # Load the example mark model
#' file_path <- system.file("extdata", "example_mark_model.rds", package = "ldmppr")
#' example_mark_model <- readRDS(file_path)
#'
#' # Unbundle the model
#' mark_model <- bundle::unbundle(example_mark_model)
#'
#' # Predict the mark values
#' predict_marks(
#'   sim_realization = generated_locs$thinned,
#'   raster_list = scaled_raster_list,
#'   scaled_rasters = TRUE,
#'   mark_model = mark_model,
#'   xy_bounds = c(0, 25, 0, 25),
#'   include_comp_inds = TRUE,
#'   competition_radius = 10,
#'   correction = "none"
#' )
#'
predict_marks <- function(sim_realization,
                          raster_list = NULL,
                          scaled_rasters = FALSE,
                          mark_model = NULL,
                          xy_bounds = NULL,
                          include_comp_inds = FALSE,
                          competition_radius = 15,
                          correction = "none") {
  # Check the arguments
  if (!is.data.frame(sim_realization)) stop("Provide a thinned or unthinned simulation realization from simulate sc for the sim_realization argument.", .call = FALSE)
  if (is.null(raster_list) | !is.list(raster_list)) stop("Provide a list of rasters for the raster_list argument.", .call = FALSE)
  if (is.null(mark_model)) stop("Provide an unbundled mark model for the mark_model argument.", .call = FALSE)
  if (is.null(xy_bounds) | !(length(xy_bounds) == 4)) stop("Provide (x,y) bounds in the form (a_x, b_x, a_y, b_y) for the xy_bounds argument.", .call = FALSE)
  if (xy_bounds[2] < xy_bounds[1] | xy_bounds[4] < xy_bounds[3]) stop("Provide (x,y) bounds in the form (a_x, b_x, a_y, b_y) for the xy_bounds argument.", .call = FALSE)
  if (!correction %in% c("none", "toroidal")) stop("Provide a valid correction type.", .call = FALSE)
  if (include_comp_inds == TRUE & (is.null(competition_radius) | competition_radius < 0)) stop("Provide the desired radius for competition indices.", .call = FALSE)
  if (!is.logical(include_comp_inds)) stop("Provide a logical value for include_comp_inds.", .call = FALSE)
  if (!is.logical(scaled_rasters)) stop("Provide a logical value for scaled_rasters.", .call = FALSE)

  # Obtain a matrix of (x, y) locations
  s <- sim_realization[, c("x", "y")]

  # Scale the rasters if not already scaled
  if (scaled_rasters == FALSE) {
    raster_list <- scale_rasters(raster_list)
  }

  # Obtain the location specific covariate values from the scaled rasters
  X <- extract_covars(locations = s, raster_list = raster_list)
  X$x <- sim_realization$x
  X$y <- sim_realization$y
  X$time <- sim_realization$time

  if (include_comp_inds == TRUE) {
    # Calculate competition indices in a 15 unit radius
    X$near_nbr_dist <- NA
    X$near_nbr_num <- NA
    X$avg_nbr_dist <- NA
    X$near_nbr_time <- NA
    X$near_nbr_time_all <- NA
    X$near_nbr_time_dist_ratio <- NA
    colnames(s) <- c("x", "y")

    # Calculate distance matrices for selected correction method
    if (correction == "none") {
      distance_matrix <- base::as.matrix(stats::dist(s, method = "euclidean"))
    } else if (correction == "toroidal") {
      distance_matrix <- toroidal_dist_matrix_optimized(s, xy_bounds[2] - xy_bounds[1], xy_bounds[4] - xy_bounds[3])
    }


    for (i in 1:base::nrow(X)) {
      close_points <- base::unique(base::which(distance_matrix[i, ] < competition_radius & distance_matrix[i, ] != 0))
      close_times <- X$time[close_points]
      X$near_nbr_dist[i] <- base::min(distance_matrix[i, ][-i])
      X$near_nbr_num[i] <- base::length(close_points)
      X$avg_nbr_dist[i] <- base::mean(distance_matrix[i, ][close_points])
      if (base::length(close_points) == 0) {
        X$avg_nbr_dist[i] <- base::min(distance_matrix[i, ][-i])
      }
      X$near_nbr_time[i] <- X$time[base::unique(base::which(distance_matrix[i, ] == X$near_nbr_dist[i]))]
      X$near_nbr_time_all[i] <- mean(close_times)
      if (base::length(close_points) == 0) {
        X$near_nbr_time_all[i] <- X$time[base::unique(base::which(distance_matrix[i, ] == X$near_nbr_dist[i]))]
      }
      X$near_nbr_time_dist_ratio[i] <- X$near_nbr_time[i] / X$near_nbr_dist[i]
    }
  }

  return(parsnip::predict.model_fit(mark_model, X))
}
