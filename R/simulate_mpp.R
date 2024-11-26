#' Simulate a realization of a location dependent marked point process
#'
#' @param sc_params vector of parameter values corresponding to (alpha_1, beta_1, gamma_1, alpha_2, beta_2, alpha_3, beta_3, gamma_3).
#' @param t_min minimum value for time.
#' @param t_max maximum value for time.
#' @param anchor_point vector of (x,y) coordinates of point to condition on.
#' @param raster_list list of raster objects.
#' @param mark_model a model object (typically from \code{train_mark_model}).
#' @param xy_bounds a vector of domain bounds (2 for x, 2 for y).
#' @param include_comp_inds `TRUE` or `FALSE` indicating whether to generate and use competition indices as covariates.
#' @param competition_radius distance for competition radius if \code{include_comp_inds} is `TRUE`.
#' @param correction type of correction to apply ("none" or "toroidal").
#' @param thinning `TRUE` or `FALSE` indicating whether to thin the realization.
#'
#' @return a list containing the marked point process realization and the data frame of the realization.
#' @export
#'
#' @examples
#' # Specify the generating parameters of the self-correcting process
#' generating_parameters <- c(2, 8, .02, 2.5, 3, 1, 2.5, .2)
#'
#' # Specify an anchor point
#' M_n <- matrix(c(10, 14), ncol = 1)
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
#' mark_model <- bundle::unbundle(readRDS(file_path))
#'
#' # Simulate a realization
#' example_mpp <- simulate_mpp(
#'   sc_params = generating_parameters,
#'   t_min = 0,
#'   t_max = 1,
#'   anchor_point = M_n,
#'   raster_list = scaled_raster_list,
#'   mark_model = mark_model,
#'   xy_bounds = c(0, 25, 0, 25),
#'   include_comp_inds = TRUE,
#'   competition_radius = 10,
#'   correction = "none",
#'   thinning = TRUE
#' )
#'
#' # Plot the realization
#' plot_mpp(example_mpp$mpp, pattern_type = "simulated")
#'
simulate_mpp <- function(sc_params = NULL,
                         t_min = 0,
                         t_max = 1,
                         anchor_point = NULL,
                         raster_list = NULL,
                         mark_model = NULL,
                         xy_bounds = NULL,
                         include_comp_inds = FALSE,
                         competition_radius = 15,
                         correction = "none",
                         thinning = TRUE) {
  # Check the arguments
  if (length(sc_params) != 8 | anyNA(sc_params) | any(sc_params[2:8] < 0)) stop("Provide a valid set of parameter values for the sc_params argument.", .call = FALSE)
  if (t_min < 0 | t_min >= t_max | is.null(t_min)) stop("Provide a value greater than 0 and less than t_max for the t_min argument.", .call = FALSE)
  if (t_max > 1 | t_min >= t_max | is.null(t_max)) stop("Provide a value greater than t_min and less than 1 for the t_max argument.", .call = FALSE)
  if (length(anchor_point) != 2) stop("Provide a vector of (x,y) coordinates for the anchor_point argument.", .call = FALSE)
  if (is.null(raster_list) | !is.list(raster_list)) stop("Provide a list of rasters for the raster_list argument.", .call = FALSE)
  if (is.null(mark_model)) stop("Provide an unbundled mark model for the mark_model argument.", .call = FALSE)
  if (is.null(xy_bounds) | !(length(xy_bounds) == 4)) stop("Provide (x,y) bounds in the form (a_x, b_x, a_y, b_y) for the xy_bounds argument.", .call = FALSE)
  if (xy_bounds[2] < xy_bounds[1] | xy_bounds[4] < xy_bounds[3]) stop("Provide (x,y) bounds in the form (a_x, b_x, a_y, b_y) for the xy_bounds argument.", .call = FALSE)
  if (!correction %in% c("none", "toroidal")) stop("Provide a valid correction type.", .call = FALSE)
  if (include_comp_inds == TRUE & (is.null(competition_radius) | competition_radius < 0)) stop("Provide the desired radius for competition indices.", .call = FALSE)
  if (!is.logical(thinning)) stop("Provide a logical value for the thinning argument.", .call = FALSE)

  # Simulate times and locations
  sim_times <- stats::na.omit(sim_temporal_sc(t_min, t_max, sc_params[1:3]))
  sim_locs <- sim_spatial_sc(anchor_point, sc_params[4:5], length(sim_times), xy_bounds)
  sim_times[1] <- 0
  txy_sim <- base::cbind(sim_times, sim_locs)

  # Perform the thinning process
  thin_vals <- (stats::runif(base::nrow(txy_sim), 0, 1) < interaction_st(txy_sim, sc_params[6:8]))
  txy_sim_thin <- txy_sim[thin_vals, ]

  # Compile the thinned and unthinned results
  sim_df <- base::data.frame(time = txy_sim[, 1], x = txy_sim[, 2], y = txy_sim[, 3])
  sim_thin_df <- base::data.frame(time = txy_sim_thin[, 1], x = txy_sim_thin[, 2], y = txy_sim_thin[, 3])

  if (thinning) {
    # Obtain a matrix of (x, y) locations
    s <- sim_thin_df[, c("x", "y")]

    # Predict the mark values
    marks <- predict_marks(
      sim_realization = sim_thin_df,
      raster_list = raster_list,
      mark_model = mark_model,
      xy_bounds = xy_bounds,
      include_comp_inds = include_comp_inds,
      competition_radius = competition_radius,
      correction = correction
    )

    # Compile the results
    sim_realization <- data.frame(time = sim_thin_df$time, x = sim_thin_df$x, y = sim_thin_df$y, marks = marks$.pred)
    sim_mpp <- generate_mpp(locations = sim_realization[, c("x", "y")], marks = sim_realization$marks, xy_bounds = xy_bounds)
  } else {
    # Obtain a matrix of (x, y) locations
    s <- sim_df[, c("x", "y")]

    # Predict the mark values
    marks <- predict_marks(
      sim_realization = sim_df,
      raster_list = raster_list,
      mark_model = mark_model,
      xy_bounds = xy_bounds,
      include_comp_inds = include_comp_inds,
      competition_radius = competition_radius,
      correction = correction
    )

    # Compile the results
    sim_realization <- data.frame(time = sim_df$time, x = sim_df$x, y = sim_df$y, marks = marks$.pred)
    sim_mpp <- generate_mpp(locations = sim_realization[, c("x", "y")], marks = sim_realization$marks, xy_bounds = xy_bounds)
  }

  return(list(
    mpp = sim_mpp,
    data = sim_realization
  ))
}
