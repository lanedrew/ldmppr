#' Simulate a realization of a location dependent marked point process
#'
#' @param sc_params vector of parameter values corresponding to
#'   (alpha_1, beta_1, gamma_1, alpha_2, beta_2, alpha_3, beta_3, gamma_3).
#' @param t_min minimum value for time.
#' @param t_max maximum value for time.
#' @param anchor_point vector (or 1x2 matrix) of (x,y) coordinates to condition on.
#' @param raster_list list of raster objects.
#' @param scaled_rasters `TRUE` or `FALSE` indicating whether the rasters have been scaled.
#' @param mark_model a mark model (e.g., from \code{train_mark_model()}), or a path to a saved mark model.
#' @param xy_bounds a vector of domain bounds (a_x, b_x, a_y, b_y).
#' @param include_comp_inds `TRUE` or `FALSE` indicating whether to generate competition indices as covariates.
#' @param competition_radius distance for competition radius if \code{include_comp_inds} is `TRUE`.
#' @param correction type of correction to apply ("none" or "toroidal").
#' @param thinning `TRUE` or `FALSE` indicating whether to thin the realization.
#'
#' @return An object of class \code{"ldmppr_sim"}.
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
#' raster_paths <- raster_paths[!grepl("_med\\.tif$", raster_paths)]
#' rasters <- lapply(raster_paths, terra::rast)
#'
#' # Scale the rasters
#' scaled_raster_list <- scale_rasters(rasters)
#'
#' # Load the example mark model
#' file_path <- system.file("extdata", "example_mark_model.rds", package = "ldmppr")
#' mark_model <- load_mark_model(file_path)
#'
#' # Simulate a realization
#' example_mpp <- simulate_mpp(
#'   sc_params = generating_parameters,
#'   t_min = 0,
#'   t_max = 1,
#'   anchor_point = M_n,
#'   raster_list = scaled_raster_list,
#'   scaled_rasters = TRUE,
#'   mark_model = mark_model,
#'   xy_bounds = c(0, 25, 0, 25),
#'   include_comp_inds = TRUE,
#'   competition_radius = 10,
#'   correction = "none",
#'   thinning = TRUE
#' )
#'
#' # Plot the realization and provide a summary
#' plot(example_mpp, pattern_type = "simulated")
#' summary(example_mpp)
simulate_mpp <- function(sc_params = NULL,
                         t_min = 0,
                         t_max = 1,
                         anchor_point = NULL,
                         raster_list = NULL,
                         scaled_rasters = FALSE,
                         mark_model = NULL,
                         xy_bounds = NULL,
                         include_comp_inds = FALSE,
                         competition_radius = 15,
                         correction = "none",
                         thinning = TRUE) {
  # ---- Checks ----
  if (is.null(sc_params) || length(sc_params) != 8 || anyNA(sc_params) || any(sc_params[2:8] < 0)) {
    stop("Provide a valid set of parameter values for the sc_params argument.", call. = FALSE)
  }

  if (is.null(t_min) || t_min < 0 || t_min >= t_max) stop("Provide t_min >= 0 and < t_max.", call. = FALSE)
  if (is.null(t_max) || t_max > 1 || t_min >= t_max) stop("Provide t_max > t_min and <= 1.", call. = FALSE)

  anchor_point <- as.numeric(anchor_point)
  if (length(anchor_point) != 2) stop("Provide a vector/matrix of (x,y) coordinates for anchor_point.", call. = FALSE)

  if (is.null(raster_list) || !is.list(raster_list)) stop("Provide a list of rasters for raster_list.", call. = FALSE)
  if (!is.logical(scaled_rasters)) stop("Provide a logical value for scaled_rasters.", call. = FALSE)

  if (is.null(xy_bounds) || length(xy_bounds) != 4) stop("Provide xy_bounds = c(a_x, b_x, a_y, b_y).", call. = FALSE)
  if (xy_bounds[2] < xy_bounds[1] || xy_bounds[4] < xy_bounds[3]) stop("Invalid xy_bounds ordering.", call. = FALSE)

  if (!correction %in% c("none", "toroidal")) stop("Provide correction = 'none' or 'toroidal'.", call. = FALSE)
  if (!is.logical(thinning)) stop("Provide a logical value for thinning.", call. = FALSE)
  if (!is.logical(include_comp_inds)) stop("Provide a logical value for include_comp_inds.", call. = FALSE)
  if (isTRUE(include_comp_inds) && (is.null(competition_radius) || competition_radius < 0)) {
    stop("Provide a nonnegative competition_radius.", call. = FALSE)
  }

  # Accept new or legacy mark model formats
  mark_model <- as_mark_model(mark_model)

  # ---- Simulate times and locations ----
  sim_times <- stats::na.omit(sim_temporal_sc(t_min, t_max, sc_params[1:3]))
  sim_locs <- sim_spatial_sc(anchor_point, sc_params[4:5], length(sim_times), xy_bounds)

  # Ensure first time is exactly 0 (as in your current implementation)
  if (length(sim_times) >= 1) sim_times[1] <- 0

  txy_sim <- cbind(sim_times, sim_locs)

  # ---- Thinning ----
  thin_keep <- stats::runif(nrow(txy_sim), 0, 1) < interaction_st(txy_sim, sc_params[6:8])
  txy_sim_thin <- txy_sim[thin_keep, , drop = FALSE]

  sim_df <- data.frame(time = txy_sim[, 1], x = txy_sim[, 2], y = txy_sim[, 3])
  sim_thin_df <- data.frame(time = txy_sim_thin[, 1], x = txy_sim_thin[, 2], y = txy_sim_thin[, 3])

  used_df <- if (isTRUE(thinning)) sim_thin_df else sim_df

  if (nrow(used_df) == 0) {
    realization <- data.frame(time = numeric(0), x = numeric(0), y = numeric(0), marks = numeric(0))
    sim_mpp <- generate_mpp(locations = matrix(numeric(0), ncol = 2), marks = numeric(0), xy_bounds = xy_bounds)
    bounds <- list(t_min = t_min, t_max = t_max, xy_bounds = xy_bounds)
    return(new_ldmppr_sim(
      process = "self_correcting",
      mpp = sim_mpp,
      realization = realization,
      params = sc_params,
      bounds = bounds,
      anchor_point = anchor_point,
      thinning = thinning,
      correction = correction,
      include_comp_inds = include_comp_inds,
      competition_radius = competition_radius,
      call = match.call(),
      meta = list(unthinned = sim_df, thinned = sim_thin_df)
    ))
  }


  # ---- Predict marks ----
  marks <- predict_marks(
    sim_realization = used_df,
    raster_list = raster_list,
    scaled_rasters = scaled_rasters,
    mark_model = mark_model,
    xy_bounds = xy_bounds,
    include_comp_inds = include_comp_inds,
    competition_radius = competition_radius,
    correction = correction
  )

  # Accept either a tibble/data.frame with `.pred` OR a numeric vector
  marks_vec <- if (is.data.frame(marks)) {
    if (!(".pred" %in% names(marks))) stop("predict_marks() returned a data.frame without a `.pred` column.", call. = FALSE)
    marks[[".pred"]]
  } else {
    as.numeric(marks)
  }

  realization <- data.frame(time = used_df$time, x = used_df$x, y = used_df$y, marks = marks_vec)
  sim_mpp <- generate_mpp(
    locations = realization[, c("x", "y")],
    marks = realization$marks,
    xy_bounds = xy_bounds
  )

  bounds <- list(t_min = t_min, t_max = t_max, xy_bounds = xy_bounds)

  new_ldmppr_sim(
    process = "self_correcting",
    mpp = sim_mpp,
    realization = realization,
    params = sc_params,
    bounds = bounds,
    anchor_point = anchor_point,
    thinning = thinning,
    correction = correction,
    include_comp_inds = include_comp_inds,
    competition_radius = competition_radius,
    call = match.call(),
    meta = list(unthinned = sim_df, thinned = sim_thin_df)
  )
}
