#' Simulate a realization of a location dependent marked point process
#'
#' @param process type of process used (currently supports \code{"self_correcting"}).
#' @param process_fit either (1) a \code{ldmppr_fit} object returned by
#'   \code{\link{estimate_process_parameters}}, or (2) a numeric vector of
#'   length 8 giving self-correcting parameters:
#'   \eqn{(\alpha_1,\beta_1,\gamma_1,\alpha_2,\beta_2,\alpha_3,\beta_3,\gamma_3)} (alpha_1, beta_1, gamma_1, alpha_2, beta_2, alpha_3, beta_3, gamma_3).
#' @param t_min minimum value for time.
#' @param t_max maximum value for time.
#' @param anchor_point (optional) vector of (x,y) coordinates of the point to condition on.
#'   If \code{NULL}, inferred from the reference data (largest mark if available) or from
#'   \code{process_fit$data_original} (largest size).
#' @param raster_list a list of raster objects used for predicting marks.
#' @param scaled_rasters \code{TRUE} or \code{FALSE} indicating whether the rasters have already been scaled.
#' @param mark_model a mark model object. May be a \code{ldmppr_mark_model} or a legacy model.
#' @param xy_bounds (optional) vector of bounds as \code{c(a_x, b_x, a_y, b_y)}. If \code{NULL}, will be
#'   inferred from \code{reference_data}'s window when \code{reference_data} is provided,
#'   otherwise from \code{ldmppr_fit} with lower bounds assumed to be 0.
#' @param include_comp_inds \code{TRUE} or \code{FALSE} indicating whether to compute competition indices.
#' @param competition_radius distance for competition radius if \code{include_comp_inds = TRUE}.
#' @param thinning \code{TRUE} or \code{FALSE} indicating whether to use the thinned simulated values.
#' @param edge_correction type of edge correction to apply (\code{"none"} or \code{"toroidal"}).
#' @param seed integer seed for reproducibility.
#'
#' @return an object of class \code{"ldmppr_sim"}.
#'
#' @examples
#' # Specify the generating parameters of the self-correcting process
#' generating_parameters <- c(2, 8, .02, 2.5, 3, 1, 2.5, .2)
#'
#' # Specify an anchor point
#' M_n <- c(10, 14)
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
#'   process = "self_correcting",
#'   process_fit = generating_parameters,
#'   t_min = 0,
#'   t_max = 1,
#'   anchor_point = M_n,
#'   raster_list = scaled_raster_list,
#'   scaled_rasters = TRUE,
#'   mark_model = mark_model,
#'   xy_bounds = c(0, 25, 0, 25),
#'   include_comp_inds = TRUE,
#'   competition_radius = 10,
#'   edge_correction = "none",
#'   thinning = TRUE
#' )
#'
#' # Plot the realization and provide a summary
#' plot(example_mpp, pattern_type = "simulated")
#' summary(example_mpp)
#' @export
simulate_mpp <- function(process = c("self_correcting"),
                         process_fit = NULL,
                         t_min = 0,
                         t_max = 1,
                         anchor_point = NULL,
                         raster_list = NULL,
                         scaled_rasters = FALSE,
                         mark_model = NULL,
                         xy_bounds = NULL,
                         include_comp_inds = FALSE,
                         competition_radius = 15,
                         edge_correction = "none",
                         thinning = TRUE,
                         seed = NULL) {

  process <- match.arg(process)

  # ---- checks ----
  if (is.null(process_fit)) {
    stop("Provide `process_fit` as an ldmppr_fit object or a numeric parameter vector.", call. = FALSE)
  }

  sc_params <- .as_sc_params(process_fit)
  if (length(sc_params) != 8 || anyNA(sc_params) || any(sc_params[2:8] < 0)) {
    stop("Provide a valid set of self-correcting parameters (length 8; entries 2:8 >= 0), or an ldmppr_fit object.", call. = FALSE)
  }

  if (is.null(t_min) || t_min < 0 || t_min >= t_max) stop("Provide t_min >= 0 and < t_max.", call. = FALSE)
  if (is.null(t_max) || t_max > 1 || t_min >= t_max) stop("Provide t_max > t_min and <= 1.", call. = FALSE)

  if (!is.null(seed)) {
    if (is.na(seed) || seed < 0 || seed != as.integer(seed)) stop("Provide a nonnegative integer `seed`.", call. = FALSE)
    set.seed(as.integer(seed))
  }

  if (is.null(xy_bounds)) {
    xy_bounds <- .infer_xy_bounds(process_fit)
  }
  if (is.null(xy_bounds) || length(xy_bounds) != 4) {
    stop("Provide xy_bounds = c(a_x, b_x, a_y, b_y), or pass an ldmppr_fit with grid$upper_bounds.", call. = FALSE)
  }
  if (xy_bounds[2] < xy_bounds[1] || xy_bounds[4] < xy_bounds[3]) stop("Invalid xy_bounds ordering.", call. = FALSE)

  if (is.null(anchor_point)) {
    anchor_point <- .infer_anchor_point(process_fit)
  }
  anchor_point <- as.numeric(anchor_point)
  if (length(anchor_point) != 2) {
    stop("Provide anchor_point = c(x,y), or pass an ldmppr_fit that contains data_original/data to infer it.", call. = FALSE)
  }

  if (is.null(raster_list) || !is.list(raster_list)) stop("Provide a list of rasters for raster_list.", call. = FALSE)
  if (!is.logical(scaled_rasters)) stop("Provide a logical value for scaled_rasters.", call. = FALSE)

  if (!edge_correction %in% c("none", "toroidal")) stop("Provide correction = 'none' or 'toroidal'.", call. = FALSE)
  if (!is.logical(thinning)) stop("Provide a logical value for thinning.", call. = FALSE)
  if (!is.logical(include_comp_inds)) stop("Provide a logical value for include_comp_inds.", call. = FALSE)
  if (isTRUE(include_comp_inds) && (is.null(competition_radius) || competition_radius < 0)) {
    stop("Provide a nonnegative competition_radius.", call. = FALSE)
  }

  # Accept new or legacy mark model formats
  mark_model <- as_mark_model(mark_model)

  # Scale the rasters if necessary
  if (!isTRUE(scaled_rasters)) {
    raster_list <- scale_rasters(raster_list)
    scaled_rasters <- TRUE
  }

  # ---- Simulate times and locations ----
  # NOTE: these are your existing internal helpers
  sim_times <- stats::na.omit(sim_temporal_sc(t_min, t_max, sc_params[1:3]))
  sim_locs  <- sim_spatial_sc(anchor_point, sc_params[4:5], length(sim_times), xy_bounds)

  if (length(sim_times) >= 1) sim_times[1] <- 0
  txy_sim <- cbind(sim_times, sim_locs)

  # ---- Thinning ----
  thin_keep <- stats::runif(nrow(txy_sim), 0, 1) < interaction_st(txy_sim, sc_params[6:8])
  txy_sim_thin <- txy_sim[thin_keep, , drop = FALSE]

  sim_df <- data.frame(time = txy_sim[, 1], x = txy_sim[, 2], y = txy_sim[, 3])
  sim_thin_df <- data.frame(time = txy_sim_thin[, 1], x = txy_sim_thin[, 2], y = txy_sim_thin[, 3])

  used_df <- if (isTRUE(thinning)) sim_thin_df else sim_df

  bounds <- list(t_min = t_min, t_max = t_max, xy_bounds = xy_bounds)

  if (nrow(used_df) == 0) {
    realization <- data.frame(time = numeric(0), x = numeric(0), y = numeric(0), marks = numeric(0))
    sim_mpp <- generate_mpp(locations = matrix(numeric(0), ncol = 2), marks = numeric(0), xy_bounds = xy_bounds)

    return(new_ldmppr_sim(
      process = "self_correcting",
      mpp = sim_mpp,
      realization = realization,
      params = sc_params,
      bounds = bounds,
      anchor_point = anchor_point,
      thinning = thinning,
      edge_correction = edge_correction,
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
    scaled_rasters = TRUE,
    mark_model = mark_model,
    xy_bounds = xy_bounds,
    include_comp_inds = include_comp_inds,
    competition_radius = competition_radius,
    edge_correction = edge_correction
  )

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

  new_ldmppr_sim(
    process = "self_correcting",
    mpp = sim_mpp,
    realization = realization,
    params = sc_params,
    bounds = bounds,
    anchor_point = anchor_point,
    thinning = thinning,
    edge_correction = edge_correction,
    include_comp_inds = include_comp_inds,
    competition_radius = competition_radius,
    call = match.call(),
    meta = list(unthinned = sim_df, thinned = sim_thin_df)
  )
}
