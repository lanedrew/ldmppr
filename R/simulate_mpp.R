#' Simulate a realization of a location dependent marked point process
#'
#' @param process type of process used (currently supports \code{"self_correcting"}).
#' @param process_fit either (1) a \code{ldmppr_fit} object returned by
#'   \code{\link{estimate_process_parameters}}, or (2) a numeric vector of
#'   length 8 giving the self-correcting process parameters:
#'   \eqn{(\alpha_1,\beta_1,\gamma_1,\alpha_2,\beta_2,\alpha_3,\beta_3,\gamma_3)} (alpha_1, beta_1, gamma_1, alpha_2, beta_2, alpha_3, beta_3, gamma_3).
#' @param t_min minimum value for time.
#' @param t_max maximum value for time.
#' @param anchor_point (optional) vector of (x,y) coordinates of the point to condition on.
#'   If \code{NULL}, inferred from the reference data (largest mark if available) or from
#'   \code{process_fit$data_original} (largest size).
#' @param raster_list (optional) list of raster objects used for mark prediction.
#'   Required when \code{mark_mode='mark_model'} unless rasters are stored in \code{mark_model}.
#' @param scaled_rasters \code{TRUE} or \code{FALSE} indicating whether rasters are already scaled.
#'   Ignored when \code{mark_mode='time_to_size'}.
#' @param mark_model a mark model object used when \code{mark_mode='mark_model'}.
#'   May be an \code{ldmppr_mark_model}, \code{model_fit}, or \code{workflow}.
#' @param xy_bounds (optional) vector of bounds as \code{c(a_x, b_x, a_y, b_y)}. If \code{NULL},
#'   bounds are inferred from \code{process_fit} when available.
#' @param include_comp_inds \code{TRUE} or \code{FALSE} indicating whether to compute competition indices.
#' @param competition_radius positive numeric distance used when \code{include_comp_inds = TRUE}.
#' @param thinning \code{TRUE} or \code{FALSE} indicating whether to use the thinned simulated values.
#' @param edge_correction type of edge correction to apply (\code{"none"} or \code{"toroidal"}).
#' @param seed integer seed for reproducibility.
#' @param mark_mode (optional) mark generation mode: \code{"mark_model"} uses
#'   \code{predict()} on a mark model, while \code{"time_to_size"} maps simulated
#'   times back to sizes via \code{delta}. If \code{NULL}, inferred as
#'   \code{"mark_model"} when \code{mark_model} is provided, otherwise
#'   \code{"time_to_size"}.
#' @param size_range numeric vector \code{c(smin, smax)} used for \code{mark_mode='time_to_size'}.
#'   If \code{NULL}, inferred from \code{process_fit} when possible.
#' @param delta positive scalar used for \code{mark_mode='time_to_size'}.
#'   If \code{NULL}, inferred from \code{process_fit} when possible.
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
#'   thinning = TRUE,
#'   seed = 90210
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
                         seed = NULL,
                         mark_mode = NULL,
                         size_range = NULL,
                         delta = NULL) {

  process <- match.arg(process)
  if (is.null(mark_mode)) {
    mark_mode <- if (!is.null(mark_model)) "mark_model" else "time_to_size"
  } else {
    mark_mode <- match.arg(mark_mode, c("mark_model", "time_to_size"))
  }

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

  if (!edge_correction %in% c("none", "toroidal")) stop("Provide correction = 'none' or 'toroidal'.", call. = FALSE)
  if (!is.logical(thinning)) stop("Provide a logical value for thinning.", call. = FALSE)
  if (!is.logical(include_comp_inds)) stop("Provide a logical value for include_comp_inds.", call. = FALSE)
  if (isTRUE(include_comp_inds) && (is.null(competition_radius) || competition_radius < 0)) {
    stop("Provide a nonnegative competition_radius.", call. = FALSE)
  }

  # ---- helpers for mark_mode == "time_to_size" ----
  .infer_size_range <- function(process_fit) {
    # Try common places: process_fit$data_original, process_fit$data, or marks in embedded ppp
    # We look for columns: size, marks, dbh (fallback).
    candidates <- NULL

    if (inherits(process_fit, "ldmppr_fit")) {
      if (is.data.frame(process_fit$data_original)) candidates <- process_fit$data_original
      if (is.null(candidates) && is.matrix(process_fit$data)) {
        candidates <- as.data.frame(process_fit$data)
      }
    }

    if (is.data.frame(candidates)) {
      nm <- names(candidates)
      if ("size" %in% nm) {
        v <- candidates[["size"]]
      } else if ("marks" %in% nm) {
        v <- candidates[["marks"]]
      } else if ("dbh" %in% nm) {
        v <- candidates[["dbh"]]
      } else {
        v <- NULL
      }
      if (!is.null(v)) {
        v <- as.numeric(v)
        v <- v[is.finite(v)]
        if (length(v)) return(range(v))
      }
    }

    # last resort: if process_fit contains a ppp-like object in data_original, try marks()
    if (inherits(process_fit, "ldmppr_fit") && !is.null(process_fit$ppp) &&
        spatstat.geom::is.ppp(process_fit$ppp) && !is.null(spatstat.geom::marks(process_fit$ppp))) {
      v <- as.numeric(spatstat.geom::marks(process_fit$ppp))
      v <- v[is.finite(v)]
      if (length(v)) return(range(v))
    }

    NULL
  }

  .infer_delta <- function(process_fit) {
    if (inherits(process_fit, "ldmppr_fit")) {
      # common place in your code: fit_obj$mapping$delta
      if (!is.null(process_fit$mapping) && !is.null(process_fit$mapping$delta)) {
        d <- as.numeric(process_fit$mapping$delta)
        if (length(d) == 1L && is.finite(d)) return(d)
      }
      # sometimes stored as attribute on data matrix (from .build_sc_matrix)
      if (!is.null(process_fit$data) && !is.null(attr(process_fit$data, "ldmppr_delta"))) {
        d <- as.numeric(attr(process_fit$data, "ldmppr_delta"))
        if (length(d) == 1L && is.finite(d)) return(d)
      }
    }
    NULL
  }

  .time_to_size <- function(t, smin, smax, delta) {
    t <- as.numeric(t)
    # clamp to [0,1] to be safe
    t <- pmin(pmax(t, 0), 1)
    if (!is.finite(delta) || delta <= 0) stop("delta must be a positive finite scalar for time_to_size.", call. = FALSE)
    # inverse of: t = 1 - ((size - smin)/(smax - smin))^delta
    u <- (1 - t)^(1 / delta)
    smin + u * (smax - smin)
  }

  # ---- mark_model mode requires mark_model + rasters; time_to_size does not ----
  if (mark_mode == "mark_model") {
    mark_model <- as_mark_model(mark_model)

    if (!is.null(raster_list) && !is.list(raster_list)) {
      stop("Provide `raster_list` as a list of rasters (or leave NULL to infer from `mark_model`).", call. = FALSE)
    }
    if (!is.logical(scaled_rasters) || length(scaled_rasters) != 1L) {
      stop("Provide a single logical value for `scaled_rasters`.", call. = FALSE)
    }

    # ---- infer rasters / scaled flag from mark_model if not provided ----
    if (is.null(raster_list)) {
      mm_rasters <- NULL
      if (is.list(mark_model) && !is.null(mark_model$rasters)) mm_rasters <- mark_model$rasters

      if (!is.null(mm_rasters)) {
        raster_list <- mm_rasters
        mm_scaled <- NULL
        if (is.list(mark_model) && !is.null(mark_model$info) && !is.null(mark_model$info$scaled_rasters)) {
          mm_scaled <- mark_model$info$scaled_rasters
        }
        if (isTRUE(mm_scaled)) scaled_rasters <- TRUE
      }
    }

    if (is.null(raster_list) || !is.list(raster_list)) {
      stop(
        "mark_mode='mark_model' requires `raster_list` (or an `ldmppr_mark_model` that contains `rasters`).",
        call. = FALSE
      )
    }

    if (!isTRUE(scaled_rasters)) {
      raster_list <- scale_rasters(raster_list)
      scaled_rasters <- TRUE
    }
  } else {
    # time_to_size mode: infer size_range and delta if not given
    if (is.null(size_range)) {
      size_range <- .infer_size_range(process_fit)
    }
    if (is.null(size_range) || length(size_range) != 2L || any(!is.finite(size_range))) {
      stop(
        "mark_mode='time_to_size' requires `size_range = c(smin, smax)` or an ldmppr_fit containing original sizes/marks.",
        call. = FALSE
      )
    }
    if (is.null(delta)) {
      delta <- .infer_delta(process_fit)
    }
    if (is.null(delta) || length(delta) != 1L || !is.finite(delta) || delta <= 0) {
      stop(
        "mark_mode='time_to_size' requires a positive scalar `delta`, or an ldmppr_fit with mapping$delta / ldmppr_delta.",
        call. = FALSE
      )
    }
    size_range <- as.numeric(size_range)
    delta <- as.numeric(delta)
  }

  # ---- Simulate (delegate to simulate_sc so thinning uses interaction_st_fast, and empty cases are safe) ----
  sim_obj <- simulate_sc(
    t_min = t_min,
    t_max = t_max,
    sc_params = sc_params,
    anchor_point = anchor_point,
    xy_bounds = xy_bounds
  )

  sim_df      <- sim_obj$unthinned
  sim_thin_df <- sim_obj$thinned

  used_df <- if (isTRUE(thinning)) sim_thin_df else sim_df

  bounds <- list(t_min = t_min, t_max = t_max, xy_bounds = xy_bounds)

  # ---- If empty, return empty sim object (still valid ldmppr_sim) ----
  if (is.null(used_df) || nrow(used_df) == 0) {
    realization <- data.frame(time = numeric(0), x = numeric(0), y = numeric(0), marks = numeric(0))

    sim_mpp <- generate_mpp(
      locations = matrix(numeric(0), ncol = 2),
      marks = numeric(0),
      xy_bounds = xy_bounds
    )

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
      meta = list(unthinned = sim_df, thinned = sim_thin_df, mark_mode = mark_mode)
    ))
  }

  # ---- Marks ----
  if (mark_mode == "mark_model") {

    marks <- if (inherits(mark_model, "ldmppr_mark_model")) {
      stats::predict(
        mark_model,
        sim_realization = used_df,
        raster_list = raster_list,
        scaled_rasters = TRUE,
        xy_bounds = xy_bounds,
        include_comp_inds = include_comp_inds,
        competition_radius = competition_radius,
        edge_correction = edge_correction,
        seed = seed
      )
    } else {
      predict_marks(
        sim_realization = used_df,
        raster_list = raster_list,
        scaled_rasters = TRUE,
        mark_model = mark_model,
        xy_bounds = xy_bounds,
        include_comp_inds = include_comp_inds,
        competition_radius = competition_radius,
        edge_correction = edge_correction,
        seed = seed
      )
    }

    marks_vec <- if (is.data.frame(marks)) {
      if (!(".pred" %in% names(marks))) stop("predict_marks() returned a data.frame without a `.pred` column.", call. = FALSE)
      as.numeric(marks[[".pred"]])
    } else {
      as.numeric(marks)
    }

  } else {
    # location-independent marks by inverting the time mapping
    smin <- size_range[1]; smax <- size_range[2]
    if (!is.finite(smin) || !is.finite(smax) || smax <= smin) {
      stop("Invalid size_range: must be finite with smax > smin.", call. = FALSE)
    }
    marks_vec <- .time_to_size(used_df$time, smin = smin, smax = smax, delta = delta)
  }

  realization <- data.frame(
    time = used_df$time,
    x = used_df$x,
    y = used_df$y,
    marks = marks_vec
  )

  sim_mpp <- generate_mpp(
    locations = realization[, c("x", "y"), drop = FALSE],
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
    meta = list(
      unthinned = sim_df,
      thinned = sim_thin_df,
      mark_mode = mark_mode,
      mapping = if (mark_mode == "time_to_size") list(delta = delta, size_range = size_range) else NULL
    )
  )
}
