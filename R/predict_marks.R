#' Predict values from the mark distribution
#'
#' @description
#' Prefer using the S3 method \code{predict()} on an \code{ldmppr_mark_model}:
#' \code{predict(mark_model, sim_realization = ..., xy_bounds = ...)}.
#' This wrapper is retained for backward compatibility and is deprecated.
#'
#' @param sim_realization a data.frame containing a thinned or unthinned realization from
#'   \code{simulate_mpp} (or \code{simulate_sc}). Must contain \code{x}, \code{y}, and \code{time}.
#' @param raster_list list of raster objects used for mark prediction.
#' @param scaled_rasters \code{TRUE} or \code{FALSE} indicating whether rasters are already scaled.
#' @param mark_model a mark model object. May be an \code{ldmppr_mark_model}, \code{model_fit}, or \code{workflow}.
#' @param xy_bounds vector of domain bounds as \code{c(a_x, b_x, a_y, b_y)}.
#' @param include_comp_inds \code{TRUE} or \code{FALSE} indicating whether to generate and use competition indices as covariates.
#' @param competition_radius positive numeric distance used when \code{include_comp_inds = TRUE}.
#' @param edge_correction type of edge correction to apply (\code{"none"} or \code{"toroidal"}).
#' @param seed optional nonnegative integer seed for reproducibility.
#'
#' @return a vector of predicted mark values.
#' @export
predict_marks <- function(sim_realization,
                          raster_list = NULL,
                          scaled_rasters = FALSE,
                          mark_model = NULL,
                          xy_bounds = NULL,
                          include_comp_inds = FALSE,
                          competition_radius = 15,
                          edge_correction = "none",
                          seed = NULL) {

  .Deprecated(msg = "`predict_marks()` is deprecated. Use `predict(mark_model, sim_realization = ..., ...)` instead.")

  mark_model <- as_mark_model(mark_model)

  if (inherits(mark_model, "ldmppr_mark_model")) {
    return(stats::predict(
      mark_model,
      sim_realization = sim_realization,
      raster_list = raster_list,
      scaled_rasters = scaled_rasters,
      xy_bounds = xy_bounds,
      include_comp_inds = include_comp_inds,
      competition_radius = competition_radius,
      edge_correction = edge_correction,
      seed = seed
    ))
  }

  .predict_marks_legacy(
    sim_realization = sim_realization,
    raster_list = raster_list,
    scaled_rasters = scaled_rasters,
    mark_model = mark_model,
    xy_bounds = xy_bounds,
    include_comp_inds = include_comp_inds,
    competition_radius = competition_radius,
    edge_correction = edge_correction,
    seed = seed
  )
}


#' @rdname ldmppr-internal
#' @keywords internal
.build_mark_predictors <- function(sim_realization,
                                   raster_list,
                                   scaled_rasters,
                                   xy_bounds,
                                   include_comp_inds,
                                   competition_radius,
                                   edge_correction) {

  if (!is.data.frame(sim_realization)) stop("Provide a data.frame for sim_realization.", call. = FALSE)
  if (nrow(sim_realization) == 0) return(data.frame())

  req_cols <- c("x", "y", "time")
  if (!all(req_cols %in% names(sim_realization))) {
    stop("sim_realization must contain columns: x, y, time.", call. = FALSE)
  }

  if (is.null(raster_list) || !is.list(raster_list)) {
    stop("Provide raster_list, or pass an ldmppr_mark_model that contains stored rasters.", call. = FALSE)
  }
  if (!is.logical(scaled_rasters) || length(scaled_rasters) != 1L) {
    stop("Provide a single logical value for scaled_rasters.", call. = FALSE)
  }
  if (is.null(xy_bounds) || length(xy_bounds) != 4) {
    stop("Provide xy_bounds = c(a_x, b_x, a_y, b_y).", call. = FALSE)
  }
  if (xy_bounds[2] < xy_bounds[1] || xy_bounds[4] < xy_bounds[3]) stop("Invalid xy_bounds ordering.", call. = FALSE)
  if (!edge_correction %in% c("none", "toroidal")) {
    stop("Provide edge_correction = 'none' or 'toroidal'.", call. = FALSE)
  }
  if (!is.logical(include_comp_inds) || length(include_comp_inds) != 1L) {
    stop("Provide a single logical value for include_comp_inds.", call. = FALSE)
  }
  if (isTRUE(include_comp_inds) && (is.null(competition_radius) || !is.numeric(competition_radius) || competition_radius <= 0)) {
    stop("Provide a positive numeric competition_radius when include_comp_inds=TRUE.", call. = FALSE)
  }

  if (!isTRUE(scaled_rasters)) {
    raster_list <- scale_rasters(raster_list)
  }

  s <- sim_realization[, c("x", "y"), drop = FALSE]
  X <- extract_covars(locations = s, raster_list = raster_list)
  X <- as.data.frame(X)
  if ("ID" %in% names(X)) X[["ID"]] <- NULL

  names(X) <- make.unique(names(X), sep = "__")
  X$x <- sim_realization$x
  X$y <- sim_realization$y
  X$time <- sim_realization$time
  names(X) <- make.unique(names(X), sep = "__")

  if (isTRUE(include_comp_inds)) {
    n <- nrow(X)

    X$near_nbr_dist <- NA_real_
    X$near_nbr_num <- NA_integer_
    X$avg_nbr_dist <- NA_real_
    X$near_nbr_time <- NA_real_
    X$near_nbr_time_all <- NA_real_
    X$near_nbr_time_dist_ratio <- NA_real_

    if (n >= 2) {
      s_mat <- as.matrix(s)
      colnames(s_mat) <- c("x", "y")

      distance_matrix <- if (edge_correction == "none") {
        as.matrix(stats::dist(s_mat))
      } else {
        toroidal_dist_matrix_optimized(
          s_mat,
          xy_bounds[2] - xy_bounds[1],
          xy_bounds[4] - xy_bounds[3]
        )
      }

      for (i in seq_len(n)) {
        d_i <- distance_matrix[i, ]
        d_i[i] <- Inf

        nn_dist <- min(d_i)
        nn_idx <- which(d_i == nn_dist)[1]

        close_points <- which(distance_matrix[i, ] < competition_radius & seq_len(n) != i)
        close_times <- X$time[close_points]

        X$near_nbr_dist[i] <- nn_dist
        X$near_nbr_num[i] <- length(close_points)
        X$avg_nbr_dist[i] <- if (length(close_points) > 0) mean(d_i[close_points]) else nn_dist

        X$near_nbr_time[i] <- X$time[nn_idx]
        X$near_nbr_time_all[i] <- if (length(close_times) > 0) mean(close_times) else X$time[nn_idx]
        X$near_nbr_time_dist_ratio[i] <- X$near_nbr_time[i] / X$near_nbr_dist[i]
      }
    }
  }

  X
}


#' @rdname ldmppr-internal
#' @keywords internal
.predict_marks_legacy <- function(sim_realization,
                                  raster_list,
                                  scaled_rasters,
                                  mark_model,
                                  xy_bounds,
                                  include_comp_inds,
                                  competition_radius,
                                  edge_correction,
                                  seed = NULL) {

  if (is.null(mark_model)) stop("Provide a mark model for mark_model.", call. = FALSE)

  if (!is.null(seed)) {
    if (is.na(seed) || seed < 0 || seed != as.integer(seed)) {
      stop("Provide a nonnegative integer `seed`.", call. = FALSE)
    }
    set.seed(as.integer(seed))
  }

  X <- .build_mark_predictors(
    sim_realization = sim_realization,
    raster_list = raster_list,
    scaled_rasters = scaled_rasters,
    xy_bounds = xy_bounds,
    include_comp_inds = include_comp_inds,
    competition_radius = competition_radius,
    edge_correction = edge_correction
  )

  pred <- stats::predict(mark_model, new_data = X)
  if (is.data.frame(pred)) {
    if (".pred" %in% names(pred)) return(as.numeric(pred[[".pred"]]))
    if (ncol(pred) == 1) return(as.numeric(pred[[1]]))
    stop("Prediction returned a data.frame with unexpected columns.", call. = FALSE)
  }
  as.numeric(pred)
}
