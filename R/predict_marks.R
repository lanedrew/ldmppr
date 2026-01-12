#' Predict values from the mark distribution
#'
#' @param sim_realization a data.frame containing a thinned or unthinned realization from \code{simulate_mpp} (or \code{simulate_sc}).
#' @param raster_list a list of raster objects.
#' @param scaled_rasters \code{TRUE} or \code{FALSE} indicating whether the rasters have been scaled.
#' @param mark_model a mark model object. May be a \code{ldmppr_mark_model} or a legacy model.
#' @param xy_bounds  a vector of domain bounds (2 for x, 2 for y).
#' @param include_comp_inds \code{TRUE} or \code{FALSE} indicating whether to generate and use competition indices as covariates.
#' @param competition_radius distance for competition radius if \code{include_comp_inds} is \code{TRUE}.
#' @param edge_correction type of edge correction to apply (\code{"none"} or \code{"toroidal"}).
#'
#' @return a vector of predicted mark values.
#' @export
#'
#' @examples
#' # Simulate a realization
#' generating_parameters <- c(2, 8, .02, 2.5, 3, 1, 2.5, .2)
#' M_n <- c(10, 14)
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
#' # Predict the mark values
#' predict_marks(
#'   sim_realization = generated_locs$thinned,
#'   raster_list = scaled_raster_list,
#'   scaled_rasters = TRUE,
#'   mark_model = mark_model,
#'   xy_bounds = c(0, 25, 0, 25),
#'   include_comp_inds = TRUE,
#'   competition_radius = 10,
#'   edge_correction = "none"
#' )
#'
predict_marks <- function(sim_realization,
                          raster_list = NULL,
                          scaled_rasters = FALSE,
                          mark_model = NULL,
                          xy_bounds = NULL,
                          include_comp_inds = FALSE,
                          competition_radius = 15,
                          edge_correction = "none") {
  if (!is.data.frame(sim_realization)) stop("Provide a data.frame for sim_realization.", call. = FALSE)
  if (is.null(raster_list) || !is.list(raster_list)) stop("Provide a list of rasters for raster_list.", call. = FALSE)
  if (is.null(mark_model)) stop("Provide a mark model for mark_model.", call. = FALSE)
  if (is.null(xy_bounds) || length(xy_bounds) != 4) stop("Provide xy_bounds = c(a_x, b_x, a_y, b_y).", call. = FALSE)
  if (xy_bounds[2] < xy_bounds[1] || xy_bounds[4] < xy_bounds[3]) stop("Invalid xy_bounds ordering.", call. = FALSE)
  if (!edge_correction %in% c("none", "toroidal")) stop("Provide correction = 'none' or 'toroidal'.", call. = FALSE)
  if (!is.logical(include_comp_inds)) stop("Provide a logical value for include_comp_inds.", call. = FALSE)
  if (!is.logical(scaled_rasters)) stop("Provide a logical value for scaled_rasters.", call. = FALSE)
  if (isTRUE(include_comp_inds) && (is.null(competition_radius) || competition_radius < 0)) {
    stop("Provide a nonnegative competition_radius.", call. = FALSE)
  }

  # Empty input => empty output
  if (nrow(sim_realization) == 0) return(numeric(0))

  # ---- Load mark_model path if needed ----
  if (is.character(mark_model) && length(mark_model) == 1 && file.exists(mark_model)) {
    # Prefer ldmppr loader if available; else fallback
    loader <- get0("load_mark_model", mode = "function")
    mark_model <- if (is.function(loader)) loader(mark_model) else readRDS(mark_model)
  }

  # ---- Legacy bundle support (optional) ----
  if (any(grepl("^bundle", class(mark_model)))) {
    if (!requireNamespace("bundle", quietly = TRUE)) {
      stop("mark_model is a bundled object but the 'bundle' package is not available.", call. = FALSE)
    }
    mark_model <- bundle::unbundle(mark_model)
  }

  # ---- Covariates ----
  s <- sim_realization[, c("x", "y")]
  if (!scaled_rasters) raster_list <- scale_rasters(raster_list)

  X <- extract_covars(locations = s, raster_list = raster_list)

  # Drop terra ID column if present; ensure data.frame
  X <- as.data.frame(X)
  if ("ID" %in% names(X)) X[["ID"]] <- NULL

  # Add required columns
  X$x <- sim_realization$x
  X$y <- sim_realization$y
  X$time <- sim_realization$time

  # Ensure unique names AFTER adding x/y/time to avoid collisions
  names(X) <- make.unique(names(X), sep = "__")

  # ---- Competition indices (optional) ----
  if (isTRUE(include_comp_inds)) {
    n <- nrow(X)

    X$near_nbr_dist <- NA_real_
    X$near_nbr_num  <- NA_integer_
    X$avg_nbr_dist  <- NA_real_
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
        # Exclude self
        d_i <- distance_matrix[i, ]
        d_i[i] <- Inf

        nn_dist <- min(d_i)
        nn_idx <- which(d_i == nn_dist)[1]  # pick one if ties

        close_points <- which(distance_matrix[i, ] < competition_radius & seq_len(n) != i)
        close_times <- X$time[close_points]

        X$near_nbr_dist[i] <- nn_dist
        X$near_nbr_num[i]  <- length(close_points)
        X$avg_nbr_dist[i]  <- if (length(close_points) > 0) mean(d_i[close_points]) else nn_dist

        X$near_nbr_time[i] <- X$time[nn_idx]
        X$near_nbr_time_all[i] <- if (length(close_times) > 0) mean(close_times) else X$time[nn_idx]

        X$near_nbr_time_dist_ratio[i] <- X$near_nbr_time[i] / X$near_nbr_dist[i]
      }
    }
  }

  # ---- Predict ----
  if (inherits(mark_model, "ldmppr_mark_model")) {
    pred <- stats::predict(mark_model, new_data = X)
    if (is.data.frame(pred)) {
      if (".pred" %in% names(pred)) return(as.numeric(pred[[".pred"]]))
      if ("pred" %in% names(pred)) return(as.numeric(pred[["pred"]]))
      if (ncol(pred) == 1) return(as.numeric(pred[[1]]))
      stop("Prediction returned a data.frame with unexpected columns.", call. = FALSE)
    }
    return(as.numeric(pred))
  }

  if (inherits(mark_model, "model_fit") || inherits(mark_model, "workflow")) {
    # pred_tbl <- parsnip::predict(mark_model, new_data = X)
    pred_tbl <- stats::predict(mark_model, new_data = X)
    if (is.data.frame(pred_tbl) && ".pred" %in% names(pred_tbl)) return(as.numeric(pred_tbl[[".pred"]]))
    return(as.numeric(pred_tbl[[1]]))
  }

  stop("Unsupported mark_model type: ", paste(class(mark_model), collapse = ", "), call. = FALSE)
}
