#' Train a flexible model for the mark distribution
#'
#' @description
#' Trains a predictive model for the mark distribution of a spatio-temporal process.
#' \code{data} may be either (1) a data.frame containing columns \code{x}, \code{y}, \code{size} and \code{time},
#' (2) a data.frame containing \code{x}, \code{y}, \code{size} (time will be derived via \code{delta}),
#' or (3) a \code{ldmppr_fit} object returned by \code{\link{estimate_process_parameters}}.
#' Allows the user to incorporate location specific information and competition indices as covariates in the mark model.
#'
#' @param data a data.frame or a \code{ldmppr_fit} object. See Description.
#' @param raster_list a list of raster objects.
#' @param scaled_rasters \code{TRUE} or \code{FALSE} indicating whether the rasters have been scaled.
#' @param model_type the machine learning model type (\code{"xgboost"} or \code{"random_forest"}).
#' @param xy_bounds a vector of domain bounds (2 for x, 2 for y). If \code{data} is an \code{ldmppr_fit}
#'   and \code{xy_bounds} is \code{NULL}, defaults to \code{c(0, b_x, 0, b_y)} derived from fit.
#' @param delta (optional) numeric scalar used only when \code{data} contains \code{(x,y,size)} but not \code{time}.
#'   If \code{data} is an \code{ldmppr_fit} and time is missing, the function will infer the \code{delta} value from the fit.
#' @param save_model \code{TRUE} or \code{FALSE} indicating whether to save the generated model.
#' @param save_path path for saving the generated model.
#' @param parallel \code{TRUE} or \code{FALSE} indicating whether to use parallelization in model training.
#' @param n_cores number of cores to use in parallel model training (if \code{parallel} is \code{TRUE}).
#' @param include_comp_inds \code{TRUE} or \code{FALSE} indicating whether to generate and use competition indices as covariates.
#' @param competition_radius distance for competition radius if \code{include_comp_inds} is \code{TRUE}.
#' @param edge_correction type of edge correction to apply (\code{"none"}, \code{"toroidal"}, or \code{"truncation"}).
#' @param selection_metric metric to use for identifying the optimal model (\code{"rmse"}, \code{"mae"}, or \code{"rsq"}).
#' @param cv_folds number of cross-validation folds to use in model training.
#'   If \code{cv_folds <= 1}, tuning is skipped and the model is fit once with default hyperparameters.
#' @param tuning_grid_size size of the tuning grid for hyperparameter tuning.
#' @param verbose \code{TRUE} or \code{FALSE} indicating whether to show progress of model training.
#'
#' @return an object of class \code{"ldmppr_mark_model"} containing the trained mark model.
#' @export
#'
#' @examples
#' # Load the small example data
#' data(small_example_data)
#'
#' # Load example raster data
#' raster_paths <- list.files(system.file("extdata", package = "ldmppr"),
#'   pattern = "\\.tif$", full.names = TRUE
#' )
#' raster_paths <- raster_paths[!grepl("_med\\.tif$", raster_paths)]
#' rasters <- lapply(raster_paths, terra::rast)
#'
#' # Scale the rasters
#' scaled_raster_list <- scale_rasters(rasters)
#'
#'
#' # Train the model
#' mark_model <- train_mark_model(
#'   data = small_example_data,
#'   raster_list = scaled_raster_list,
#'   scaled_rasters = TRUE,
#'   model_type = "xgboost",
#'   xy_bounds = c(0, 25, 0, 25),
#'   delta = 1,
#'   parallel = FALSE,
#'   include_comp_inds = FALSE,
#'   competition_radius = 10,
#'   edge_correction = "none",
#'   selection_metric = "rmse",
#'   cv_folds = 3,
#'   tuning_grid_size = 2,
#'   verbose = TRUE
#' )
#'
#' print(mark_model)
#'
train_mark_model <- function(data,
                             raster_list = NULL,
                             scaled_rasters = FALSE,
                             model_type = "xgboost",
                             xy_bounds = NULL,
                             delta = NULL,
                             save_model = FALSE,
                             save_path = NULL,
                             parallel = TRUE,
                             n_cores = NULL,
                             include_comp_inds = FALSE,
                             competition_radius = 15,
                             edge_correction = "none",
                             selection_metric = "rmse",
                             cv_folds = 5,
                             tuning_grid_size = 200,
                             verbose = TRUE) {

  # -------------------------
  # argument checks
  # -------------------------
  if (is.null(raster_list) || !is.list(raster_list)) {
    stop("Provide a list of rasters for the raster_list argument.", call. = FALSE)
  }
  if (!is.logical(scaled_rasters)) stop("Provide a logical value for scaled_rasters.", call. = FALSE)

  if (!model_type %in% c("xgboost", "random_forest")) {
    stop("Provide a valid model type for model_type ('xgboost' or 'random_forest').", call. = FALSE)
  }
  if (!edge_correction %in% c("none", "toroidal", "truncation")) {
    stop("Provide a valid correction type for correction ('none', 'toroidal', 'truncation').", call. = FALSE)
  }
  if (!selection_metric %in% c("rmse", "mae", "rsq")) {
    stop("Provide a valid metric for selection_metric ('rmse', 'mae', 'rsq').", call. = FALSE)
  }
  if (!is.logical(parallel)) stop("Provide a logical value for parallel.", call. = FALSE)
  if (isTRUE(parallel) && !is.null(n_cores)) {
    if (!is.numeric(n_cores) || n_cores < 1) stop("Provide n_cores >= 1.", call. = FALSE)
  }
  if (!is.logical(include_comp_inds)) stop("Provide a logical value for include_comp_inds.", call. = FALSE)
  if (!is.logical(save_model)) stop("Provide a logical value for save_model.", call. = FALSE)
  if (isTRUE(save_model) && is.null(save_path)) stop("Provide save_path when save_model=TRUE.", call. = FALSE)
  if (!is.logical(verbose)) stop("Provide a logical value for verbose.", call. = FALSE)

  # Coerce training data (handles ldmppr_fit)
  coerced <- .coerce_training_df(data, delta = delta, xy_bounds = xy_bounds)
  df <- coerced$df
  xy_bounds <- coerced$xy_bounds

  if (is.null(xy_bounds) || length(xy_bounds) != 4) {
    stop("Provide xy_bounds = c(a_x, b_x, a_y, b_y), or supply an ldmppr_fit with grid$upper_bounds.", call. = FALSE)
  }
  if (xy_bounds[2] < xy_bounds[1] || xy_bounds[4] < xy_bounds[3]) {
    stop("Provide xy_bounds in the form (a_x, b_x, a_y, b_y) with b_x >= a_x and b_y >= a_y.", call. = FALSE)
  }

  cv_folds <- as.integer(cv_folds)
  if (is.na(cv_folds) || cv_folds < 1L) stop("cv_folds must be >= 1.", call. = FALSE)

  tuning_grid_size <- as.integer(tuning_grid_size)
  if (is.na(tuning_grid_size) || tuning_grid_size < 1L) stop("tuning_grid_size must be >= 1.", call. = FALSE)

  # -------------------------
  # parallel backend (foreach)
  # -------------------------
  cl <- NULL
  if (isTRUE(parallel)) {
    n_workers <- if (!is.null(n_cores)) as.integer(n_cores) else max(1L, floor(parallel::detectCores() / 2))
    cl <- parallel::makePSOCKcluster(n_workers)
    doParallel::registerDoParallel(cl)
    on.exit({
      parallel::stopCluster(cl)
      foreach::registerDoSEQ()
    }, add = TRUE)
  } else {
    foreach::registerDoSEQ()
  }

  # -------------------------
  # raster covariates
  # -------------------------
  if (!isTRUE(scaled_rasters)) {
    raster_list <- scale_rasters(raster_list)
  }

  # IMPORTANT: keep x,y in the SAME coordinate system the user provides;
  # do NOT shift by lower bounds unless you do it everywhere (training + predict).
  s <- as.matrix(df[, c("x", "y")])

  X <- extract_covars(locations = s, raster_list = raster_list)

  # terra::extract commonly returns ID; drop and repair names to avoid tibble collisions
  if ("ID" %in% names(X)) {
    X <- X[, names(X) != "ID", drop = FALSE]
  }
  names(X) <- make.unique(names(X), sep = "__")

  X$x <- df$x
  X$y <- df$y
  X$time <- df$time

  if (isTRUE(include_comp_inds)) {
    X$near_nbr_dist <- NA_real_
    X$near_nbr_num <- NA_real_
    X$avg_nbr_dist <- NA_real_
    X$near_nbr_time <- NA_real_
    X$near_nbr_time_all <- NA_real_
    X$near_nbr_time_dist_ratio <- NA_real_

    # distance matrix
    if (edge_correction %in% c("none", "truncation")) {
      distance_matrix <- as.matrix(stats::dist(s, method = "euclidean"))
    } else {
      distance_matrix <- toroidal_dist_matrix_optimized(
        s,
        xy_bounds[2] - xy_bounds[1],
        xy_bounds[4] - xy_bounds[3]
      )
    }

    for (i in seq_len(nrow(X))) {
      close_points <- unique(which(distance_matrix[i, ] < competition_radius & distance_matrix[i, ] != 0))
      close_times <- X$time[close_points]

      X$near_nbr_dist[i] <- min(distance_matrix[i, ][-i])
      X$near_nbr_num[i] <- length(close_points)
      X$avg_nbr_dist[i] <- if (length(close_points)) mean(distance_matrix[i, close_points]) else min(distance_matrix[i, ][-i])

      nn_idx <- unique(which(distance_matrix[i, ] == X$near_nbr_dist[i]))
      X$near_nbr_time[i] <- X$time[nn_idx][1]
      X$near_nbr_time_all[i] <- if (length(close_points)) mean(close_times) else X$time[nn_idx][1]
      X$near_nbr_time_dist_ratio[i] <- X$near_nbr_time[i] / X$near_nbr_dist[i]
    }
  }

  model_data <- data.frame(size = df$size, X)

  if (edge_correction == "truncation") {
    # Apply truncation in ORIGINAL coordinate scale
    ax <- xy_bounds[1]; bx <- xy_bounds[2]
    ay <- xy_bounds[3]; by <- xy_bounds[4]
    model_data <- model_data[
      model_data$x > (ax + 15) &
        model_data$x < (bx - 15) &
        model_data$y > (ay + 15) &
        model_data$y < (by - 15),
      ,
      drop = FALSE
    ]
  }

  if (nrow(model_data) < 2) stop("Not enough observations to train a model after filtering.", call. = FALSE)

  if (isTRUE(verbose)) message("Processing data...")

  # -------------------------
  # tidymodels spec
  # -------------------------
  metric_set <- if (selection_metric == "rsq") {
    yardstick::metric_set(yardstick::rmse, yardstick::mae, yardstick::rsq)
  } else {
    yardstick::metric_set(yardstick::rmse, yardstick::mae)
  }

  ctrl <- tune::control_grid(
    verbose = FALSE,
    parallel_over = if (isTRUE(parallel)) "resamples" else NULL
  )

  # If parallel resampling, force engine threads to 1 to avoid nested parallelism.
  engine_threads <- if (isTRUE(parallel)) 1L else max(1L, parallel::detectCores() - 1L)

  # Use an unprepped recipe in the workflow; we'll prep once at the end for storage.
  recipe_spec <- recipes::recipe(size ~ ., data = model_data)

  do_tuning <- (cv_folds >= 2L) && (nrow(model_data) >= cv_folds)

  if (!do_tuning && isTRUE(verbose)) {
    message("Skipping CV tuning (cv_folds <= 1 or not enough rows). Fitting a single model with default hyperparameters.")
  }

  if (model_type == "xgboost") {
    if (isTRUE(verbose)) message("Training XGBoost model...")

    spec <- parsnip::boost_tree(
      mode = "regression",
      trees = if (do_tuning) hardhat::tune() else 500,
      min_n = if (do_tuning) hardhat::tune() else 5,
      tree_depth = if (do_tuning) hardhat::tune() else 6,
      learn_rate = if (do_tuning) hardhat::tune() else 0.05,
      loss_reduction = if (do_tuning) hardhat::tune() else 0
    ) %>%
      parsnip::set_engine(
        "xgboost",
        objective = "reg:squarederror",
        nthread = engine_threads,
        verbose = 0
      )

    wf <- workflows::workflow() %>%
      workflows::add_model(spec) %>%
      workflows::add_recipe(recipe_spec)

    if (do_tuning) {
      folds <- rsample::vfold_cv(model_data, v = cv_folds)

      params <- dials::parameters(
        dials::trees(),
        dials::min_n(),
        dials::tree_depth(),
        dials::learn_rate(),
        dials::loss_reduction()
      )

      grid <- dials::grid_space_filling(params, size = tuning_grid_size)

      tuned <- tune::tune_grid(
        object = wf,
        resamples = folds,
        grid = grid,
        metrics = metric_set,
        control = ctrl
      )

      best <- tune::select_best(tuned, metric = selection_metric)
      wf_final <- tune::finalize_workflow(wf, best)
      wf_fit <- parsnip::fit(wf_final, data = model_data)
    } else {
      wf_fit <- parsnip::fit(wf, data = model_data)
    }

    fit_engine <- workflows::extract_fit_engine(wf_fit)
    engine <- "xgboost"

  } else {
    if (isTRUE(verbose)) message("Training Random Forest model...")

    spec <- parsnip::rand_forest(
      mode = "regression",
      trees = if (do_tuning) hardhat::tune() else 500,
      min_n = if (do_tuning) hardhat::tune() else 5
    ) %>%
      parsnip::set_engine("ranger", num.threads = engine_threads)

    wf <- workflows::workflow() %>%
      workflows::add_model(spec) %>%
      workflows::add_recipe(recipe_spec)

    if (do_tuning) {
      folds <- rsample::vfold_cv(model_data, v = cv_folds)

      params <- dials::parameters(dials::trees(), dials::min_n())
      grid <- dials::grid_space_filling(params, size = tuning_grid_size)

      tuned <- tune::tune_grid(
        object = wf,
        resamples = folds,
        grid = grid,
        metrics = metric_set,
        control = ctrl
      )

      best <- tune::select_best(tuned, metric = selection_metric)
      wf_final <- tune::finalize_workflow(wf, best)
      wf_fit <- parsnip::fit(wf_final, data = model_data)
    } else {
      wf_fit <- parsnip::fit(wf, data = model_data)
    }

    fit_engine <- workflows::extract_fit_engine(wf_fit)
    engine <- "ranger"
  }

  if (isTRUE(verbose)) message("Training complete!")

  # Prep recipe for storage + stable feature order during predict()
  preprocessing_recipe <- recipes::prep(recipe_spec, training = model_data, retain = TRUE)
  baked <- recipes::bake(preprocessing_recipe, new_data = model_data)
  baked$size <- NULL
  feature_names <- names(baked)

  mm <- ldmppr_mark_model(
    engine = engine,
    fit_engine = fit_engine,
    recipe = preprocessing_recipe,
    outcome = "size",
    feature_names = feature_names,
    info = list(
      model_type = model_type,
      edge_correction = edge_correction,
      include_comp_inds = include_comp_inds,
      competition_radius = competition_radius,
      selection_metric = selection_metric,
      cv_folds = cv_folds,
      tuning_grid_size = tuning_grid_size
    )
  )

  if (isTRUE(save_model)) {
    save_mark_model(mm, save_path)
  }

  mm
}
