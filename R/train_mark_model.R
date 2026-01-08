#' Train a flexible model for the mark distribution
#'
#' @description
#' Trains a predictive model for the mark distribution of a spatio-temporal process.
#' Allows the user to incorporate location specific information and competition indices as
#' covariates in the mark model.
#'
#'
#' @param data a data frame containing named vectors x, y, size, and time.
#' @param raster_list a list of raster objects.
#' @param scaled_rasters `TRUE` or `FALSE` indicating whether the rasters have been scaled.
#' @param model_type the machine learning model type ("xgboost" or "random_forest").
#' @param xy_bounds a vector of domain bounds (2 for x, 2 for y).
#' @param save_model `TRUE` or `FALSE` indicating whether to save the generated model.
#' @param save_path path for saving the generated model.
#' @param parallel `TRUE` or `FALSE` indicating whether to use parallelization in model training.
#' @param n_cores number of cores to use in parallel model training (if `parallel` is `TRUE`).
#' @param include_comp_inds `TRUE` or `FALSE` indicating whether to generate and use competition indices as covariates.
#' @param competition_radius distance for competition radius if \code{include_comp_inds} is `TRUE`.
#' @param correction type of correction to apply ("none", "toroidal", or "truncation").
#' @param selection_metric metric to use for identifying the optimal model ("rmse" or "mae").
#' @param cv_folds number of cross-validation folds to use in model training.
#' @param tuning_grid_size size of the tuning grid for hyperparameter tuning.
#' @param verbose `TRUE` or `FALSE` indicating whether to show progress of model training.
#'
#' @return an ldmppr_mark_model object.
#' @export
#'
#' @examples
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
#' # Load example locations
#' locations <- small_example_data %>%
#'   dplyr::mutate(time = power_law_mapping(size, .5))
#'
#' # Train the model
#' mark_model <- train_mark_model(
#'   data = locations,
#'   raster_list = scaled_raster_list,
#'   scaled_rasters = TRUE,
#'   model_type = "xgboost",
#'   xy_bounds = c(0, 25, 0, 25),
#'   parallel = FALSE,
#'   include_comp_inds = FALSE,
#'   competition_radius = 10,
#'   correction = "none",
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
                             save_model = FALSE,
                             save_path = NULL,
                             parallel = TRUE,
                             n_cores = NULL,
                             include_comp_inds = FALSE,
                             competition_radius = 15,
                             correction = "none",
                             selection_metric = "rmse",
                             cv_folds = 5,
                             tuning_grid_size = 200,
                             verbose = TRUE) {
  if (!is.data.frame(data)) stop("Provide a data frame of the form (x, y, size, time) for the data argument.", .call = FALSE)
  if (is.null(raster_list) | !is.list(raster_list)) stop("Provide a list of rasters for the raster_list argument.", .call = FALSE)
  if (is.null(xy_bounds) | !(length(xy_bounds) == 4)) stop("Provide (x,y) bounds in the form (a_x, b_x, a_y, b_y) for the xy_bounds argument.", .call = FALSE)
  if (xy_bounds[2] < xy_bounds[1] | xy_bounds[4] < xy_bounds[3]) stop("Provide (x,y) bounds in the form (a_x, b_x, a_y, b_y) for the xy_bounds argument.", .call = FALSE)
  if (isTRUE(save_model) && is.null(save_path)) stop("Provide a path for saving the mark model object.", .call = FALSE)
  if (!correction %in% c("none", "toroidal", "truncation")) stop("Provide a valid correction type for the correction argument.", .call = FALSE)
  if (include_comp_inds == TRUE & (is.null(competition_radius) | competition_radius < 0)) stop("Provide the desired radius for competition_indices argument.", .call = FALSE)
  if (!selection_metric %in% c("rmse", "mae", "rsq")) stop("Provide a valid metric for selection_metric ('rmse', 'mae', or 'rsq').", .call = FALSE)
  if (!model_type %in% c("xgboost", "random_forest")) stop("Provide a valid model type for the model_type argument.", .call = FALSE)
  if (!is.logical(parallel)) stop("Provide a logical value for the parallel argument.", .call = FALSE)
  if (parallel == TRUE & !is.null(n_cores)) {
    if (!is.numeric(n_cores) | n_cores < 1) stop("Provide a numeric value greater than 0 for the n_cores argument.", .call = FALSE)
  }
  if (!is.logical(include_comp_inds)) stop("Provide a logical value for the include_comp_inds argument.", .call = FALSE)
  if (!is.logical(save_model)) stop("Provide a logical value for the save_model argument.", .call = FALSE)
  if (!is.logical(verbose)) stop("Provide a logical value for the verbose argument.", .call = FALSE)
  if (!is.numeric(cv_folds) | cv_folds < 2) stop("Provide a numeric value greater than 1 for the cv_folds argument.", .call = FALSE)
  if (!is.numeric(tuning_grid_size) | tuning_grid_size < 1) stop("Provide a numeric value greater than 0 for the tuning_grid_size argument.", .call = FALSE)
  if (!is.logical(scaled_rasters)) stop("Provide a logical value for the scaled_rasters argument.", .call = FALSE)

  # Initialize parallelization for model training
  cl <- NULL
  if (isTRUE(parallel)) {
    if( !is.null(n_cores)) {
      n_workers <- n_cores
    } else {
      n_workers <- max(1L, floor(parallel::detectCores()/2))
    }
    cl <- parallel::makePSOCKcluster(n_workers)
    doParallel::registerDoParallel(cl)
    on.exit({
      parallel::stopCluster(cl)
      foreach::registerDoSEQ()
    }, add = TRUE)
  } else {
    foreach::registerDoSEQ()
  }

  # Obtain a matrix of (x, y) locations
  s <- base::as.matrix(base::cbind(data$x - xy_bounds[1], data$y - xy_bounds[3]))

  # Scale the rasters if not already scaled
  if (scaled_rasters == FALSE) {
    raster_list <- scale_rasters(raster_list)
  }

  # Obtain the location specific covariate values from the scaled rasters
  X <- extract_covars(locations = s, raster_list = raster_list)
  X$x <- s[, 1]
  X$y <- s[, 2]
  X$time <- data$time

  if (include_comp_inds == TRUE) {
    # Calculate competition indices in a given radius
    X$near_nbr_dist <- NA_real_
    X$near_nbr_num <- NA_real_
    X$avg_nbr_dist <- NA_real_
    X$near_nbr_time <- NA_real_
    X$near_nbr_time_all <- NA_real_
    X$near_nbr_time_dist_ratio <- NA_real_

    # Calculate distance matrices for selected correction method
    if ((correction == "none") | (correction == "truncation")) {
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

  ## Define the model data for fitting the size model
  model_data <- data.frame(size = data$size, X)

  if (correction == "truncation") {
    model_data <- model_data[(model_data$x > 15 &
                                model_data$x < (xy_bounds[2] - xy_bounds[1]) - 15 &
                                model_data$y > 15 &
                                model_data$y < (xy_bounds[4] - xy_bounds[3]) - 15), ]
  }

  if (verbose) {
    base::message("Processing data...")
  }

  # Create the data split
  data_split <- rsample::initial_split(
    data = model_data,
    prop = 0.8,
    strata = size
  )

  # Specify the recipe (retain = TRUE makes the prepped object more robust across sessions)
  preprocessing_recipe <-
    recipes::recipe(size ~ ., data = rsample::training(data_split)) %>%
    recipes::prep(retain = TRUE)

  # Generate the cross validation folds
  data_cv_folds <-
    recipes::bake(
      preprocessing_recipe,
      new_data = rsample::training(data_split)
    ) %>%
    rsample::vfold_cv(v = cv_folds)

  # Define the metric set
  metric_set <- if (selection_metric == "rsq") {
    yardstick::metric_set(yardstick::rmse, yardstick::mae, yardstick::rsq)
  } else {
    yardstick::metric_set(yardstick::rmse, yardstick::mae)
  }

  # Define control object for tuning
  ctrl <- tune::control_grid(
    verbose = FALSE,
    parallel_over = if (isTRUE(parallel)) "resamples" else NULL
  )


  if (model_type == "xgboost") {
    if (verbose) base::message("Training XGBoost model...")

    # Choose threads safely: if we're parallelizing resamples, force xgboost to 1 thread
    nthread <- if (isTRUE(parallel)) 1L else max(1L, parallel::detectCores() - 1L)

    xgboost_model <-
      parsnip::boost_tree(
        mode = "regression",
        trees = hardhat::tune(),
        min_n = hardhat::tune(),
        tree_depth = hardhat::tune(),
        learn_rate = hardhat::tune(),
        loss_reduction = hardhat::tune()
      )

    xgboost_engine_args <- list(
      objective = "reg:squarederror",
      nthread   = nthread,
      verbose   = 0
    )

    xgboost_model <- do.call(
      parsnip::set_engine,
      c(list(object = xgboost_model, engine = "xgboost"), xgboost_engine_args)
    )

    xgboost_params <-
      dials::parameters(
        dials::trees(),
        dials::min_n(),
        dials::tree_depth(),
        dials::learn_rate(),
        dials::loss_reduction()
      )

    xgboost_grid <-
      dials::grid_space_filling(
        xgboost_params,
        size = tuning_grid_size
      )

    xgboost_wf <-
      workflows::workflow() %>%
      workflows::add_model(xgboost_model) %>%
      workflows::add_formula(size ~ .)

    xgboost_tuned <- tune::tune_grid(
      object = xgboost_wf,
      resamples = data_cv_folds,
      grid = xgboost_grid,
      metrics = yardstick::metric_set(yardstick::rmse, yardstick::mae, yardstick::rsq),
      control = tune::control_grid(verbose = FALSE)
    )

    xgboost_best_params <- xgboost_tuned %>%
      tune::select_best(metric = selection_metric)

    xgboost_model_final <- xgboost_model %>%
      tune::finalize_model(xgboost_best_params)

    mark_model <- xgboost_model_final %>%
      parsnip::fit(
        formula = size ~ .,
        data    = model_data
      )

  } else if (model_type == "random_forest") {
    if (verbose) base::message("Training Random Forest model...")

    rf_model <-
      parsnip::rand_forest(
        mode = "regression",
        trees = hardhat::tune(),
        min_n = hardhat::tune()
      )

    rf_engine_args <- list(
      num.threads = if (isTRUE(parallel)) 1L else max(1L, parallel::detectCores() - 1L)
    )

    rf_model <- do.call(
      parsnip::set_engine,
      c(list(object = rf_model, engine = "ranger"), rf_engine_args)
    )

    rf_params <-
      dials::parameters(
        dials::trees(),
        dials::min_n()
      )

    rf_grid <-
      dials::grid_space_filling(
        rf_params,
        size = tuning_grid_size
      )

    rf_wf <-
      workflows::workflow() %>%
      workflows::add_model(rf_model) %>%
      workflows::add_formula(size ~ .)

    rf_tuned <- tune::tune_grid(
      object = rf_wf,
      resamples = data_cv_folds,
      grid = rf_grid,
      metrics = yardstick::metric_set(yardstick::rmse, yardstick::mae, yardstick::rsq),
      control = tune::control_grid(verbose = FALSE)
    )

    rf_best_params <- rf_tuned %>%
      tune::select_best(metric = selection_metric)

    rf_model_final <- rf_model %>%
      tune::finalize_model(rf_best_params)

    mark_model <- rf_model_final %>%
      parsnip::fit(
        formula = size ~ .,
        data    = model_data
      )

  } else {
    stop("Please provide xgboost or random_forest argument for model_type.")
  }

  if (verbose) base::message("Training complete!")

  # Identify engine + extract engine fit
  engine <- if (model_type == "xgboost") "xgboost" else "ranger"
  fit_engine <- mark_model$fit

  # Derive feature_names in baked order (for stable prediction)
  baked_train <- recipes::bake(preprocessing_recipe, new_data = rsample::training(data_split))
  baked_train$size <- NULL
  feature_names <- names(baked_train)

  mm <- ldmppr_mark_model(
    engine = engine,
    fit_engine = fit_engine,
    recipe = preprocessing_recipe,
    outcome = "size",
    feature_names = feature_names,
    info = list(
      model_type = model_type,
      correction = correction,
      include_comp_inds = include_comp_inds,
      cv_folds = cv_folds
    )
  )

  if (isTRUE(save_model)) {
    save_mark_model(mm, save_path)
  }

  return(mm)
}

