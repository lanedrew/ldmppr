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
#' @param model_type the machine learning model type ("xgboost" or "random_forest").
#' @param xy_bounds a vector of domain bounds (2 for x, 2 for y).
#' @param save_model `TRUE` or `FALSE` indicating whether to save the generated model.
#' @param save_path path for saving the generated model.
#' @param parallel `TRUE` or `FALSE` indicating whether to use parallelization in model training.
#' @param include_comp_inds `TRUE` or `FALSE` indicating whether to generate and use competition indices as covariates.
#' @param competition_radius distance for competition radius if \code{include_comp_inds} is `TRUE`.
#' @param correction type of correction to apply ("none", "toroidal", or "truncation").
#' @param verbose `TRUE` or `FALSE` indicating whether to show progress of model training.
#'
#' @return a bundled model object.
#' @export
#'
train_mark_model <- function(data,
                             raster_list = NULL,
                             model_type = "xgboost",
                             xy_bounds = NULL,
                             save_model = FALSE,
                             save_path = NULL,
                             parallel = TRUE,
                             include_comp_inds = FALSE,
                             competition_radius = 15,
                             correction = "none",
                             verbose = TRUE){

  if(!is.data.frame(data)) stop("Provide a data frame of the form (x, y, size, time) for the data argument.", .call = FALSE)
  if(is.null(raster_list) | !is.list(raster_list)) stop("Provide a list of rasters for the raster_list argument.", .call = FALSE)
  if(is.null(xy_bounds) | !(length(xy_bounds) == 4)) stop("Provide (x,y) bounds in the form (a_x, b_x, a_y, b_y) for the xy_bounds argument.", .call = FALSE)
  if(xy_bounds[2] < xy_bounds[1] | xy_bounds[4] < xy_bounds[3]) stop("Provide (x,y) bounds in the form (a_x, b_x, a_y, b_y) for the xy_bounds argument.", .call = FALSE)
  if(save_model == TRUE & is.null(save_path)) stop("Provide a path for saving the bundled model object.", .call = FALSE)
  if(!correction %in% c("none", "toroidal", "truncation")) stop("Provide a valid correction type for the correction argument.", .call = FALSE)
  if(include_comp_inds == TRUE & (is.null(competition_radius) | competition_radius < 0)) stop("Provide the desired radius for competition_indices argument.", .call = FALSE)


  # Initialize parallelization for model training
  if(parallel) {
    doParallel::registerDoParallel()
  }

  # Obtain a matrix of (x, y) locations
  s <- base::as.matrix(base::cbind(data$x - xy_bounds[1], data$y - xy_bounds[3]))

  # Scale the provided rasters in the raster_list
  raster_trans <- scale_rasters(raster_list)

  # Obtain the location specific covariate values from the scaled rasters
  X <- extract_covars(x = s, raster_list = raster_trans)
  X$x <- s[,1]
  X$y <- s[,2]
  X$time <- data$time

  if(include_comp_inds == TRUE){

    # Calculate competition indices in a 15 unit radius
    X$near_nbr_dist <- NA
    X$near_nbr_num <- NA
    X$avg_nbr_dist <- NA
    X$near_nbr_time <- NA
    X$near_nbr_time_all <- NA
    X$near_nbr_time_dist_ratio <- NA
    colnames(s) <- c("x", "y")

    # Calculate distance matrices for selected correction method
    if((correction == "none") | (correction == "truncation")) {
      distance_matrix <- base::as.matrix(stats::dist(s, method = "euclidean"))
    }else if(correction == "toroidal") {
      distance_matrix <- toroidal_dist_matrix_optimized(s, xy_bounds[2] - xy_bounds[1], xy_bounds[4] - xy_bounds[3])
    }


    for(i in 1:base::nrow(X)){
      close_points <- base::unique(base::which(distance_matrix[i,] < competition_radius & distance_matrix[i,] != 0))
      close_times <- X$time[close_points]
      X$near_nbr_dist[i] <- base::min(distance_matrix[i,][-i])
      X$near_nbr_num[i] <- base::length(close_points)
      X$avg_nbr_dist[i] <- base::mean(distance_matrix[i,][close_points])
      if(base::length(close_points) == 0){
        X$avg_nbr_dist[i] <- base::min(distance_matrix[i,][-i])
      }
      X$near_nbr_time[i] <- X$time[base::unique(base::which(distance_matrix[i,] == X$near_nbr_dist[i]))]
      X$nea_.nbr_time_all[i] <- mean(close_times)
      if(base::length(close_points) == 0){
        X$near_nbr_time_all[i] <- X$time[base::unique(base::which(distance_matrix[i,] == X$near_nbr_dist[i]))]
      }
      X$near_nbr_time_dist_ratio[i] <- X$near_nbr_time[i]/X$near_nbr_dist[i]
    }
  }


  ## Define the model data for fitting the size model
  model_data <- data.frame(size = data$size, X)

  if(correction == "truncation"){
    model_data <- model_data[(model_data$x > 15 &
                              model_data$x < (xy_bounds[2] - xy_bounds[1]) - 15 &
                              model_data$y > 15 &
                              model_data$y < (xy_bounds[4] - xy_bounds[3]) - 15),]
  }

  if(verbose) {
    base::message("Processing data...")
  }

  # Create the data split
  data_split <- rsample::initial_split(
    data = model_data,
    prop = 0.8,
    strata = size
  )

  # Specify the recipe
  preprocessing_recipe <-
    recipes::recipe(size ~ ., data = rsample::training(data_split)) %>%
    recipes::prep()

  # Generate the cross validation folds
  data_cv_folds <-
    recipes::bake(
      preprocessing_recipe,
      new_data = rsample::training(data_split)
    ) %>%
    rsample::vfold_cv(v = 5)


  if(model_type == "xgboost"){

    if(verbose) {
      base::message("Training XGBoost model...")
    }

    # Specify the XGBoost model
    xgboost_model <-
      parsnip::boost_tree(
        mode = "regression",
        trees = 1000,
        min_n = hardhat::tune(),
        tree_depth = hardhat::tune(),
        learn_rate = hardhat::tune(),
        loss_reduction = hardhat::tune()
      ) %>%
      parsnip::set_engine("xgboost", objective = "reg:squarederror")

    # Specify the XGBoost hyperparameters for tuning
    xgboost_params <-
      dials::parameters(
        dials::min_n(),
        dials::tree_depth(),
        dials::learn_rate(),
        dials::loss_reduction()
      )

    # Create the tuning grid
    xgboost_grid <-
      dials::grid_space_filling(
        xgboost_params,
        size = 200
      )

    # Establish the model fitting workflow
    xgboost_wf <-
      workflows::workflow() %>%
      workflows::add_model(xgboost_model) %>%
      workflows::add_formula(size ~ .)

    # Perform the hyperparameter tuning
    xgboost_tuned <- tune::tune_grid(
      object = xgboost_wf,
      resamples = data_cv_folds,
      grid = xgboost_grid,
      metrics = yardstick::metric_set(yardstick::rmse, yardstick::rsq, yardstick::mae),
      control = tune::control_grid(verbose = FALSE)
    )

    # Obtain the best model by "mae"
    xgboost_best_params <- xgboost_tuned %>%
      tune::select_best(metric = "mae")

    # Obtain the finalized model
    xgboost_model_final <- xgboost_model %>%
      tune::finalize_model(xgboost_best_params)

    # Fit the model on the full training data
    mark_model <- xgboost_model_final %>%
      parsnip::fit(
        formula = size ~ .,
        data    = model_data
      )
  }else if(model_type == "random_forest"){

    if(verbose) {
      base::message("Training Random Forest model...")
    }
    # Specify the Ranger random forest model
    rf_model <-
      parsnip::rand_forest(
        mode = "regression",
        trees = hardhat::tune(),
        min_n = hardhat::tune()
      ) %>%
      parsnip::set_engine("ranger")

    # Specify the random forest hyperparameters for tuning
    rf_params <-
      dials::parameters(
        dials::trees(),
        dials::min_n()
      )

    # Create the tuning grid
    rf_grid <-
      dials::grid_space_filling(
        rf_params,
        size = 200
      )

    # Establish the model fitting workflow
    rf_wf <-
      workflows::workflow() %>%
      workflows::add_model(rf_model) %>%
      workflows::add_formula(size ~ .)

    # Perform the hyperparameter tuning
    rf_tuned <- tune::tune_grid(
      object = rf_wf,
      resamples = data_cv_folds,
      grid = rf_grid,
      metrics = yardstick::metric_set(yardstick::rmse, yardstick::rsq, yardstick::mae),
      control = tune::control_grid(verbose = FALSE)
    )

    # Obtain the best model by "mae"
    rf_best_params <- rf_tuned %>%
      tune::select_best(metric = "mae")

    # Obtain the finalized model
    rf_model_final <- rf_model %>%
      tune::finalize_model(rf_best_params)

    # Fit the model on the full training data
    mark_model <- rf_model_final %>%
      parsnip::fit(
        formula = size ~ .,
        data    = model_data
      )
  }else {
    stop("Please provide xgboost or random_forest argument for model_type.")
  }

  if(verbose) {
    base::message("Training complete!")
  }

  # Bundle the trained model object
  bundled_mod <-
    mark_model %>%
    bundle::bundle()

  if(save_model == TRUE){
    # Save the bundled model object
    base::saveRDS(bundled_mod, file = save_path)
  }

  # Return a list containing the raw model and bundled model
  return(base::list(raw_model = mark_model, bundled_model = bundled_mod))

}
