#' Train an ML model for the mark distribution
#'
#' @param df a data.frame containing named vectors size, x, y, and time
#' @param raster_list the list of raster objects
#' @param model_type the machine learning model type
#' @param save_model determines whether to save the generated model
#' @param save_path path for saving the generated model
#' @param parallel TRUE or FALSE indicating whether to use parallelization in model training
#' @param verbose TRUE or FALSE indicating whether to show progress of model training
#'
#' @return a bundled model object
#' @export
#'
train_mark_model <- function(df, raster_list, model_type = "xgboost",
                             save_model = FALSE, save_path,
                             parallel = TRUE, verbose = TRUE){

  if(parallel) {
    doParallel::registerDoParallel()
  }

  s <- base::as.matrix(base::cbind(df$x, df$y))
  raster_trans <- scale_rasters(raster_list)

  X <- extract_covars(x = s, raster_list = raster_trans)
  X$x <- s[,1]
  X$y <- s[,2]
  X$time <- df$time

  ## Fit the size model
  model_data <- data.frame(size = df$size, X)

  if(verbose) {
    base::print("Processing data...")
  }

  data_split <- rsample::initial_split(
    data = model_data,
    prop = 0.8,
    strata = size
  )

  preprocessing_recipe <-
    recipes::recipe(size ~ ., data = rsample::training(data_split)) %>%
    recipes::prep()


  data_cv_folds <-
    recipes::bake(
      preprocessing_recipe,
      new_data = rsample::training(data_split)
    ) %>%
    rsample::vfold_cv(v = 5)


  if(model_type == "xgboost"){

    if(verbose) {
      base::print("Training XGBoost model...")
    }
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

    xgboost_params <-
      dials::parameters(
        dials::min_n(),
        dials::tree_depth(),
        dials::learn_rate(),
        dials::loss_reduction()
      )

    xgboost_grid <-
      dials::grid_max_entropy(
        xgboost_params,
        size = 200
      )

    xgboost_wf <-
      workflows::workflow() %>%
      workflows::add_model(xgboost_model) %>%
      workflows::add_formula(size ~ .)

    xgboost_tuned <- tune::tune_grid(
      object = xgboost_wf,
      resamples = data_cv_folds,
      grid = xgboost_grid,
      metrics = yardstick::metric_set(yardstick::rmse, yardstick::rsq, yardstick::mae),
      control = tune::control_grid(verbose = FALSE)
    )

    xgboost_best_params <- xgboost_tuned %>%
      tune::select_best(metric = "mae")

    xgboost_model_final <- xgboost_model %>%
      tune::finalize_model(xgboost_best_params)

    mark_model <- xgboost_model_final %>%
      # fit the model on all the training data
      parsnip::fit(
        formula = size ~ .,
        data    = model_data
      )
  }else if(model_type == "random_forest"){

    if(verbose) {
      base::print("Training Random Forest model...")
    }
    rf_model <-
      parsnip::rand_forest(
        mode = "regression",
        # mtry = hardhat::tune(),
        trees = hardhat::tune(),
        min_n = hardhat::tune()
      ) %>%
      parsnip::set_engine("ranger")

    rf_params <-
      dials::parameters(
        # dials::mtry(c(1, base::ncol(model_data) - 1)),
        dials::trees(),
        dials::min_n()
      )

    rf_grid <-
      dials::grid_max_entropy(
        rf_params,
        size = 200
      )

    rf_wf <-
      workflows::workflow() %>%
      workflows::add_model(rf_model) %>%
      workflows::add_formula(size ~ .)

    rf_tuned <- tune::tune_grid(
      object = rf_wf,
      resamples = data_cv_folds,
      grid = rf_grid,
      metrics = yardstick::metric_set(yardstick::rmse, yardstick::rsq, yardstick::mae),
      control = tune::control_grid(verbose = FALSE)
    )

    rf_best_params <- rf_tuned %>%
      tune::select_best(metric = "mae")

    rf_model_final <- rf_model %>%
      tune::finalize_model(rf_best_params)

    mark_model <- rf_model_final %>%
      # fit the model on all the training data
      parsnip::fit(
        formula = size ~ .,
        data    = model_data
      )
  }else {
    return(print("Please input xboost or random_forest for model_type."))
  }

  if(verbose) {
    base::print("Training complete!")
  }

  bundled_mod <-
    mark_model %>%
    # butcher::butcher() %>%
    bundle::bundle()

  if(save_model == TRUE){
    base::saveRDS(bundled_mod, file = save_path)
  }

  return(base::list(raw_model = mark_model, bundled_model = bundled_mod))

}
