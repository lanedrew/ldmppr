#' Train an XGBoost model for the mark distribution
#'
#' @param x a data.frame containing size, x, y, and age
#' @param raster_list the list of raster objects
#'
#' @return a bundled model object
#' @export
#'
train_gbm_mark_model <- function(df, raster_list, save_model = FALSE, save_path){

  s <- base::as.matrix(base::cbind(df$x, df$y))
  raster_trans <- scale_rasters(raster_list)
  covars <- extract_covars(x = s, raster_list = raster_trans)

  X <- covars
  X$x <- s[,1]
  X$y <- s[,2]
  X$age <- df$t

  ## Fit the size model
  mod.data.size <- data.frame(size = df$size, X)

  sprl_split <- rsample::initial_split(
    mod.data.size,
    prop = 0.8,
    strata = size
  )

  preprocessing_recipe <-
    recipes::recipe(size ~ ., data = rsample::training(sprl_split)) %>%
    recipes::prep()


  sprl_cv_folds <-
    recipes::bake(
      preprocessing_recipe,
      new_data = rsample::training(sprl_split)
    ) %>%
    rsample::vfold_cv(v = 5)

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
    resamples = sprl_cv_folds,
    grid = xgboost_grid,
    metrics = yardstick::metric_set(yardstick::rmse, yardstick::rsq, yardstick::mae),
    control = tune::control_grid(verbose = TRUE)
  )

  xgboost_best_params <- xgboost_tuned %>%
    tune::select_best(metric = "mae")

  xgboost_model_final <- xgboost_model %>%
    tune::finalize_model(xgboost_best_params)

  train_processed <- recipes::bake(preprocessing_recipe,  new_data = rsample::training(sprl_split))
  train_prediction <- xgboost_model_final %>%
    # fit the model on all the training data
    parsnip::fit(
      formula = size ~ .,
      data    = train_processed
    ) %>%
    # predict the sale prices for the training data
    stats::predict(new_data = train_processed) %>%
    dplyr::bind_cols(rsample::training(sprl_split))
  xgboost_score_train <-
    train_prediction %>%
    yardstick::metrics(size, .pred) %>%
    dplyr::mutate(.estimate = base::format(base::round(.estimate, 2), big.mark = ","))
  knitr::kable(xgboost_score_train)


  test_processed  <- recipes::bake(preprocessing_recipe, new_data = rsample::testing(sprl_split))
  test_prediction <- xgboost_model_final %>%
    # fit the model on all the training data
    parsnip::fit(
      formula = size ~ .,
      data    = train_processed
    ) %>%
    # use the training model fit to predict the test data
    stats::predict(new_data = test_processed) %>%
    dplyr::bind_cols(rsample::testing(sprl_split))
  # measure the accuracy of our model using `yardstick`
  xgboost_score <-
    test_prediction %>%
    yardstick::metrics(size, .pred) %>%
    dplyr::mutate(.estimate = base::format(base::round(.estimate, 2), big.mark = ","))
  knitr::kable(xgboost_score)

  size.mod <- xgboost_model_final %>%
    # fit the model on all the training data
    parsnip::fit(
      formula = size ~ .,
      data    = mod.data.size
    )

  all_prediction <- size.mod %>%
    # use the training model fit to predict the test data
    stats::predict(new_data = mod.data.size) %>%
    dplyr::bind_cols(mod.data.size)
  # measure the accuracy of our model using `yardstick`
  xgboost_score <-
    all_prediction %>%
    yardstick::metrics(size, .pred) %>%
    dplyr::mutate(.estimate = base::format(base::round(.estimate, 2), big.mark = ","))
  knitr::kable(xgboost_score)

  bundled_mod <- bundle::bundle(size.mod)

  if(save_model == TRUE){
    base::saveRDS(bundled_mod, file = save_path)
  }

  return(bundled_mod)

}
