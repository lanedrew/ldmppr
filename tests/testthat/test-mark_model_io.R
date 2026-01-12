test_that("save_mark_model/load_mark_model roundtrip preserves predictions (xgboost)", {
  skip_if_not_installed("xgboost")

  # Load example rasters (same pattern as your docs)
  raster_paths <- list.files(system.file("extdata", package = "ldmppr"),
                             pattern = "\\.tif$", full.names = TRUE)
  raster_paths <- raster_paths[!grepl("_med\\.tif$", raster_paths)]
  rasters <- lapply(raster_paths, terra::rast)
  scaled_raster_list <- scale_rasters(rasters)

  # Example locations
  locations <- small_example_data |>
    dplyr::mutate(time = power_law_mapping(size, 0.5))

  # Train a *tiny* model (fast settings)
  set.seed(1)
  mm <- train_mark_model(
    data = locations,
    raster_list = scaled_raster_list,
    scaled_rasters = TRUE,
    model_type = "xgboost",
    xy_bounds = c(0, 25, 0, 25),
    parallel = FALSE,
    include_comp_inds = FALSE,
    edge_correction = "none",
    selection_metric = "rmse",
    cv_folds = 2,
    tuning_grid_size = 1,
    verbose = FALSE
  )

  # Use a small "sim_realization"-shaped object (x/y/time columns)
  sim_real <- locations[, c("x", "y", "time")]

  p1 <- predict_marks(
    sim_realization = sim_real,
    raster_list = scaled_raster_list,
    scaled_rasters = TRUE,
    mark_model = mm,
    xy_bounds = c(0, 25, 0, 25),
    include_comp_inds = FALSE,
    edge_correction = "none"
  )

  tmp <- tempfile(fileext = ".rds")
  save_mark_model(mm, tmp)
  mm2 <- load_mark_model(tmp)

  p2 <- predict_marks(
    sim_realization = sim_real,
    raster_list = scaled_raster_list,
    scaled_rasters = TRUE,
    mark_model = mm2,
    xy_bounds = c(0, 25, 0, 25),
    include_comp_inds = FALSE,
    edge_correction = "none"
  )

  expect_type(p1, "double")
  expect_length(p1, nrow(sim_real))
  expect_equal(p1, p2, tolerance = 1e-10)
})

test_that("save_mark_model/load_mark_model roundtrip preserves predictions (ranger)", {
  skip_if_not_installed("ranger")

  raster_paths <- list.files(system.file("extdata", package = "ldmppr"),
                             pattern = "\\.tif$", full.names = TRUE)
  raster_paths <- raster_paths[!grepl("_med\\.tif$", raster_paths)]
  rasters <- lapply(raster_paths, terra::rast)
  scaled_raster_list <- scale_rasters(rasters)

  locations <- small_example_data |>
    dplyr::mutate(time = power_law_mapping(size, 0.5))

  set.seed(1)
  mm <- train_mark_model(
    data = locations,
    raster_list = scaled_raster_list,
    scaled_rasters = TRUE,
    model_type = "random_forest",
    xy_bounds = c(0, 25, 0, 25),
    parallel = FALSE,
    include_comp_inds = FALSE,
    edge_correction = "none",
    selection_metric = "rmse",
    cv_folds = 2,
    tuning_grid_size = 1,
    verbose = FALSE
  )

  sim_real <- locations[, c("x", "y", "time")]

  p1 <- predict_marks(
    sim_realization = sim_real,
    raster_list = scaled_raster_list,
    scaled_rasters = TRUE,
    mark_model = mm,
    xy_bounds = c(0, 25, 0, 25),
    include_comp_inds = FALSE,
    edge_correction = "none"
  )

  tmp <- tempfile(fileext = ".rds")
  save_mark_model(mm, tmp)
  mm2 <- load_mark_model(tmp)

  p2 <- predict_marks(
    sim_realization = sim_real,
    raster_list = scaled_raster_list,
    scaled_rasters = TRUE,
    mark_model = mm2,
    xy_bounds = c(0, 25, 0, 25),
    include_comp_inds = FALSE,
    edge_correction = "none"
  )

  expect_type(p1, "double")
  expect_length(p1, nrow(sim_real))
  expect_equal(p1, p2, tolerance = 1e-10)
})
