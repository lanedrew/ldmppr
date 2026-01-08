test_that("check_model_fit returns ldmppr_model_check", {
  skip_on_cran()
  skip_if_not_installed("GET")
  skip_if_not_installed("spatstat.geom")
  skip_if_not_installed("spatstat.explore")

  # Minimal run: small n_sim
  # (Use your existing small example data + rasters + saved example mark model)
  file_path <- system.file("extdata", "example_mark_model.rds", package = "ldmppr")
  mm <- load_mark_model(file_path)

  data(small_example_data, package = "ldmppr")
  reference_data <- generate_mpp(
    locations = small_example_data[, c("x", "y")],
    marks = small_example_data$size,
    xy_bounds = c(0, 25, 0, 25)
  )

  raster_paths <- list.files(system.file("extdata", package = "ldmppr"),
                             pattern = "\\.tif$", full.names = TRUE)
  raster_paths <- raster_paths[!grepl("_med\\.tif$", raster_paths)]
  rasters <- lapply(raster_paths, terra::rast)
  scaled_raster_list <- scale_rasters(rasters)

  M_n <- as.numeric(small_example_data[1, c("x", "y")])

  est <- c(1.4, 8.6, 0.02, 1.9, 2.3, 1.1, 2.6, 0.16)

  res <- check_model_fit(
    reference_data = reference_data,
    t_min = 0, t_max = 1,
    sc_params = est,
    anchor_point = M_n,
    raster_list = scaled_raster_list,
    scaled_rasters = TRUE,
    mark_model = mm,
    xy_bounds = c(0, 25, 0, 25),
    include_comp_inds = TRUE,
    thinning = TRUE,
    correction = "none",
    competition_radius = 10,
    n_sim = 100,
    save_sims = FALSE,
    verbose = FALSE,
    seed = 1
  )

  expect_s3_class(res, "ldmppr_model_check")
  expect_true(!is.null(res$combined_env))
  expect_true(all(c("L","F","G","J","E","V") %in% names(res$envs)))
})
