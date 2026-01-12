test_that("simulate_mpp returns ldmppr_sim", {
  skip_on_cran()

  raster_paths <- list.files(system.file("extdata", package = "ldmppr"),
                             pattern = "\\.tif$", full.names = TRUE)
  raster_paths <- raster_paths[!grepl("_med\\.tif$", raster_paths)]
  rasters <- lapply(raster_paths, terra::rast)
  scaled_raster_list <- scale_rasters(rasters)

  mm <- load_mark_model(system.file("extdata", "example_mark_model.rds", package = "ldmppr"))

  sim <- simulate_mpp(
    process = "self_correcting",
    process_fit = c(2, 8, .02, 2.5, 3, 1, 2.5, .2),
    t_min = 0, t_max = 1,
    anchor_point = c(10, 14),
    raster_list = scaled_raster_list,
    scaled_rasters = TRUE,
    mark_model = mm,
    xy_bounds = c(0, 25, 0, 25),
    include_comp_inds = TRUE,
    edge_correction = "none",
    competition_radius = 10,
    thinning = TRUE
  )

  expect_s3_class(sim, "ldmppr_sim")
  df <- as.data.frame(sim)
  expect_true(all(c("time","x","y","marks") %in% names(df)))
  expect_type(df$marks, "double")
})
