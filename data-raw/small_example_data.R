# This code generates the small_example_data data object and trained bundled XGBoost model object. Useful for examples and testing code.


# NOTE: since this requires the package to generate the result, we recommend this file
# be run *after* the development version of ldmppr has been installed locally,
# potentially with updated functions. To avoid any conflicts with documentation
# or vignettes which use this object, install the development version with devtools::install(quick = TRUE).
library(ldmppr)

# Read in the raster files
raster_paths <- list.files(system.file("extdata", package = "ldmppr"), pattern = "\\.tif$", full.names = TRUE)
rasters <- lapply(raster_paths, terra::rast)
scaled_raster_list <- scale_rasters(rasters)

# Set a seed for reproducibility
set.seed(90210)

# Set a point to condition on
M_n <- matrix(c(10, 14), ncol = 1)

# Specify the generating parameters of the self-correcting process
generating_parameters <- c(2, 8, .015, 2.5, 3, 1, 2.5, .2)

# Simulate the self-correcting process
generated_locs <- simulate_sc(t_min = 0, t_max = 1, sc_params = generating_parameters, anchor_point = M_n, xy_bounds = c(0, 25, 0, 25))

# Specify the size coefficients
size_coeffs <- c(2, 1, -.5, -1)

# Generate the covariates
covars <- as.matrix(cbind(generated_locs$thinned$time, extract_covars(as.matrix(generated_locs$thinned[, 2:3]), scaled_raster_list)))

# Generate the size values given arrival times and the location specific covariates
size_preds <- 900 * exp(-3 * covars[, 1]) + covars[, -1] %*% size_coeffs + rnorm(nrow(covars), sd = 2)

# Create the example dataset and save it
small_example_data <- generated_locs$thinned %>%
  dplyr::select(x, y) %>%
  dplyr::mutate(size = c(size_preds)) %>%
  dplyr::arrange(dplyr::desc(size))

usethis::use_data(small_example_data, overwrite = TRUE)

# Define the grid values
ngridx <- 10
ngridy <- 10
ngridt <- 10
Xgrid <- seq(0, 25, length.out = ngridx)
Ygrid <- seq(0, 25, length.out = ngridy)
Tgrid <- seq(0, 1, length.out = ngridt)

# Define the parameter initialization values
x_0 <- c(1.5, 8.5, .015, 1.5, 3.2, .75, 3, .08)

# Estimate the parameters in parallel
estimate_demo_params <- estimate_parameters_sc_parallel(
  x_grid = Xgrid,
  y_grid = Ygrid,
  t_grid = Tgrid,
  data = small_example_data,
  delta_values = c(.35, .5, .65, .9, 1),
  parameter_inits = x_0,
  upper_bounds = c(1, 25, 25),
  opt_algorithm = "NLOPT_LN_SBPLX",
  nloptr_options = list(maxeval = 300, xtol_rel = 1e-03),
  verbose = TRUE,
  set_future_plan = TRUE
)

# Generate the time values using the power law mapping with optimal delta
# selected from the parameter estimation above
model_training_data <- small_example_data %>%
  dplyr::mutate(time = power_law_mapping(size, .5))

# Train the model
train_mark_mod <- train_mark_model(
  data = model_training_data,
  raster_list = scaled_raster_list,
  model_type = "xgboost",
  xy_bounds = c(0, 25, 0, 25),
  parallel = TRUE,
  include_comp_inds = TRUE,
  competition_radius = 10,
  correction = "none",
  selection_metric = "rmse"
)

# Save the trained model
example_mark_model <- train_mark_mod$bundled_model
saveRDS(example_mark_model, "inst/extdata/example_mark_model.rds")
