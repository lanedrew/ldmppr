# This code generates:
#  (1) the small_example_data dataset
#  (2) an example saved mark model in inst/extdata/example_mark_model.rds
#
# Run after installing the dev version of ldmppr locally.
# devtools::install(quick = TRUE)

library(ldmppr)
library(dplyr)
library(terra)
library(usethis)

# ------------------------------------------------------------
# 0) Setup
# ------------------------------------------------------------
set.seed(90210)

xy_bounds <- c(0, 25, 0, 25)
anchor_point <- c(10, 14)  # prefer vector; avoids matrix-shape ambiguity

# Read in raster files shipped with the package
raster_paths <- list.files(
  system.file("extdata", package = "ldmppr"),
  pattern = "\\.tif$",
  full.names = TRUE
)

# If you have "_med" rasters you want to exclude (as in vignette), uncomment:
raster_paths <- raster_paths[!grepl("_med\\.tif$", raster_paths)]

rasters <- lapply(raster_paths, terra::rast)

# Scale rasters once
scaled_raster_list <- scale_rasters(rasters)

# ------------------------------------------------------------
# 1) Simulate a self-correcting process (locations + times)
# ------------------------------------------------------------
# Generating parameters (alpha1, beta1, gamma1, alpha2, beta2, alpha3, beta3, gamma3)
generating_parameters <- c(2, 8, 0.015, 2.5, 3, 1, 2.5, 0.2)

generated_locs <- simulate_sc(
  t_min = 0,
  t_max = 1,
  sc_params = generating_parameters,
  anchor_point = anchor_point,
  xy_bounds = xy_bounds
)

# We'll use the thinned realization as the example reference data
sim_df <- generated_locs$thinned
stopifnot(all(c("time", "x", "y") %in% names(sim_df)))

# ------------------------------------------------------------
# 2) Create synthetic marks ("size") from time + covariates
# ------------------------------------------------------------
# Extract covariates at simulated locations
Xcov <- extract_covars(
  locations = as.matrix(sim_df[, c("x", "y")]),
  raster_list = scaled_raster_list
) %>%
  as.matrix()


# Coefficients for synthetic mark generation (intercept handled separately)
size_coeffs <- c(2, 1, -.5, -1)

# Synthetic mean structure: decays with time + linear in covariates
mu <- 900 * exp(-3 * sim_df$time) + drop(Xcov %*% size_coeffs)

size_vals <- mu + rnorm(nrow(sim_df), sd = 2)

small_example_data <- sim_df %>%
  select(x, y) %>%
  mutate(size = as.numeric(size_vals)) %>%
  arrange(desc(size))

# Save the dataset into data/
usethis::use_data(small_example_data, overwrite = TRUE)

# ------------------------------------------------------------
# 3) Estimate process parameters (optional, but useful for examples)
# ------------------------------------------------------------
# Grid for likelihood approximation (keep modest for data-generation script)
Xgrid <- seq(0, 25, length.out = 30)
Ygrid <- seq(0, 25, length.out = 30)
Tgrid <- seq(0, 1,  length.out = 30)

parameter_inits <- c(1.5, 8.5, 0.015, 1.5, 3.2, 0.75, 3, 0.08)
upper_bounds <- c(1, 25, 25)

# Use the new unified estimator. Supply delta_values to perform a delta search.
fit_sc <- estimate_process_parameters(
  data = small_example_data,
  process = "self_correcting",
  x_grid = Xgrid,
  y_grid = Ygrid,
  t_grid = Tgrid,
  upper_bounds = upper_bounds,
  parameter_inits = parameter_inits,
  delta_values = c(.35, .5, .65, .9, 1),  # search over these delta values
  parallel = FALSE,
  strategy = c("global_local"),
  global_algorithm = "NLOPT_GN_CRS2_LM",
  local_algorithm = "NLOPT_LN_BOBYQA",
  global_options = list(maxeval = 100),
  local_options = list(maxeval = 300, xtol_rel = 1e-5, maxtime = NULL),
  verbose = TRUE
)

sc_params_hat <- coef(fit_sc)

# Pull the delta used
delta_hat <- fit_sc$mapping$delta

# ------------------------------------------------------------
# 4) Train the mark model using the chosen delta mapping
# ------------------------------------------------------------
model_training_data <- small_example_data %>%
  mutate(time = power_law_mapping(size, delta = delta_hat))

mark_model <- train_mark_model(
  data = model_training_data,
  raster_list = scaled_raster_list,
  scaled_rasters = TRUE,
  model_type = "xgboost",
  xy_bounds = xy_bounds,
  parallel = TRUE,
  include_comp_inds = TRUE,
  competition_radius = 10,
  edge_correction = "none",
  selection_metric = "rmse",
  cv_folds = 5,
  tuning_grid_size = 20,
  verbose = TRUE
)


# ------------------------------------------------------------
# 5) Save the mark model for use in examples/vignettes
# ------------------------------------------------------------
out_path <- file.path("inst", "extdata", "example_mark_model.rds")
dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)

save_mark_model(mark_model, out_path)

message("Wrote: ", out_path)
