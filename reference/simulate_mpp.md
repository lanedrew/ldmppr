# Simulate a realization of a location dependent marked point process

Simulate a realization of a location dependent marked point process

## Usage

``` r
simulate_mpp(
  process = c("self_correcting"),
  process_fit = NULL,
  t_min = 0,
  t_max = 1,
  anchor_point = NULL,
  raster_list = NULL,
  scaled_rasters = FALSE,
  mark_model = NULL,
  xy_bounds = NULL,
  include_comp_inds = FALSE,
  competition_radius = 15,
  edge_correction = "none",
  thinning = TRUE,
  seed = NULL,
  mark_mode = c("mark_model", "time_to_size"),
  size_range = NULL,
  delta = NULL
)
```

## Arguments

- process:

  type of process used (currently supports `"self_correcting"`).

- process_fit:

  either (1) a `ldmppr_fit` object returned by
  [`estimate_process_parameters`](https://lanedrew.github.io/ldmppr/reference/estimate_process_parameters.md),
  or (2) a numeric vector of length 8 giving the self-correcting process
  parameters:
  \\(\alpha_1,\beta_1,\gamma_1,\alpha_2,\beta_2,\alpha_3,\beta_3,\gamma_3)\\
  (alpha_1, beta_1, gamma_1, alpha_2, beta_2, alpha_3, beta_3, gamma_3).

- t_min:

  minimum value for time.

- t_max:

  maximum value for time.

- anchor_point:

  (optional) vector of (x,y) coordinates of the point to condition on.
  If `NULL`, inferred from the reference data (largest mark if
  available) or from `process_fit$data_original` (largest size).

- raster_list:

  (optional) list of raster objects used for mark prediction. Required
  when `mark_mode='mark_model'` unless rasters are stored in
  `mark_model`.

- scaled_rasters:

  `TRUE` or `FALSE` indicating whether rasters are already scaled.
  Ignored when `mark_mode='time_to_size'`.

- mark_model:

  a mark model object used when `mark_mode='mark_model'`. May be an
  `ldmppr_mark_model`, `model_fit`, or `workflow`.

- xy_bounds:

  (optional) vector of bounds as `c(a_x, b_x, a_y, b_y)`. If `NULL`,
  bounds are inferred from `process_fit` when available.

- include_comp_inds:

  `TRUE` or `FALSE` indicating whether to compute competition indices.

- competition_radius:

  positive numeric distance used when `include_comp_inds = TRUE`.

- edge_correction:

  type of edge correction to apply (`"none"` or `"toroidal"`).

- thinning:

  `TRUE` or `FALSE` indicating whether to use the thinned simulated
  values.

- seed:

  integer seed for reproducibility.

- mark_mode:

  mark generation mode: `"mark_model"` uses
  [`predict()`](https://rdrr.io/r/stats/predict.html) on a mark model,
  while `"time_to_size"` maps simulated times back to sizes via `delta`.

- size_range:

  numeric vector `c(smin, smax)` used for `mark_mode='time_to_size'`. If
  `NULL`, inferred from `process_fit` when possible.

- delta:

  positive scalar used for `mark_mode='time_to_size'`. If `NULL`,
  inferred from `process_fit` when possible.

## Value

an object of class `"ldmppr_sim"`.

## Examples

``` r
# Specify the generating parameters of the self-correcting process
generating_parameters <- c(2, 8, .02, 2.5, 3, 1, 2.5, .2)

# Specify an anchor point
M_n <- c(10, 14)

# Load the raster files
raster_paths <- list.files(system.file("extdata", package = "ldmppr"),
  pattern = "\\.tif$", full.names = TRUE
)
raster_paths <- raster_paths[!grepl("_med\\.tif$", raster_paths)]
rasters <- lapply(raster_paths, terra::rast)

# Scale the rasters
scaled_raster_list <- scale_rasters(rasters)

# Load the example mark model
file_path <- system.file("extdata", "example_mark_model.rds", package = "ldmppr")
mark_model <- load_mark_model(file_path)

# Simulate a realization
example_mpp <- simulate_mpp(
  process = "self_correcting",
  process_fit = generating_parameters,
  t_min = 0,
  t_max = 1,
  anchor_point = M_n,
  raster_list = scaled_raster_list,
  scaled_rasters = TRUE,
  mark_model = mark_model,
  xy_bounds = c(0, 25, 0, 25),
  include_comp_inds = TRUE,
  competition_radius = 10,
  edge_correction = "none",
  thinning = TRUE,
  seed = 90210
)

# Plot the realization and provide a summary
plot(example_mpp, pattern_type = "simulated")

summary(example_mpp)
#> Summary: ldmppr Simulation
#>   process:          self_correcting
#>   n_points:         107
#>   mark_range:       [57.15, 639.738]
#>   time_range:       [0, 0.993054]
#>   thinning:         TRUE
#>   edge_correction:  none
#>   xy_bounds:        [0, 25, 0, 25]
```
