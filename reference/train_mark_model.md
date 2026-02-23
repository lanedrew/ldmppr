# Train a flexible model for the mark distribution

Trains a predictive model for the mark distribution of a spatio-temporal
process. `data` may be either (1) a data.frame containing columns `x`,
`y`, `size` and `time`, (2) a data.frame containing `x`, `y`, `size`
(time will be derived via `delta`), or (3) a `ldmppr_fit` object
returned by
[`estimate_process_parameters`](https://lanedrew.github.io/ldmppr/reference/estimate_process_parameters.md).
Allows the user to incorporate location specific information and
competition indices as covariates in the mark model.

## Usage

``` r
train_mark_model(
  data,
  raster_list = NULL,
  scaled_rasters = FALSE,
  model_type = "xgboost",
  xy_bounds = NULL,
  delta = NULL,
  save_model = FALSE,
  save_path = NULL,
  parallel = FALSE,
  num_cores = NULL,
  include_comp_inds = FALSE,
  competition_radius = 15,
  edge_correction = "none",
  selection_metric = "rmse",
  cv_folds = 5,
  tuning_grid_size = 200,
  seed = 0,
  verbose = TRUE
)
```

## Arguments

- data:

  a data.frame or a `ldmppr_fit` object. See Description.

- raster_list:

  list of raster objects used for mark-model training.

- scaled_rasters:

  `TRUE` or `FALSE` indicating whether rasters are already scaled.

- model_type:

  the machine learning model type (`"xgboost"` or `"random_forest"`).

- xy_bounds:

  a vector of domain bounds (2 for x, 2 for y). If `data` is an
  `ldmppr_fit` and `xy_bounds` is `NULL`, defaults to
  `c(0, b_x, 0, b_y)` derived from fit.

- delta:

  (optional) numeric scalar used only when `data` contains `(x,y,size)`
  but not `time`. If `data` is an `ldmppr_fit` and time is missing, the
  function will infer the `delta` value from the fit.

- save_model:

  `TRUE` or `FALSE` indicating whether to save the generated model.

- save_path:

  path for saving the generated model.

- parallel:

  `TRUE` or `FALSE`. If `TRUE`, tuning is parallelized over resamples.
  For small datasets, parallel overhead may outweigh speed gains.

- num_cores:

  number of workers to use when `parallel=TRUE`. Ignored when
  `parallel=FALSE`.

- include_comp_inds:

  `TRUE` or `FALSE` indicating whether to generate and use competition
  indices as covariates.

- competition_radius:

  positive numeric distance used when `include_comp_inds = TRUE`.

- edge_correction:

  type of edge correction to apply (`"none"`, `"toroidal"`, or
  `"truncation"`).

- selection_metric:

  metric to use for identifying the optimal model (`"rmse"`, `"mae"`, or
  `"rsq"`).

- cv_folds:

  number of cross-validation folds to use in model training. If
  `cv_folds <= 1`, tuning is skipped and the model is fit once with
  default hyperparameters.

- tuning_grid_size:

  size of the tuning grid for hyperparameter tuning.

- seed:

  integer seed for reproducible resampling/tuning/model fitting.

- verbose:

  `TRUE` or `FALSE` indicating whether to show progress of model
  training.

## Value

an object of class `"ldmppr_mark_model"` containing the trained mark
model.

## Examples

``` r
# Load the small example data
data(small_example_data)

# Load example raster data
raster_paths <- list.files(system.file("extdata", package = "ldmppr"),
  pattern = "\\.tif$", full.names = TRUE
)
raster_paths <- raster_paths[!grepl("_med\\.tif$", raster_paths)]
rasters <- lapply(raster_paths, terra::rast)

# Scale the rasters
scaled_raster_list <- scale_rasters(rasters)


# Train the model
mark_model <- train_mark_model(
  data = small_example_data,
  raster_list = scaled_raster_list,
  scaled_rasters = TRUE,
  model_type = "xgboost",
  xy_bounds = c(0, 25, 0, 25),
  delta = 1,
  parallel = FALSE,
  include_comp_inds = FALSE,
  competition_radius = 10,
  edge_correction = "none",
  selection_metric = "rmse",
  cv_folds = 3,
  tuning_grid_size = 2,
  verbose = TRUE
)
#> [ldmppr::train_mark_model]
#> Training mark model
#>   Model type: xgboost
#>   Selection metric: rmse
#>   CV folds: 3
#>   Tuning grid size: 2
#>   Include competition indices: no
#>   Edge correction: none
#> Step 1/6: Preparing training data...
#>   Rows: 121
#>   Seed: 0
#>   Done in 0.0s.
#> Step 2/6: Configuring parallel backend...
#>   Parallel: off
#>   Model engine threads: 1
#>   Done in 0.0s.
#> Step 3/6: Extracting raster covariates...
#>   Using pre-scaled rasters (scaled_rasters = TRUE).
#>   Extracted 4 raster feature(s).
#>   Done in 0.0s.
#> Step 4/6: Building model matrix (and competition indices if requested)...
#>   Final training rows: 121
#>   Final feature columns (incl x,y,time): 7
#>   Done in 0.0s.
#> Step 5/6: Fitting model (with optional CV tuning)...
#>   foreach backend: doSEQ | workers=1
#>   Done in 2.5s.
#> Step 6/6: Finalizing output object...
#>   Residual bootstrap stored (source=oos, transform=sqrt, bins=6, min/bin=8).
#>   Done in 0.1s.
#> Training complete. Total time: 2.6s.

print(mark_model)
#> ldmppr Mark Model
#>   engine:           xgboost
#>   has_fit_engine:   TRUE
#>   has_xgb_raw:      FALSE
#>   n_features:       7
#>   n_rasters:        4
#>   raster_names:     Snodgrass_DEM_1m, Snodgrass_aspect_southness_1m,
#>                     Snodgrass_slope_1m, Snodgrass_wetness_index_1m
#>   scaled_rasters:   TRUE
#>   comp_indices:     FALSE
```
