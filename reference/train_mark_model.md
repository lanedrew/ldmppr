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
  parallel = TRUE,
  num_cores = NULL,
  include_comp_inds = FALSE,
  competition_radius = 15,
  edge_correction = "none",
  selection_metric = "rmse",
  cv_folds = 5,
  tuning_grid_size = 200,
  verbose = TRUE
)
```

## Arguments

- data:

  a data.frame or a `ldmppr_fit` object. See Description.

- raster_list:

  a list of raster objects.

- scaled_rasters:

  `TRUE` or `FALSE` indicating whether the rasters have been scaled.

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

  `TRUE` or `FALSE` indicating whether to use parallelization in model
  training.

- num_cores:

  number of cores to use in parallel model training (if `parallel` is
  `TRUE`).

- include_comp_inds:

  `TRUE` or `FALSE` indicating whether to generate and use competition
  indices as covariates.

- competition_radius:

  distance for competition radius if `include_comp_inds` is `TRUE`.

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
#> Training mark model
#>   Model type: xgboost
#>   Selection metric: rmse
#>   CV folds: 3
#>   Tuning grid size: 2
#>   Include competition indices: no
#>   Edge correction: none
#> Step 1/6: Preparing training data...
#>   Rows: 121
#>   Done in 0.0s.
#> Step 2/6: Configuring parallel backend...
#>   Parallel: off
#>   Model engine threads: 3
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
#>   Tuning enabled: 3-fold CV with 2 candidate(s).
#>   Total model fits: 6 (3 folds x 2 grid).
#>   Done in 2.8s.
#> Step 6/6: Finalizing output object...
#>   Done in 0.0s.
#> Training complete. Total time: 2.8s.

print(mark_model)
#> <ldmppr_mark_model>
#>   engine: xgboost
#>   has fit_engine: TRUE
#>   has xgb_raw: FALSE
#>   n_features: 7
```
