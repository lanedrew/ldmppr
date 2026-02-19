# Predict values from the mark distribution

Prefer using the S3 method
[`predict()`](https://rdrr.io/r/stats/predict.html) on an
`ldmppr_mark_model`:
`predict(mark_model, sim_realization = ..., xy_bounds = ...)`. This
wrapper is retained for backward compatibility and is deprecated.

## Usage

``` r
predict_marks(
  sim_realization,
  raster_list = NULL,
  scaled_rasters = FALSE,
  mark_model = NULL,
  xy_bounds = NULL,
  include_comp_inds = FALSE,
  competition_radius = 15,
  edge_correction = "none",
  seed = NULL
)
```

## Arguments

- sim_realization:

  a data.frame containing a thinned or unthinned realization from
  `simulate_mpp` (or `simulate_sc`). Must contain `x`, `y`, and `time`.

- raster_list:

  list of raster objects used for mark prediction.

- scaled_rasters:

  `TRUE` or `FALSE` indicating whether rasters are already scaled.

- mark_model:

  a mark model object. May be an `ldmppr_mark_model`, `model_fit`, or
  `workflow`.

- xy_bounds:

  vector of domain bounds as `c(a_x, b_x, a_y, b_y)`.

- include_comp_inds:

  `TRUE` or `FALSE` indicating whether to generate and use competition
  indices as covariates.

- competition_radius:

  positive numeric distance used when `include_comp_inds = TRUE`.

- edge_correction:

  type of edge correction to apply (`"none"` or `"toroidal"`).

- seed:

  optional nonnegative integer seed for reproducibility.

## Value

a vector of predicted mark values.
