# Mark model object

`ldmppr_mark_model` objects store a fitted mark model and preprocessing
information used to predict marks at new locations and times. These
objects are typically returned by
[`train_mark_model`](https://lanedrew.github.io/ldmppr/reference/train_mark_model.md)
and can be saved/loaded with `save_mark_model` and `load_mark_model`.

## Usage

``` r
ldmppr_mark_model(
  engine,
  fit_engine = NULL,
  xgb_raw = NULL,
  recipe = NULL,
  outcome = "size",
  feature_names = NULL,
  rasters = NULL,
  info = list()
)

# S3 method for class 'ldmppr_mark_model'
print(x, ...)

# S3 method for class 'ldmppr_mark_model'
summary(object, ...)

# S3 method for class 'summary.ldmppr_mark_model'
print(x, ...)

# S3 method for class 'ldmppr_mark_model'
predict(
  object,
  new_data = NULL,
  sim_realization = NULL,
  raster_list = NULL,
  scaled_rasters = FALSE,
  xy_bounds = NULL,
  include_comp_inds = FALSE,
  competition_radius = 15,
  edge_correction = "none",
  seed = NULL,
  ...
)

save_mark_model(object, path, ...)

# S3 method for class 'ldmppr_mark_model'
save_mark_model(object, path, ...)

load_mark_model(path)
```

## Arguments

- engine:

  character string (currently `"xgboost"` and `"ranger"`).

- fit_engine:

  fitted engine object (e.g. `xgb.Booster` or a ranger fit).

- xgb_raw:

  raw xgboost payload (e.g. UBJ) used for rehydration.

- recipe:

  a prepped recipes object used for preprocessing new data.

- outcome:

  outcome column name (default `"size"`).

- feature_names:

  (optional) vector of predictor names required at prediction time.

- rasters:

  (optional) list of rasters used for prediction (e.g. for spatial
  covariates).

- info:

  (optional) list of metadata.

- x:

  an object of class `summary.ldmppr_mark_model`.

- ...:

  passed to methods.

- object:

  a `ldmppr_mark_model` object.

- new_data:

  a data frame of predictors (and possibly outcome columns). Ignored
  when `sim_realization` is supplied.

- sim_realization:

  optional simulation realization containing `x`, `y`, and `time`. When
  supplied, predictors are built from rasters and optional competition
  indices.

- raster_list:

  optional list of rasters used when `sim_realization` is supplied. If
  omitted, uses rasters stored in `object` when available.

- scaled_rasters:

  `TRUE` or `FALSE`; whether supplied rasters are pre-scaled.

- xy_bounds:

  domain bounds `c(a_x, b_x, a_y, b_y)` used for competition indices.

- include_comp_inds:

  `TRUE` or `FALSE`; include competition-index features.

- competition_radius:

  positive numeric distance used when `include_comp_inds = TRUE`.

- edge_correction:

  edge correction for competition indices (`"none"` or `"toroidal"`).

- seed:

  optional nonnegative integer seed.

- path:

  path to an `.rds` created by `save_mark_model` (or legacy objects).

## Value

- [`print()`](https://rdrr.io/r/base/print.html):

  prints a brief summary.

- [`predict()`](https://rdrr.io/r/stats/predict.html):

  returns numeric predictions for new data.

an object of class `"ldmppr_mark_model"`.

## Details

The model may be backed by different engines (currently `"xgboost"` and
`"ranger"`). For `"xgboost"`, the object can store a serialized booster
payload to make saving/loading robust across R sessions.

## Methods (by generic)

- `print(ldmppr_mark_model)`: Print a brief summary of the mark model.

- `summary(ldmppr_mark_model)`: Summarize a mark model.

- `predict(ldmppr_mark_model)`: Predict marks for new data.

- `save_mark_model(ldmppr_mark_model)`: Save method for
  `ldmppr_mark_model`.

## Functions

- `ldmppr_mark_model()`: Create a mark model container.

- `print(summary.ldmppr_mark_model)`: Print a summary produced by
  `summary.ldmppr_mark_model`.

- `save_mark_model()`: Save a mark model to disk.

- `load_mark_model()`: Load a saved mark model from disk.
