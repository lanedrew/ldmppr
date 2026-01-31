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
  info = list()
)

# S3 method for class 'ldmppr_mark_model'
print(x, ...)

# S3 method for class 'ldmppr_mark_model'
predict(object, new_data, ...)

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

- info:

  (optional) list of metadata.

- x:

  a `ldmppr_mark_model` object.

- ...:

  passed to methods.

- object:

  a `ldmppr_mark_model` object.

- new_data:

  a data frame of predictors (and possibly outcome columns).

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

- `predict(ldmppr_mark_model)`: Predict marks for new data.

- `save_mark_model(ldmppr_mark_model)`: Save method for
  `ldmppr_mark_model`.

## Functions

- `ldmppr_mark_model()`: Create a mark model container.

- `save_mark_model()`: Save a mark model to disk.

- `load_mark_model()`: Load a saved mark model from disk.
