# Fitted point-process model object

Objects of class `ldmppr_fit` are returned by
[`estimate_process_parameters`](https://lanedrew.github.io/ldmppr/reference/estimate_process_parameters.md).
They contain the best-fitting optimization result (and optionally
multiple fits, e.g. from a delta search) along with metadata used to
reproduce the fit.

## Usage

``` r
# S3 method for class 'ldmppr_fit'
print(x, ...)

# S3 method for class 'ldmppr_fit'
coef(object, ...)

# S3 method for class 'ldmppr_fit'
logLik(object, ...)

# S3 method for class 'ldmppr_fit'
summary(object, ...)

# S3 method for class 'summary.ldmppr_fit'
print(x, ...)

# S3 method for class 'ldmppr_fit'
plot(x, ...)

as_nloptr(x, ...)

# S3 method for class 'ldmppr_fit'
as_nloptr(x, ...)
```

## Arguments

- x:

  an object of class `ldmppr_fit`.

- ...:

  additional arguments (not used).

- object:

  an object of class `ldmppr_fit`.

## Value

- [`print()`](https://rdrr.io/r/base/print.html):

  prints a brief summary of the fit.

- [`coef()`](https://rdrr.io/r/stats/coef.html):

  returns the estimated parameter vector.

- [`logLik()`](https://rdrr.io/r/stats/logLik.html):

  returns the log-likelihood at the optimum.

- [`summary()`](https://rdrr.io/r/base/summary.html):

  returns a `summary.ldmppr_fit`.

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html):

  plots diagnostics for multi-fit runs, if available.

## Details

A `ldmppr_fit` is a list with (at minimum):

- `process`: process name (e.g. `"self_correcting"`)

- `fit`: best optimization result (currently an `nloptr` object)

- `mapping`: mapping information (e.g. chosen `delta`, objectives)

- `grid`: grid definitions used by likelihood approximation

## Methods (by generic)

- `print(ldmppr_fit)`: Print a brief summary of a fitted model.

- `coef(ldmppr_fit)`: Extract the estimated parameter vector.

- `logLik(ldmppr_fit)`: Log-likelihood at the optimum.

- `summary(ldmppr_fit)`: Summarize a fitted model.

- `plot(ldmppr_fit)`: Plot diagnostics for a fitted model.

- `as_nloptr(ldmppr_fit)`: Extract the underlying `nloptr` result.

## Functions

- `print(summary.ldmppr_fit)`: Print a summary produced by
  `summary.ldmppr_fit`.

- `as_nloptr()`: Extract the underlying `nloptr` result.
