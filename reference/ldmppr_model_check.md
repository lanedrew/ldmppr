# Model fit diagnostic object

Objects of class `ldmppr_model_check` are returned by
[`check_model_fit`](https://lanedrew.github.io/ldmppr/reference/check_model_fit.md).
They contain global envelope test results and curve sets for multiple
summary functions/statistics.

## Usage

``` r
# S3 method for class 'ldmppr_model_check'
print(x, ...)

# S3 method for class 'ldmppr_model_check'
summary(object, ...)

# S3 method for class 'summary.ldmppr_model_check'
print(x, ...)

# S3 method for class 'ldmppr_model_check'
plot(x, which = c("combined", "L", "F", "G", "J", "E", "V"), ...)
```

## Arguments

- x:

  an object of class `ldmppr_model_check`.

- ...:

  additional arguments passed to the underlying
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method (e.g.,
  from \*\*GET\*\*).

- object:

  an object of class `ldmppr_model_check`.

- which:

  which envelope to plot. `"combined"` plots the global envelope;
  otherwise one of `"L"`, `"F"`, `"G"`, `"J"`, `"E"`, `"V"`.

## Value

- [`print()`](https://rdrr.io/r/base/print.html):

  prints a brief summary of the diagnostic object.

- [`summary()`](https://rdrr.io/r/base/summary.html):

  returns a `summary.ldmppr_model_check` object.

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html):

  plots the combined envelope or a selected statistic envelope.

## Details

An `ldmppr_model_check` is a list with components such as:

- `combined_env`: a global envelope test object (typically from
  \*\*GET\*\*)

- `envs`: named list of envelope test objects (e.g., `L`, `F`, `G`, `J`,
  `E`, `V`)

- `curve_sets`: named list of curve set objects

- `settings`: list of settings used when generating envelopes (e.g.,
  `n_sim`, `thinning`)

## Methods (by generic)

- `print(ldmppr_model_check)`: Print a brief summary of the diagnostic
  results.

- `summary(ldmppr_model_check)`: Summarize p-values from the combined
  and individual envelopes.

- `plot(ldmppr_model_check)`: Plot the combined envelope or a selected
  statistic.

## Functions

- `print(summary.ldmppr_model_check)`: Print a summary produced by
  `summary.ldmppr_model_check`.
