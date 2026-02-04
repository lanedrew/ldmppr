# Optimization budget specification object

Objects of class `ldmppr_budgets` define optimizer budget/options used
by
[`estimate_process_parameters`](https://lanedrew.github.io/ldmppr/reference/estimate_process_parameters.md).

## Usage

``` r
# S3 method for class 'ldmppr_budgets'
print(x, ...)

# S3 method for class 'ldmppr_budgets'
summary(object, ...)

# S3 method for class 'summary.ldmppr_budgets'
print(x, ...)

# S3 method for class 'ldmppr_budgets'
as.data.frame(x, ...)

# S3 method for class 'ldmppr_budgets'
length(x)

# S3 method for class 'ldmppr_budgets'
x[i, ...]

# S3 method for class 'ldmppr_budgets'
as.list(x, ...)
```

## Arguments

- x:

  an object of class `ldmppr_budgets`.

- ...:

  unused.

- object:

  an object of class `ldmppr_budgets`.

- i:

  indices of local stages to keep: 1 = first level, 2 = refinement
  levels.

## Value

- [`print()`](https://rdrr.io/r/base/print.html):

  prints a brief description of configured budgets.

- [`summary()`](https://rdrr.io/r/base/summary.html):

  returns a `summary.ldmppr_budgets`.

- [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html):

  a compact table of the global/local budget entries (when present).

- [`length()`](https://rdrr.io/r/base/length.html):

  number of available local budget stages (1 or 2).

- `[ ]`:

  subset local budget stages (keeps global options).

- [`as.list()`](https://rdrr.io/r/base/list.html):

  returns the underlying list structure.

## Details

A `ldmppr_budgets` is a list with (at minimum):

- `global_options`: list of NLopt options for the global stage (e.g.,
  `maxeval`, `maxtime`).

- `local_budget_first_level`: list of NLopt options for the local stage
  at the first/coarsest grid level.

- `local_budget_refinement_levels`: optional list of NLopt options for
  local refinement levels (used only when
  `estimate_process_parameters(strategy = "multires_global_local")`).

## Functions

- `print(ldmppr_budgets)`: Print a brief summary of optimization
  budgets.

- `summary(ldmppr_budgets)`: Summarize an optimization budget
  specification.

- `print(summary.ldmppr_budgets)`: Print a summary produced by
  `summary.ldmppr_budgets()`.

- `as.data.frame(ldmppr_budgets)`: Convert budgets to a data.frame.

- `length(ldmppr_budgets)`: Number of available local budget stages.

- `[`: Subset local budget stages (keeps global options).

- `as.list(ldmppr_budgets)`: Extract the underlying list representation.
