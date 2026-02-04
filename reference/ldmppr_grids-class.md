# Grid schedule object

Objects of class `ldmppr_grids` define one or more grid "levels" used by
[`estimate_process_parameters`](https://lanedrew.github.io/ldmppr/reference/estimate_process_parameters.md).
Each level contains numeric vectors `x`, `y`, and `t` defining the
approximation grid. Levels are typically ordered from coarse to fine.

## Usage

``` r
# S3 method for class 'ldmppr_grids'
print(x, ...)

# S3 method for class 'ldmppr_grids'
summary(object, ...)

# S3 method for class 'summary.ldmppr_grids'
print(x, ...)

# S3 method for class 'ldmppr_grids'
as.data.frame(x, ...)

# S3 method for class 'ldmppr_grids'
length(x)

# S3 method for class 'ldmppr_grids'
x[i, ...]

# S3 method for class 'ldmppr_grids'
as.list(x, ...)
```

## Arguments

- x:

  an object of class `ldmppr_grids`.

- ...:

  unused.

- object:

  an object of class `ldmppr_grids`.

- i:

  indices of levels to keep.

## Value

- [`print()`](https://rdrr.io/r/base/print.html):

  prints a brief description of bounds and grid levels.

- [`summary()`](https://rdrr.io/r/base/summary.html):

  returns a `summary.ldmppr_grids`.

- [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html):

  returns one row per level with dimensions and ranges.

- [`length()`](https://rdrr.io/r/base/length.html):

  returns the number of levels.

- `[ ]`:

  subsets levels, preserving class.

- [`as.list()`](https://rdrr.io/r/base/list.html):

  returns the underlying list structure.

## Details

A `ldmppr_grids` is a list with (at minimum):

- `levels`: list of levels; each level is a list with `x`, `y`, `t`

- `upper_bounds`: numeric `c(b_t, b_x, b_y)`

- `labels`: optional labels used only for printing

- `include_endpoints`: logical

## Functions

- `print(ldmppr_grids)`: Print a brief summary of a grid schedule.

- `summary(ldmppr_grids)`: Summarize a grid schedule.

- `print(summary.ldmppr_grids)`: Print a summary produced by
  `summary.ldmppr_grids()`.

- `as.data.frame(ldmppr_grids)`: Convert a grid schedule to a
  data.frame.

- `length(ldmppr_grids)`: Number of levels in a grid schedule.

- `[`: Subset grid levels.

- `as.list(ldmppr_grids)`: Extract the underlying list representation.
