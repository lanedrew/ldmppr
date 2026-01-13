# Simulated marked point process object

`ldmppr_sim` objects are returned by
[`simulate_mpp`](https://lanedrew.github.io/ldmppr/reference/simulate_mpp.md).
They contain the simulated realization, an associated marked point
pattern object, and metadata used to reproduce or inspect the
simulation.

## Usage

``` r
# S3 method for class 'ldmppr_sim'
print(x, ...)

# S3 method for class 'ldmppr_sim'
as.data.frame(x, ...)

# S3 method for class 'ldmppr_sim'
nobs(object, ...)

# S3 method for class 'ldmppr_sim'
plot(x, pattern_type = "simulated", ...)

mpp.ldmppr_sim(x, ...)
```

## Arguments

- x:

  a `ldmppr_sim` object.

- ...:

  additional arguments (not used).

- object:

  a `ldmppr_sim` object.

- pattern_type:

  type of pattern to plot `"simulated"` (default).

## Value

For methods:

- [`print()`](https://rdrr.io/r/base/print.html):

  prints a summary of the simulation.

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html):

  returns a ggplot visualization of the marked point pattern.

- [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html):

  returns the simulated realization as a data.frame.

- [`nobs()`](https://rdrr.io/r/stats/nobs.html):

  returns the number of points in the realization.

- `mpp()`:

  returns the marked point pattern object.

## Details

An `ldmppr_sim` is a list with at least:

- `process`: process name (e.g. `"self_correcting"`)

- `mpp`: a marked point pattern object

- `realization`: data.frame with columns `time`, `x`, `y`, `marks`

- `params`, `bounds`, and other metadata

## Methods (by generic)

- `print(ldmppr_sim)`: Print a brief summary of the simulation.

- `as.data.frame(ldmppr_sim)`: Coerce to a data.frame of the simulated
  realization.

- `nobs(ldmppr_sim)`: Number of simulated points.

- `plot(ldmppr_sim)`: Plot the simulated marked point pattern.

## Functions

- `mpp.ldmppr_sim()`: Extract the underlying marked point pattern
  object.
