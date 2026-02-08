# Create a grid schedule for estimate_process_parameters()

`ldmppr_grids()` constructs a multi-resolution grid schedule used by
[`estimate_process_parameters`](https://lanedrew.github.io/ldmppr/reference/estimate_process_parameters.md).
The returned object is an S3 class with helper methods; see
[`ldmppr_grids-class`](https://lanedrew.github.io/ldmppr/reference/ldmppr_grids-class.md).

## Usage

``` r
ldmppr_grids(upper_bounds, levels, labels = NULL, include_endpoints = TRUE)
```

## Arguments

- upper_bounds:

  a vector `c(b_t, b_x, b_y)` giving the maximum bounds for time and the
  spatial domain. Grids must lie within these.

- levels:

  a list describing the grid schedule. Each entry can be either:

  - a numeric length-3 vector `c(nx, ny, nt)` (number of points per
    dimension),

  - a list with elements `nx`, `ny`, `nt`,

  - a list with explicit vectors `x`, `y`, `t`.

- labels:

  (optional) character vector of length equal to `levels`, used only for
  printing.

- include_endpoints:

  `TRUE` or `FALSE` indicating; if `TRUE` (default) each generated grid
  uses `seq(0, bound, length.out = n)` including endpoints.

## Value

an object of class `"ldmppr_grids"`.

## See also

[`ldmppr_grids-class`](https://lanedrew.github.io/ldmppr/reference/ldmppr_grids-class.md)
for methods and details.

## Examples

``` r
# A 3-level coarse-to-fine schedule (counts per dimension)
g <- ldmppr_grids(
  upper_bounds = c(1, 50, 50),
  levels = list(
    c(25, 25, 25),
    c(75, 75, 75),
    c(100, 100, 100)
  )
)
g
#> <ldmppr_grids>
#>   upper_bounds: b_t=1, b_x=50, b_y=50
#>   levels: 3
#>     - level 1: 25x25x25  [x:0..50, y:0..50, t:0..1]
#>     - level 2: 75x75x75  [x:0..50, y:0..50, t:0..1]
#>     - level 3: 100x100x100  [x:0..50, y:0..50, t:0..1]
length(g)
#> [1] 3
summary(g)
#> <summary: ldmppr_grids>
#>   levels: 3
#>   dims:
#>       nx  ny  nt
#> [1,]  25  25  25
#> [2,]  75  75  75
#> [3,] 100 100 100

# Explicit vectors (single level)
g2 <- ldmppr_grids(
  upper_bounds = c(1, 50, 50),
  levels = list(list(
    x = seq(0, 50, by = 2),
    y = seq(0, 50, by = 2),
    t = seq(0, 1,  length.out = 30)
  ))
)
as.data.frame(g2)
#>   level label nx ny nt x_min x_max y_min y_max t_min t_max
#> 1     1       26 26 30     0    50     0    50     0     1
```
