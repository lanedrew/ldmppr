# Optimized function to compute toroidal distance matrix over a rectangular domain

Optimized function to compute toroidal distance matrix over a
rectangular domain

## Usage

``` r
toroidal_dist_matrix_optimized(location_matrix, x_bound, y_bound)
```

## Arguments

- location_matrix:

  a 2 column matrix of (x,y) coordinates.

- x_bound:

  the upper bound for the x dimension.

- y_bound:

  the upper bound for the y dimension.

## Value

a matrix of toroidal distances.

## Examples

``` r
# Generate a matrix of locations
location_matrix <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 2)
x_bound <- 10
y_bound <- 10

# Compute the toroidal distance matrix
toroidal_dist_matrix_optimized(location_matrix, x_bound, y_bound)
#>          [,1]     [,2]     [,3]
#> [1,] 0.000000 1.414214 2.828427
#> [2,] 1.414214 0.000000 1.414214
#> [3,] 2.828427 1.414214 0.000000
```
