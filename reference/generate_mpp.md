# Generate a marked process given locations and marks

Creates an object of class "ppp" that represents a marked point pattern
in the two-dimensional plane.

## Usage

``` r
generate_mpp(locations, marks = NULL, xy_bounds = NULL)
```

## Arguments

- locations:

  a data.frame of (x,y) locations with names "x" and "y".

- marks:

  a vector of marks.

- xy_bounds:

  a vector of domain bounds (2 for x, 2 for y).

## Value

a ppp object with marks.

## Examples

``` r
# Load example data
data(small_example_data)

# Generate a marked point process
generate_mpp(
  locations = small_example_data %>% dplyr::select(x, y),
  marks = small_example_data$size,
  xy_bounds = c(0, 25, 0, 25)
)
#> Marked planar point pattern: 121 points
#> marks are numeric, of storage type  ‘double’
#> window: rectangle = [0, 25] x [0, 25] units
```
