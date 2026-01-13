# Plot a marked point process

Plot a marked point process

## Usage

``` r
plot_mpp(mpp_data, pattern_type = c("reference", "simulated"))
```

## Arguments

- mpp_data:

  `ppp` object with marks or data frame with columns (x, y, size).

- pattern_type:

  type of pattern to plot ("reference" or "simulated").

## Value

a `ggplot` object of the marked point process.

## Examples

``` r
# Load example data
data(small_example_data)
mpp_data <- generate_mpp(
  locations = small_example_data %>% dplyr::select(x, y),
  marks = small_example_data$size,
  xy_bounds = c(0, 25, 0, 25)
)

# Plot the marked point process
plot_mpp(mpp_data, pattern_type = "reference")

```
