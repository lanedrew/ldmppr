# Predict values from the mark distribution

Predict values from the mark distribution

## Usage

``` r
predict_marks(
  sim_realization,
  raster_list = NULL,
  scaled_rasters = FALSE,
  mark_model = NULL,
  xy_bounds = NULL,
  include_comp_inds = FALSE,
  competition_radius = 15,
  edge_correction = "none"
)
```

## Arguments

- sim_realization:

  a data.frame containing a thinned or unthinned realization from
  `simulate_mpp` (or `simulate_sc`).

- raster_list:

  a list of raster objects.

- scaled_rasters:

  `TRUE` or `FALSE` indicating whether the rasters have been scaled.

- mark_model:

  a mark model object. May be a `ldmppr_mark_model` or a legacy model.

- xy_bounds:

  a vector of domain bounds (2 for x, 2 for y).

- include_comp_inds:

  `TRUE` or `FALSE` indicating whether to generate and use competition
  indices as covariates.

- competition_radius:

  distance for competition radius if `include_comp_inds` is `TRUE`.

- edge_correction:

  type of edge correction to apply (`"none"` or `"toroidal"`).

## Value

a vector of predicted mark values.

## Examples

``` r
# Simulate a realization
generating_parameters <- c(2, 8, .02, 2.5, 3, 1, 2.5, .2)
M_n <- c(10, 14)
generated_locs <- simulate_sc(
  t_min = 0,
  t_max = 1,
  sc_params = generating_parameters,
  anchor_point = M_n,
  xy_bounds = c(0, 25, 0, 25)
)

# Load the raster files
raster_paths <- list.files(system.file("extdata", package = "ldmppr"),
  pattern = "\\.tif$", full.names = TRUE
)
raster_paths <- raster_paths[!grepl("_med\\.tif$", raster_paths)]
rasters <- lapply(raster_paths, terra::rast)

# Scale the rasters
scaled_raster_list <- scale_rasters(rasters)

# Load the example mark model
file_path <- system.file("extdata", "example_mark_model.rds", package = "ldmppr")
mark_model <- load_mark_model(file_path)

# Predict the mark values
predict_marks(
  sim_realization = generated_locs$thinned,
  raster_list = scaled_raster_list,
  scaled_rasters = TRUE,
  mark_model = mark_model,
  xy_bounds = c(0, 25, 0, 25),
  include_comp_inds = TRUE,
  competition_radius = 10,
  edge_correction = "none"
)
#>   [1] 496.69830 550.63953 477.51013 512.50079 481.28769 493.59277 490.29263
#>   [8] 475.58575 520.82904 513.86536 466.98431 477.79315 517.88696 450.52911
#>  [15] 526.16461 470.89828 517.12195 471.04663 474.07663 510.92596 486.31201
#>  [22] 430.03122 499.84094 495.75162 480.13022 471.32767 469.11176 478.07281
#>  [29] 469.30585 474.80551 511.08969 477.87561 465.47995 470.30972 528.49323
#>  [36] 498.46219 477.69534 468.89758 465.75659 517.54169 471.07852 480.81259
#>  [43] 465.67297 461.38516 479.73038 478.19086 467.98520 377.38306 335.59259
#>  [50] 337.80270 356.79865 326.40955 333.52328 326.31406 352.92627 326.24289
#>  [57] 342.62057 332.72867 332.19321 314.01047 330.65536 300.84155 329.46799
#>  [64] 325.08774 336.03821 323.76074 332.38824 312.55121 256.41049 285.32153
#>  [71] 266.51505 270.69623 273.53265 269.66666 262.36526 255.72971 266.85660
#>  [78] 233.50487 242.79892 238.11795 211.03502 222.13040 208.85739 213.08667
#>  [85] 195.47977 204.04666 191.00299 165.93750 166.26091 170.26016 172.70920
#>  [92] 182.46190 165.26004 171.35146 161.47290 156.36734 152.04909 147.67635
#>  [99] 130.55685 122.88330 124.05809 125.68738 126.46715 106.92394  84.41629
#> [106]  79.76232
```
