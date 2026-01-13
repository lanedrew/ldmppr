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
#>   [1] 518.59125 479.85587 470.27548 471.02344 478.63190 477.08408 497.45270
#>   [8] 470.51154 478.63852 511.44165 558.97974 471.11227 472.54587 471.25854
#>  [15] 482.47470 479.34631 475.19351 444.87729 470.11722 471.01309 477.76242
#>  [22] 452.61002 469.33609 478.27426 478.08994 505.18069 469.88419 478.36362
#>  [29] 479.59113 548.90814 553.19800 466.33533 471.19601 495.80777 477.65985
#>  [36] 479.37799 475.51166 468.64081 471.29211 477.84937 509.66626 490.92535
#>  [43] 465.11731 476.91147 478.11774 434.09064 478.56799 468.86932 433.20410
#>  [50] 350.38388 358.95496 370.80170 354.98636 335.49127 292.93958 319.67065
#>  [57] 338.43420 332.11450 308.02448 322.24384 331.16974 333.44318 317.78787
#>  [64] 317.22314 322.52277 329.91147 330.51523 334.02448 334.94858 333.43835
#>  [71] 322.06406 310.17987 276.25070 276.38776 275.34729 262.32510 257.54297
#>  [78] 247.45091 247.96873 223.28220 244.67004 212.57436 227.73199 225.50082
#>  [85] 223.78554 242.52878 202.63051 224.38615 226.47736 216.41750 206.86746
#>  [92] 207.76668 188.12369 199.71930 182.47678 170.85396 169.01073 160.15585
#>  [99] 182.87686 176.78508 175.67845 141.68456 135.17970 139.78154 134.10811
#> [106] 127.12403 128.58470 119.35510  98.13485  87.61633  55.28815
```
