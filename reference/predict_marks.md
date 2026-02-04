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
#>   [1] 478.60419 550.63953 477.90698 512.50079 481.19864 515.62238 447.87650
#>   [8] 475.72870 477.95331 490.79068 453.07886 477.52185 478.09671 449.84671
#>  [15] 535.62219 470.79535 516.90540 471.04663 536.38623 510.93567 465.46246
#>  [22] 430.03122 511.12341 477.66531 479.64542 511.17557 478.71515 469.14319
#>  [29] 477.75848 478.26489 474.80551 511.10114 477.87930 466.13617 479.68292
#>  [36] 470.32117 528.49323 498.46219 478.57635 469.03687 583.84149 495.85376
#>  [43] 471.07852 480.81259 466.17786 470.69772 518.74060 377.36536 377.71265
#>  [50] 337.69907 337.10513 333.81494 355.31027 327.92233 323.83075 367.51651
#>  [57] 328.76309 323.25378 357.47913 314.33759 322.86954 300.86902 329.46799
#>  [64] 325.82834 348.16464 320.03094 324.61841 323.39676 332.46262 288.25809
#>  [71] 256.05487 281.81244 263.67923 304.34985 249.28656 255.37373 266.99313
#>  [78] 255.43820 240.76312 233.22511 214.93825 200.09312 222.13016 208.15224
#>  [85] 208.34506 182.32797 220.22682 185.31569 183.82967 165.50479 167.02931
#>  [92] 173.80612 174.59676 174.06190 160.64307 156.12914 163.40846 147.63071
#>  [99] 116.98822 136.99391 122.88330 124.28596 105.56876  82.29893  80.81445
#> [106]  73.55502  72.50830  49.74652
```
