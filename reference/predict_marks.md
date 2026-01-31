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
#>   [1] 478.52460 485.39972 477.22626 523.52246 526.79694 505.02661 512.42023
#>   [8] 478.67145 517.76581 507.37003 441.06317 428.67139 464.64606 451.84088
#>  [15] 524.59839 505.44232 469.88852 469.70685 477.94470 477.92010 477.77658
#>  [22] 477.46381 477.84921 526.60779 476.07861 470.56546 477.85986 511.05966
#>  [29] 514.10797 506.83389 470.73532 508.83652 470.57806 459.89880 468.65213
#>  [36] 471.67233 469.29099 470.70377 469.93533 471.21109 469.14642 424.63986
#>  [43] 478.52399 463.88010 463.77939 466.35175 452.53641 450.49976 470.75027
#>  [50] 517.07892 460.29861 507.98456 336.72302 307.25616 366.44418 326.68359
#>  [57] 334.71240 331.42584 338.46588 368.20856 335.86038 364.17630 324.89761
#>  [64] 324.88977 310.96869 296.59030 320.27115 278.25345 275.19797 285.45914
#>  [71] 268.44699 277.08783 278.86362 276.11969 255.22427 243.73689 244.92024
#>  [78] 226.26805 239.85568 220.88206 227.04111 240.29170 205.67236 208.85669
#>  [85] 215.72606 216.32140 212.16431 204.49498 203.00606 196.46677 192.83017
#>  [92] 178.08577 179.86328 168.03270 177.74602 184.85240 133.72577 109.57582
#>  [99] 119.09148  81.94879  79.82088
```
