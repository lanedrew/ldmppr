# Extract covariate values from a set of rasters

Extract covariate values from a set of rasters

## Usage

``` r
extract_covars(locations, raster_list)
```

## Arguments

- locations:

  a matrix/data.frame of (x,y) locations.

- raster_list:

  a list of SpatRaster objects.

## Value

a data.frame of covariates (no ID column; unique names).

## Examples

``` r
# Load example raster data
raster_paths <- list.files(system.file("extdata", package = "ldmppr"),
  pattern = "\\.tif$", full.names = TRUE
)
raster_paths <- raster_paths[!grepl("_med\\.tif$", raster_paths)]
rasters <- lapply(raster_paths, terra::rast)

# Scale the rasters
scaled_raster_list <- scale_rasters(rasters)

# Load example locations
locations <- small_example_data %>%
  dplyr::select(x, y) %>%
  as.matrix()

# Extract covariates
example_covars <- extract_covars(locations, scaled_raster_list)
head(example_covars)
#>   Snodgrass_DEM_1m Snodgrass_aspect_southness_1m Snodgrass_slope_1m
#> 1       -0.3878719                    -0.2481063          0.4541215
#> 2       -1.3977397                     1.7131221          0.6044641
#> 3        1.6946886                    -1.1153459         -1.5861887
#> 4        0.7115911                    -0.9918452         -1.9867798
#> 5        1.4753020                    -0.5887282         -1.3989392
#> 6        1.5844332                     0.7520857         -2.1610400
#>   Snodgrass_wetness_index_1m
#> 1                 -0.4099315
#> 2                 -0.2639008
#> 3                  1.2917963
#> 4                  2.0485587
#> 5                  0.9809867
#> 6                  1.8836371
```
