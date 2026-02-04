# Scale a set of rasters

Scale a set of rasters

## Usage

``` r
scale_rasters(raster_list, reference_resolution = NULL)
```

## Arguments

- raster_list:

  a list of raster objects.

- reference_resolution:

  the resolution to resample the rasters to.

## Value

a list of scaled raster objects.

## Examples

``` r
# Create two example rasters
rast_a <- terra::rast(
  ncol = 10, nrow = 10,
  xmin = 0, xmax = 10,
  ymin = 0, ymax = 10,
  vals = runif(100)
)

rast_b <- terra::rast(
  ncol = 10, nrow = 10,
  xmin = 0, xmax = 10,
  ymin = 0, ymax = 10,
  vals = runif(100)
)

# Scale example rasters in a list
rast_list <- list(rast_a, rast_b)
scale_rasters(rast_list)
#> [[1]]
#> class       : SpatRaster 
#> size        : 10, 10, 1  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 10, 0, 10  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84) 
#> source(s)   : memory
#> name        :     lyr.1 
#> min value   : -2.002002 
#> max value   :  1.479849 
#> 
#> [[2]]
#> class       : SpatRaster 
#> size        : 10, 10, 1  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 10, 0, 10  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84) 
#> source(s)   : memory
#> name        :     lyr.1 
#> min value   : -1.540647 
#> max value   :  1.840785 
#> 
```
