# This code generates four raster objects from the Snodgrass dataset located at https://data.ess-dive.lbl.gov/view/doi:10.15485/2476543.
# The rasters are cropped to a smaller extent to use in examples and testing.

# Specify the spatial extent
a_x <- 326996
a_y <- 4311239
b_x <- 327021
b_y <- 4311264

new_ext <- c(a_x, b_x, a_y, b_y)

# Load the raster images
south <- terra::rast("https://data.ess-dive.lbl.gov/catalog/d1/mn/v2/object/ess-dive-c3c46ff25d50885-20240513T173925432")
wet <- terra::rast("https://data.ess-dive.lbl.gov/catalog/d1/mn/v2/object/ess-dive-dffdeec81023d23-20240513T173925427")
slope <- terra::rast("https://data.ess-dive.lbl.gov/catalog/d1/mn/v2/object/ess-dive-8a59960f4ffd550-20240513T173925429")
DEM <- terra::rast("https://data.ess-dive.lbl.gov/catalog/d1/mn/v2/object/ess-dive-f6d46b0898ecb21-20240513T173925433")
raster_list <- list(south, wet, slope, DEM)
names(raster_list[[1]]) <- "Snodgrass_aspect_southness_1m"
terra::varnames(raster_list[[1]]) <- "Snodgrass_aspect_southness_1m"
names(raster_list[[2]]) <- "Snodgrass_wetness_index_1m"
terra::varnames(raster_list[[2]]) <- "Snodgrass_wetness_index_1m"
names(raster_list[[3]]) <- "Snodgrass_slope_1m"
terra::varnames(raster_list[[3]]) <- "Snodgrass_slope_1m"
names(raster_list[[4]]) <- "Snodgrass_DEM_1m"
terra::varnames(raster_list[[4]]) <- "Snodgrass_DEM_1m"
cropped_rasters <- lapply(raster_list, terra::crop, y = new_ext)

# Save the cropped rasters
terra::writeRaster(cropped_rasters[[1]], "./inst/extdata/Snodgrass_aspect_southness_1m.tif", overwrite = TRUE)
terra::writeRaster(cropped_rasters[[2]], "./inst/extdata/Snodgrass_wetness_index_1m.tif", overwrite = TRUE)
terra::writeRaster(cropped_rasters[[3]], "./inst/extdata/Snodgrass_slope_1m.tif", overwrite = TRUE)
terra::writeRaster(cropped_rasters[[4]], "./inst/extdata/Snodgrass_DEM_1m.tif", overwrite = TRUE)


# Specify the spatial extent for the real example data
a_x <- 326496
a_y <- 4311439
b_x <- 326596 - 50
b_y <- 4311539 - 50

# Load the raster images
south <- terra::rast("https://data.ess-dive.lbl.gov/catalog/d1/mn/v2/object/ess-dive-c3c46ff25d50885-20240513T173925432")
wet <- terra::rast("https://data.ess-dive.lbl.gov/catalog/d1/mn/v2/object/ess-dive-dffdeec81023d23-20240513T173925427")
slope <- terra::rast("https://data.ess-dive.lbl.gov/catalog/d1/mn/v2/object/ess-dive-8a59960f4ffd550-20240513T173925429")
DEM <- terra::rast("https://data.ess-dive.lbl.gov/catalog/d1/mn/v2/object/ess-dive-f6d46b0898ecb21-20240513T173925433")
raster_list <- list(south, wet, slope, DEM)
names(raster_list[[1]]) <- "Snodgrass_aspect_southness_1m"
terra::varnames(raster_list[[1]]) <- "Snodgrass_aspect_southness_1m"
names(raster_list[[2]]) <- "Snodgrass_wetness_index_1m"
terra::varnames(raster_list[[2]]) <- "Snodgrass_wetness_index_1m"
names(raster_list[[3]]) <- "Snodgrass_slope_1m"
terra::varnames(raster_list[[3]]) <- "Snodgrass_slope_1m"
names(raster_list[[4]]) <- "Snodgrass_DEM_1m"
terra::varnames(raster_list[[4]]) <- "Snodgrass_DEM_1m"
cropped_rasters <- lapply(raster_list, terra::crop, y = new_ext)

# Save the cropped rasters for the real data example
terra::writeRaster(cropped_rasters[[1]], "./inst/extdata/Snodgrass_aspect_southness_1m_med.tif", overwrite = TRUE)
terra::writeRaster(cropped_rasters[[2]], "./inst/extdata/Snodgrass_wetness_index_1m_med.tif", overwrite = TRUE)
terra::writeRaster(cropped_rasters[[3]], "./inst/extdata/Snodgrass_slope_1m_med.tif", overwrite = TRUE)
terra::writeRaster(cropped_rasters[[4]], "./inst/extdata/Snodgrass_DEM_1m_med.tif", overwrite = TRUE)
