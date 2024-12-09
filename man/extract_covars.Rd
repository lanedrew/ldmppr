% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_covars.R
\name{extract_covars}
\alias{extract_covars}
\title{Extract covariate values from a set of rasters}
\usage{
extract_covars(locations, raster_list)
}
\arguments{
\item{locations}{a data frame of (x,y) locations with names "x" and "y".}

\item{raster_list}{a list of raster objects.}
}
\value{
a matrix of covariates drawn from the provided rasters.
}
\description{
Extract covariate values from a set of rasters
}
\examples{
# Load example raster data
raster_paths <- list.files(system.file("extdata", package = "ldmppr"),
  pattern = "\\\\.tif$", full.names = TRUE
)
rasters <- lapply(raster_paths, terra::rast)

# Scale the rasters
scaled_raster_list <- scale_rasters(rasters)

# Load example locations
locations <- small_example_data \%>\%
  dplyr::select(x, y) \%>\%
  as.matrix()

# Extract covariates
example_covars <- extract_covars(locations, scaled_raster_list)
head(example_covars)

}
