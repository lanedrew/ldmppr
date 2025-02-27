% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_mpp.R
\name{simulate_mpp}
\alias{simulate_mpp}
\title{Simulate a realization of a location dependent marked point process}
\usage{
simulate_mpp(
  sc_params = NULL,
  t_min = 0,
  t_max = 1,
  anchor_point = NULL,
  raster_list = NULL,
  scaled_rasters = FALSE,
  mark_model = NULL,
  xy_bounds = NULL,
  include_comp_inds = FALSE,
  competition_radius = 15,
  correction = "none",
  thinning = TRUE
)
}
\arguments{
\item{sc_params}{vector of parameter values corresponding to (alpha_1, beta_1, gamma_1, alpha_2, beta_2, alpha_3, beta_3, gamma_3).}

\item{t_min}{minimum value for time.}

\item{t_max}{maximum value for time.}

\item{anchor_point}{vector of (x,y) coordinates of point to condition on.}

\item{raster_list}{list of raster objects.}

\item{scaled_rasters}{`TRUE` or `FALSE` indicating whether the rasters have been scaled.}

\item{mark_model}{a model object (typically from \code{train_mark_model}).}

\item{xy_bounds}{a vector of domain bounds (2 for x, 2 for y).}

\item{include_comp_inds}{`TRUE` or `FALSE` indicating whether to generate and use competition indices as covariates.}

\item{competition_radius}{distance for competition radius if \code{include_comp_inds} is `TRUE`.}

\item{correction}{type of correction to apply ("none" or "toroidal").}

\item{thinning}{`TRUE` or `FALSE` indicating whether to thin the realization.}
}
\value{
a list containing the marked point process realization and the data frame of the realization.
}
\description{
Simulate a realization of a location dependent marked point process
}
\examples{
# Specify the generating parameters of the self-correcting process
generating_parameters <- c(2, 8, .02, 2.5, 3, 1, 2.5, .2)

# Specify an anchor point
M_n <- matrix(c(10, 14), ncol = 1)

# Load the raster files
raster_paths <- list.files(system.file("extdata", package = "ldmppr"),
  pattern = "\\\\.tif$", full.names = TRUE
)
raster_paths <- raster_paths[!grepl("_med\\\\.tif$", raster_paths)]
rasters <- lapply(raster_paths, terra::rast)

# Scale the rasters
scaled_raster_list <- scale_rasters(rasters)

# Load the example mark model
file_path <- system.file("extdata", "example_mark_model.rds", package = "ldmppr")
mark_model <- bundle::unbundle(readRDS(file_path))

# Simulate a realization
example_mpp <- simulate_mpp(
  sc_params = generating_parameters,
  t_min = 0,
  t_max = 1,
  anchor_point = M_n,
  raster_list = scaled_raster_list,
  scaled_rasters = TRUE,
  mark_model = mark_model,
  xy_bounds = c(0, 25, 0, 25),
  include_comp_inds = TRUE,
  competition_radius = 10,
  correction = "none",
  thinning = TRUE
)

# Plot the realization
plot_mpp(example_mpp$mpp, pattern_type = "simulated")

}
