% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_mpp.R
\name{generate_mpp}
\alias{generate_mpp}
\title{Generate a marked process given locations and marks}
\usage{
generate_mpp(locations, marks = NULL, xy_bounds = NULL)
}
\arguments{
\item{locations}{a data frame of (x,y) locations with names "x" and "y".}

\item{marks}{a vector of marks.}

\item{xy_bounds}{a vector of domain bounds (2 for x, 2 for y).}
}
\value{
a ppp object with marks.
}
\description{
Creates an object of class "ppp" that represents a marked point pattern in the two-dimensional plane.
}
\examples{
# Load example data
data(small_example_data)

# Generate a marked point process
generate_mpp(
  locations = small_example_data \%>\% dplyr::select(x, y),
  marks = small_example_data$size,
  xy_bounds = c(0, 25, 0, 25)
)

}
