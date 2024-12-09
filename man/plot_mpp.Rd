% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_mpp.R
\name{plot_mpp}
\alias{plot_mpp}
\title{Plot a marked point process}
\usage{
plot_mpp(mpp_data, pattern_type)
}
\arguments{
\item{mpp_data}{ppp object with marks or data frame with columns (x, y, size).}

\item{pattern_type}{type of pattern to plot ("reference" or "simulated").}
}
\value{
a ggplot object of the marked point process.
}
\description{
Plot a marked point process
}
\examples{
# Load example data
data(small_example_data)
mpp_data <- generate_mpp(
  locations = small_example_data \%>\% dplyr::select(x, y),
  marks = small_example_data$size,
  xy_bounds = c(0, 25, 0, 25)
)

# Plot the marked point process
plot_mpp(mpp_data, pattern_type = "reference")

}
