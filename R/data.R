#' Small Example Data
#'
#' A small example dataset for testing and examples consisting of 121 observations in a (25m x 25m) square domain.
#'
#' @format ## `small_example_data`
#' A data frame with 121 rows and 3 columns:
#' \describe{
#'   \item{x}{x coordinate}
#'   \item{y}{y coordinate}
#'   \item{size}{Size}
#'   ...
#' }
#'
#' @details
#' The dataset was generated using the example raster data and an exponential decay size function.
#'
#' The full code to generate this dataset is available in the package's `data_raw` directory.
#'
#' @source
#' Simulated dataset. Code to generate it can be found in `data_raw/small_example_data.R`.
#'
#' @usage data("small_example_data")
"small_example_data"

#' Medium Example Data
#'
#' A medium sized example dataset consisting of 111 observations in a (50m x 50m) square domain.
#'
#' @format ## `medium_example_data`
#' A data frame with 111 rows and 3 columns:
#' \describe{
#'   \item{x}{x coordinate}
#'   \item{y}{y coordinate}
#'   \item{size}{Size}
#'   ...
#' }
#'
#' @details
#' The dataset was generated using the Snodgrass dataset available at https://data.ess-dive.lbl.gov/view/doi:10.15485/2476543.
#'
#' The full code to generate this dataset is available in the package's `data_raw` directory.
#'
#' @source
#' Real example dataset. Code to generate it can be found in `data_raw/medium_example_data.R`.
#'
#' @usage data("medium_example_data")
"medium_example_data"
