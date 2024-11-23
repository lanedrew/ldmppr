#' Gentle decay (power-law) mapping function from sizes to arrival times
#'
#' @param sizes vector of sizes to be mapped to arrival times.
#' @param delta numeric value (greater than 0) for the exponent in the mapping function.
#'
#' @return vector of arrival times.
#' @export
#'
#' @examples
#' # Generate a vector of sizes
#' sizes <- runif(100, 0, 100)
#'
#' # Map the sizes to arrival times using a power-law mapping with delta = .5
#' power_law_mapping(sizes, .5)
#'
power_law_mapping <- function(sizes, delta) {
  # Check if delta value is within the appropriate support
  if (delta < 0) {
    stop("Please input a value for delta greater than 0.")
  }

  # Calculate the arrival times given sizes and delta
  t_scaled <- 1 - ((sizes - base::min(sizes)) / (base::max(sizes) - base::min(sizes)))^delta
  return(t_scaled)
}
