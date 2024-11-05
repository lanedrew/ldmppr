#' Gentle decay (power-law) mapping function from sizes to arrival times
#'
#' @param sizes vector of sizes to be mapped into arrival times
#' @param delta numeric value (greater than 0) for the exponent in the mapping function
#'
#' @return vector of arrival times
#' @export
#'
power_law_mapping <- function(sizes, delta) {
  # Check if delta value is within the appropriate support
  if(delta < 0){
    stop("Please input a value for delta greater than 0.")
  }

  # Calculate the arrival times given sizes and delta
  t_scaled <- 1 - ((sizes - base::min(sizes)) / (base::max(sizes) - base::min(sizes)))^delta
  return(t_scaled)
}
