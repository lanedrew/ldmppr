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
  if (length(delta) != 1 || is.na(delta) || delta < 0) {
    stop("Please input a value for delta >= 0.")
  }
  if (!is.numeric(sizes) || anyNA(sizes)) {
    stop("`sizes` must be a numeric vector with no missing values.")
  }

  smin <- base::min(sizes)
  smax <- base::max(sizes)

  if (smax == smin) {
    # All sizes equal: map everything to the middle (or 0.5) to avoid division by zero.
    return(rep(0.5, length(sizes)))
  }

  t_scaled <- 1 - ((sizes - smin) / (smax - smin))^delta
  return(t_scaled)
}

