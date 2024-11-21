#' Optimized function to compute toroidal distance matrix over a rectangular domain
#'
#' @param location_matrix a 2 column matrix of (x,y) coordinates.
#' @param x_bound the upper bound for the x dimension.
#' @param y_bound the upper bound for the y dimension.
#'
#' @return a matrix of toroidal distances.
#' @export
#'
#' @examples
#' # Generate a matrix of locations
#' location_matrix <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 2)
#' x_bound <- 10
#' y_bound <- 10
#'
#' # Compute the toroidal distance matrix
#' toroidal_dist_matrix_optimized(location_matrix, x_bound, y_bound)
#'
toroidal_dist_matrix_optimized <- function(location_matrix, x_bound, y_bound) {
  # Extract x and y coordinates
  x <- location_matrix[, 1]
  y <- location_matrix[, 2]

  # Number of points
  n <- length(x)

  # Compute pairwise differences for x and y using outer
  dx <- outer(x, x, function(a, b) pmin(abs(a - b), x_bound - abs(a - b)))
  dy <- outer(y, y, function(a, b) pmin(abs(a - b), y_bound - abs(a - b)))

  # Compute the distance matrix
  dist_matrix <- sqrt(dx^2 + dy^2)

  # Return the distance matrix
  return(dist_matrix)
}
