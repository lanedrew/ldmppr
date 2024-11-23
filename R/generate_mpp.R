#' Generate a marked process given locations and marks
#'
#' @description
#' Creates an object of class "ppp" that represents a marked point pattern in the two-dimensional plane.
#'
#' @param locations a data frame of (x,y) locations with names "x" and "y".
#' @param marks a vector of marks.
#' @param xy_bounds a vector of domain bounds (2 for x, 2 for y).
#'
#' @return a ppp object with marks.
#' @export
#'
#' @examples
#' # Load example data
#' data(small_example_data)
#'
#' # Generate a marked point process
#' generate_mpp(
#'   locations = small_example_data %>% dplyr::select(x, y),
#'   marks = small_example_data$size,
#'   xy_bounds = c(0, 25, 0, 25)
#' )
#'
generate_mpp <- function(locations, marks = NULL, xy_bounds = NULL) {
  # Check arguments
  if (!is.data.frame(locations)) stop("Provide a data frame with columns x and y for the locations argument.", .call = FALSE)
  if (is.null(marks)) stop("Provide a vector of marks with length matching the number of locations provided for the marks argument.", .call = FALSE)
  if (is.null(xy_bounds) | !(length(xy_bounds) == 4)) stop("Provide (x,y) bounds in the form (a_x, b_x, a_y, b_y) for the xy_bounds argument.", .call = FALSE)
  if (xy_bounds[2] < xy_bounds[1] | xy_bounds[4] < xy_bounds[3]) stop("Provide (x,y) bounds in the form (a_x, b_x, a_y, b_y) for the xy_bounds argument.", .call = FALSE)

  # Create a point process object with marks
  marked_pp <- spatstat.geom::ppp(locations$x,
    locations$y,
    xy_bounds[1:2], xy_bounds[3:4],
    marks = marks
  )

  return(marked_pp)
}
