#' Create a grid schedule for estimate_process_parameters()
#'
#' \code{ldmppr_grids()} constructs a multi-resolution grid schedule used by
#' \code{\link{estimate_process_parameters}}. The returned object is an S3 class
#' with helper methods; see \code{\link{ldmppr_grids-class}}.
#'
#' @param upper_bounds a vector \code{c(b_t, b_x, b_y)} giving the
#'   maximum bounds for time and the spatial domain. Grids must lie within these.
#' @param levels a list describing the grid schedule. Each entry can be either:
#'   \itemize{
#'     \item a numeric length-3 vector \code{c(nx, ny, nt)} (number of points per dimension),
#'     \item a list with elements \code{nx}, \code{ny}, \code{nt},
#'     \item a list with explicit vectors \code{x}, \code{y}, \code{t}.
#'   }
#' @param labels (optional) character vector of length equal to \code{levels}, used
#'   only for printing.
#' @param include_endpoints \code{TRUE} or \code{FALSE} indicating; if \code{TRUE} (default) each generated grid uses
#'   \code{seq(0, bound, length.out = n)} including endpoints.
#'
#' @return an object of class \code{"ldmppr_grids"}.
#' @export
#'
#' @seealso \code{\link{ldmppr_grids-class}} for methods and details.
#'
#' @examples
#' # A 3-level coarse-to-fine schedule (counts per dimension)
#' g <- ldmppr_grids(
#'   upper_bounds = c(1, 50, 50),
#'   levels = list(
#'     c(25, 25, 25),
#'     c(75, 75, 75),
#'     c(100, 100, 100)
#'   )
#' )
#' g
#' length(g)
#' summary(g)
#'
#' # Explicit vectors (single level)
#' g2 <- ldmppr_grids(
#'   upper_bounds = c(1, 50, 50),
#'   levels = list(list(
#'     x = seq(0, 50, by = 2),
#'     y = seq(0, 50, by = 2),
#'     t = seq(0, 1,  length.out = 30)
#'   ))
#' )
#' as.data.frame(g2)
ldmppr_grids <- function(upper_bounds,
                         levels,
                         labels = NULL,
                         include_endpoints = TRUE) {

  grids <- .ldmppr_make_grid_schedule(
    upper_bounds = upper_bounds,
    levels = levels,
    labels = labels,
    include_endpoints = include_endpoints
  )

  new_ldmppr_grids(
    levels = grids$levels,
    upper_bounds = grids$upper_bounds,
    labels = grids$labels,
    include_endpoints = grids$include_endpoints
  )
}
