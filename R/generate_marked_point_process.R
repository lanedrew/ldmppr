#' Generate a marked process given locations and marks
#'
#' @param location_data a set of locations
#' @param marks a vector of marks
#' @param window a vector of window bounds (2 for x, 2 for y)
#'
#' @return a ppp object with marks
#' @export
#'
generate_mpp <- function(location_data, marks, window){

  location_data <- base::as.data.frame(location_data)
  marked_pp <- spatstat.geom::ppp(location_data[,c("x")],
                                  location_data[,c("y")],
                                  window[1:2], window[3:4],
                                  marks = marks)

  return(marked_pp)
}
