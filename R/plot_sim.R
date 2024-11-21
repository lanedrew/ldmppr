#' Plot a marked point process
#'
#' @param mpp_data ppp object with marks or data frame with columns (x, y, size).
#' @param pattern_type type of pattern to plot ("reference" or "simulated").
#'
#' @return a ggplot object of the marked point process.
#' @export
#'
#' @examples
#' # Load example data
#' data(small_example_data)
#' mpp_data <- generate_mpp(locations = small_example_data %>% dplyr::select(x, y),
#'                          marks = small_example_data$size,
#'                          xy_bounds = c(0, 25, 0, 25))
#'
#' # Plot the marked point process
#' plot_mpp(mpp_data, pattern_type = "reference")
#'
plot_mpp <- function(mpp_data, pattern_type){

  if(pattern_type == "simulated"){
    plot_title <- "Simulated Data"
  } else if(pattern_type == "reference"){
    plot_title <- "Reference Data"
  } else {
    stop("Provide a valid pattern type ('reference' or 'simulated').", .call = FALSE)
  }

  mpp_plot <- base::as.data.frame(mpp_data) |>
    ggplot2::ggplot(ggplot2::aes(x = mpp_data$x, y = mpp_data$y, size = mpp_data$marks)) +
    ggplot2::geom_point(alpha = .5) +
    ggplot2::labs(x = "", y = "", title = plot_title, size = "Mark") +
    ggplot2::theme_bw()
  return(mpp_plot)
}
