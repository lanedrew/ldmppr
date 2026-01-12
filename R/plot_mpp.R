#' Plot a marked point process
#'
#' @param mpp_data \code{ppp} object with marks or data frame with columns (x, y, size).
#' @param pattern_type type of pattern to plot ("reference" or "simulated").
#'
#' @return a \code{ggplot} object of the marked point process.
#'
#' @examples
#' # Load example data
#' data(small_example_data)
#' mpp_data <- generate_mpp(
#'   locations = small_example_data %>% dplyr::select(x, y),
#'   marks = small_example_data$size,
#'   xy_bounds = c(0, 25, 0, 25)
#' )
#'
#' # Plot the marked point process
#' plot_mpp(mpp_data, pattern_type = "reference")
#'
plot_mpp <- function(mpp_data, pattern_type = c("reference", "simulated")) {
  pattern_type <- match.arg(pattern_type)
  plot_title <- if (pattern_type == "simulated") "Simulated Data" else "Reference Data"

  df <- base::as.data.frame(mpp_data)

  # Handle the common "size" naming too
  if (!("marks" %in% names(df)) && ("size" %in% names(df))) df$marks <- df$size

  if (!all(c("x", "y", "marks") %in% names(df))) {
    stop("`mpp_data` must contain columns x, y, and marks (or size).", call. = FALSE)
  }

  ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, size = marks)) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::labs(x = "", y = "", title = plot_title, size = "Mark")
}
