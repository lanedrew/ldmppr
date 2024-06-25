#' Plot simulated data realization
#'
#' @param sim_data ppp object with marks or dataframe with (x, y, size)
#'
#' @return ggplot object of the simulated data realization
#' @export
#'
plot_sim <- function(sim_data){
  sim_plot <- base::as.data.frame(sim_data) |>
    ggplot2::ggplot(ggplot2::aes(x = x, y = y, size = marks)) +
    ggplot2::geom_point(alpha = .5) +
    ggplot2::labs(x = "", y = "", title = "Simulated Data") +
    ggplot2::theme_bw()
  return(sim_plot)
}

#' Plot reference data
#'
#' @param ref_data ppp object with marks or dataframe with (x, y, size)
#'
#' @return ggplot object of the reference data
#' @export
#'
plot_refdata <- function(ref_data){
  ref_plot <- base::as.data.frame(ref_data) |>
    ggplot2::ggplot(ggplot2::aes(x = x, y = y, size = marks)) +
    ggplot2::geom_point(alpha = .5) +
    ggplot2::labs(x = "", y = "", title = "Reference Data") +
    ggplot2::theme_bw()
  return(ref_plot)
}
