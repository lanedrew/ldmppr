#' Plot simulated data realization
#'
#' @param sim_data ppp object with marks or dataframe with (x, y, size)
#'
#' @return ggplot object of the simulated data realization
#' @export
#'
plot_sim <- function(sim_data){
  sim_plot <- base::as.data.frame(sim_data) |>
    ggplot2::ggplot(ggplot2::aes(x = sim_data$x, y = sim_data$y, size = sim_data$marks)) +
    ggplot2::geom_point(alpha = .5) +
    ggplot2::labs(x = "", y = "", title = "Simulated Data", size = "Mark") +
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
    ggplot2::ggplot(ggplot2::aes(x = ref_data$x, y = ref_data$y, size = ref_data$marks)) +
    ggplot2::geom_point(alpha = .5) +
    ggplot2::labs(x = "", y = "", title = "Reference Data", size = "Mark") +
    ggplot2::theme_bw()
  return(ref_plot)
}
