#' Simulate from the self-correcting model
#'
#' @description
#' Allows the user to simulate a realization from the self-correcting model
#' given a set of parameters and a point to condition on.
#'
#' @param t_min minimum value for time.
#' @param t_max maximum value for time.
#' @param sc_params vector of parameter values corresponding to (alpha_1, beta_1, gamma_1, alpha_2, beta_2, alpha_3, beta_3, gamma_3).
#' @param anchor_point vector of (x,y) coordinates of point to condition on.
#' @param xy_bounds a vector of domain bounds (2 for x, 2 for y).
#'
#' @return a list containing the thinned and unthinned simulation realizations.
#' @export
#'
simulate_sc <- function(t_min = 0,
                        t_max = 1,
                        sc_params = NULL,
                        anchor_point = NULL,
                        xy_bounds = NULL){

  # Check the arguments
  if(t_min < 0 | t_min >= t_max | is.null(t_min)) stop("Provide a value greater than 0 and less than t_max for the t_min argument.")
  if(t_max > 1 | t_min >= t_max | is.null(t_max)) stop("Provide a value greater than t_min and less than 1 for the t_max argument.")
  if(length(sc_params) != 8 | anyNA(sc_params) | any(sc_params[2:8] < 0)) stop("Provide a valid set of parameter values for the sc_params argument.")
  if(length(anchor_point) != 2) stop("Provide a vector of (x,y) coordinates for the anchor_point argument.")
  if(is.null(xy_bounds) | !(length(xy_bounds) == 4)) stop("Provide (x, y) bounds in the form (a_x, b_x, a_y, b_y) for the xy_bounds argument.")
  if(xy_bounds[2] > xy_bounds[1] | xy_bounds[4] > xy_bounds[3]) stop("Provide (x, y) bounds in the form (a_x, b_x, a_y, b_y) for the xy_bounds argument.")


  # Simulate times and locations
  sim_times = stats::na.omit(sim_temporal_sc_cpp(t_min, t_max, sc_params[1:3]))
  sim_locs =  sim_spatial_sc_cpp(anchor_point, sc_params[4:5], length(sim_times), xy_bounds)
  sim_times[1] = 0
  txy_sim = base::cbind(sim_times, sim_locs)

  # Perform the thinning process
  thin_vals = (stats::runif(base::nrow(txy_sim), 0, 1) < interactionCpp_st(txy_sim, sc_params[6:8]))
  txy_sim_thin = txy_sim[thin_vals,]

  # Compile the thinned and unthinned results
  sim_df <- base::data.frame(time = txy_sim[,1], x = txy_sim[,2], y = txy_sim[,3])
  sim_thin_df <- base::data.frame(time = txy_sim_thin[,1], x = txy_sim_thin[,2], y = txy_sim_thin[,3])

  return(base::list(unthinned = sim_df, thinned = sim_thin_df))
}

