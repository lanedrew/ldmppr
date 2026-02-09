#' Simulate from the self-correcting model
#'
#' @description
#' Allows the user to simulate a realization from the self-correcting model
#' given a set of parameters and a point to condition on.
#'
#' @param t_min minimum value for time.
#' @param t_max maximum value for time.
#' @param sc_params a vector of parameter values corresponding to
#'   \eqn{(\alpha_1,\beta_1,\gamma_1,\alpha_2,\beta_2,\alpha_3,\beta_3,\gamma_3)}
#'   (i.e., alpha_1, beta_1, gamma_1, alpha_2, beta_2, alpha_3, beta_3, gamma_3).
#' @param anchor_point a vector of (x,y) coordinates of point to condition on.
#' @param xy_bounds a vector of domain bounds (2 for x, 2 for y).
#'
#' @return a list containing the thinned and unthinned simulation realizations.
#' @export
#'
#' @examples
#' # Specify the generating parameters of the self-correcting process
#' generating_parameters <- c(2, 8, .02, 2.5, 3, 1, 2.5, .2)
#'
#' # Specify an anchor point
#' M_n <- c(10, 14)
#'
#' # Simulate the self-correcting process
#' generated_locs <- simulate_sc(
#'   t_min = 0,
#'   t_max = 1,
#'   sc_params = generating_parameters,
#'   anchor_point = M_n,
#'   xy_bounds = c(0, 25, 0, 25)
#' )
#'
simulate_sc <- function(t_min = 0,
                        t_max = 1,
                        sc_params = NULL,
                        anchor_point = NULL,
                        xy_bounds = NULL) {

  # Checks
  if (is.null(t_min) || t_min < 0 || t_min >= t_max) {
    stop("Provide t_min >= 0 and < t_max.", call. = FALSE)
  }
  if (is.null(t_max) || t_max > 1 || t_min >= t_max) {
    stop("Provide t_max > t_min and <= 1.", call. = FALSE)
  }
  if (is.null(sc_params) || length(sc_params) != 8 || anyNA(sc_params) || any(sc_params[2:8] < 0)) {
    stop("Provide a valid set of parameter values for sc_params.", call. = FALSE)
  }

  anchor_point <- as.numeric(anchor_point)
  if (length(anchor_point) != 2) {
    stop("Provide a vector/matrix of (x,y) coordinates for anchor_point.", call. = FALSE)
  }

  if (is.null(xy_bounds) || length(xy_bounds) != 4) {
    stop("Provide xy_bounds = c(a_x, b_x, a_y, b_y).", call. = FALSE)
  }
  if (xy_bounds[2] < xy_bounds[1] || xy_bounds[4] < xy_bounds[3]) {
    stop("Invalid xy_bounds ordering.", call. = FALSE)
  }

  # ---- Simulate times ----
  sim_times <- stats::na.omit(sim_temporal_sc(t_min, t_max, sc_params[1:3]))
  sim_times <- as.numeric(sim_times)

  # If no times, return empty realizations with correct columns
  if (length(sim_times) == 0) {
    empty <- data.frame(time = numeric(0), x = numeric(0), y = numeric(0))
    return(list(unthinned = empty, thinned = empty))
  }

  # Match your previous convention: force first time to 0
  sim_times[1] <- 0

  # ---- Simulate locations ----
  sim_locs <- sim_spatial_sc(anchor_point, sc_params[4:5], length(sim_times), xy_bounds)

  # Ensure matrix with 2 columns
  sim_locs <- as.matrix(sim_locs)
  if (ncol(sim_locs) != 2) {
    stop("sim_spatial_sc() did not return a 2-column matrix.", call. = FALSE)
  }

  txy_sim <- cbind(sim_times, sim_locs)
  colnames(txy_sim) <- c("time", "x", "y")

  # ---- Thinning ----
  keep <- stats::runif(nrow(txy_sim), 0, 1) < interaction_st_fast(txy_sim, sc_params[6:8])

  txy_sim_thin <- txy_sim[keep, , drop = FALSE]

  sim_df      <- as.data.frame(txy_sim)
  sim_thin_df <- as.data.frame(txy_sim_thin)

  # Always return data.frames with the same columns
  sim_df      <- sim_df[, c("time", "x", "y"), drop = FALSE]
  sim_thin_df <- sim_thin_df[, c("time", "x", "y"), drop = FALSE]

  list(unthinned = sim_df, thinned = sim_thin_df)
}
