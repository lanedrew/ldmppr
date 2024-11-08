#' Check the fit of estimated self-correcting model on the reference point pattern dataset
#'
#' @description
#' Allows the user to perform global envelope tests for the L, F, G, J, E, and V functions from the \code{spatstat} package.
#' These tests serve as a goodness of fit measure for the estimated model relative to the reference dataset of interest.
#'
#'
#' @param reference_data a ppp object for the reference dataset.
#' @param t_min minimum value for time.
#' @param t_max maximum value for time.
#' @param sc_params vector of parameter values corresponding to (alpha_1, beta_1, gamma_1, alpha_2, beta_2, alpha_3, beta_3, gamma_3).
#' @param anchor_point vector of (x,y) coordinates of point to condition on.
#' @param raster_list a list of raster objects.
#' @param mark_model a model object (typically from \code{train_mark_model}).
#' @param xy_bounds a vector of domain bounds (2 for x, 2 for y).
#' @param include_comp_inds `TRUE` or `FALSE` indicating whether to generate and use competition indices as covariates.
#' @param competition_radius distance for competition radius if \code{include_comp_inds} is `TRUE`.
#' @param thinning `TRUE` or `FALSE` indicating whether to use the thinned or unthinned simulated values.
#' @param correction type of correction to apply ("none" or "toroidal").
#' @param n_sim number of simulated datasets to generate.
#' @param save_sims `TRUE` or `FALSE` indicating whether to save and return the simulated datasets.
#' @param verbose `TRUE` or `FALSE` indicating whether to show progress of model checking.
#'
#' @return a list containing various model fit summaries and
#' @export
#'
check_model_fit <- function(reference_data,
                            t_min = 0,
                            t_max = 1,
                            sc_params = NULL,
                            anchor_point = NULL,
                            raster_list = NULL,
                            mark_model = NULL,
                            xy_bounds = NULL,
                            include_comp_inds = FALSE,
                            thinning = TRUE,
                            correction = "none",
                            competition_radius = 15,
                            n_sim = 2500,
                            save_sims = TRUE,
                            verbose = TRUE){

  # Check the arguments
  if(!spatstat.geom::is.ppp(reference_data)) stop("Provide a ppp object containing the reference data pattern for the reference_data argument.")
  if(t_min < 0 | t_min >= t_max | is.null(t_min)) stop("Provide a value greater than 0 and less than t_max for the t_min argument.")
  if(t_max > 1 | t_min >= t_max | is.null(t_max)) stop("Provide a value greater than t_min and less than 1 for the t_max argument.")
  if(length(sc_params) != 8 | anyNA(sc_params) | any(sc_params[2:8] < 0)) stop("Provide a valid set of parameter values for the sc_params argument.")
  if(length(anchor_point) != 2) stop("Provide a vector of (x,y) coordinates for the anchor_point argument.")
  if(is.null(mark_model)) stop("Provide an unbundled mark model for the mark_model argument.")
  if(is.null(xy_bounds) | !(length(xy_bounds) == 4)) stop("Provide (x, y) bounds in the form (a_x, b_x, a_y, b_y) for the xy_bounds argument.")
  if(xy_bounds[2] > xy_bounds[1] | xy_bounds[4] > xy_bounds[3]) stop("Provide (x, y) bounds in the form (a_x, b_x, a_y, b_y) for the xy_bounds argument.")
  if(!correction %in% c("none", "toroidal")) stop("Provide a valid correction type for the correction argument.")
  if(include_comp_inds == TRUE & (is.null(competition_radius) | competition_radius < 0)) stop("Provide the desired radius for competition_indices argument.")


  # Obtain the radius to use in obtaining the various estimates using the reference data
  K_ref <- spatstat.explore::Kest(spatstat.geom::unmark(reference_data))
  d <- K_ref$r
  d_length <- base::length(d)

  # Create storage for the simulation results for the various metrics
  K_PP <- base::matrix(0, nrow = d_length, ncol = n_sim)
  F_PP <- base::matrix(0, nrow = d_length, ncol = n_sim)
  G_PP <- base::matrix(0, nrow = d_length, ncol = n_sim)
  J_PP <- base::matrix(0, nrow = d_length, ncol = n_sim)
  E_PP <- base::matrix(0, nrow = d_length, ncol = n_sim)
  V_PP <- base::matrix(0, nrow = d_length, ncol = n_sim)
  n_real <- base::numeric()

  if(verbose) {
    # Create the progress bar
    pb <- progress::progress_bar$new(
      format = "Data Simulation Progress: [:bar] :percent in :elapsed, ETA: :eta",
      total = n_sim, clear = FALSE, width = 100)

    base::message("Beginning Data Simulations...")
    pb$tick(0)
    base::Sys.sleep(3)

    for (j in 1:n_sim){

      # Simulate a dataset and predict the marks
      if(thinning){
        sim_j <- simulate_sc(t_min = t_min,
                             t_max = t_max,
                             sc_params = sc_params,
                             anchor_point = anchor_point,
                             xy_bounds = xy_bounds)$thinned
      } else{
        sim_j <- simulate_sc(t_min = t_min,
                             t_max = t_max,
                             sc_params = sc_params,
                             anchor_point = anchor_point,
                             xy_bounds = xy_bounds)$unthinned
      }
      pred_marks_j <- predict_marks(sim_realization = sim_j,
                                    raster_list = scale_rasters(raster_list),
                                    mark_model = mark_model,
                                    xy_bounds = xy_bounds,
                                    include_comp_inds = include_comp_inds,
                                    competition_radius = competition_radius,
                                    correction = correction)

      # Generate a ppp object from the locations and marks
      PP_xy <- generate_mpp(locations = sim_j, marks = pred_marks_j, xy_bounds = xy_bounds)

      # Calculate the point pattern metrics of interest
      K_PP[, j] <- spatstat.explore::Kest(spatstat.geom::unmark(PP_xy), correction = "isotropic", r = d)$iso
      F_PP[, j] <- spatstat.explore::Fest(spatstat.geom::unmark(PP_xy), correction = "rs",        r = d)$rs
      G_PP[, j] <- spatstat.explore::Gest(spatstat.geom::unmark(PP_xy), correction = "rs",        r = d)$rs
      J_PP[, j] <- spatstat.explore::Jest(spatstat.geom::unmark(PP_xy), correction = "rs",        r = d)$rs - 1
      E_PP[, j] <- spatstat.explore::Emark(PP_xy, correction = "isotropic", r = d)$iso
      V_PP[, j] <- spatstat.explore::Vmark(PP_xy, correction = "isotropic", r = d)$iso

      pb$tick()
      base::Sys.sleep(1 / 100)

    }
    base::message("Data Simulations Complete...")

  }else {

    for (j in 1:n_sim){

      # Simulate a dataset and predict the marks
      if(thinning){
        sim_j <- simulate_sc(t_min = t_min,
                             t_max = t_max,
                             sc_params = sc_params,
                             anchor_point = anchor_point,
                             xy_bounds = xy_bounds)$thinned
      } else{
        sim_j <- simulate_sc(t_min = t_min,
                             t_max = t_max,
                             sc_params = sc_params,
                             anchor_point = anchor_point,
                             xy_bounds = xy_bounds)$unthinned
      }
      pred_marks_j <- predict_marks(sim_realization = sim_j,
                                    raster_list = scale_rasters(raster_list),
                                    mark_model = mark_model,
                                    xy_bounds = xy_bounds,
                                    include_comp_inds = include_comp_inds,
                                    competition_radius = competition_radius,
                                    correction = correction)

      # Generate a ppp object from the locations and marks
      PP_xy <- generate_mpp(locations = sim_j, marks = pred_marks_j, xy_bounds = xy_bounds)
      n_real[j] <- base::length(sim_j[,1])

      # Calculate the point pattern metrics of interest
      K_PP[, j] <- spatstat.explore::Kest(spatstat.geom::unmark(PP_xy), correction = "isotropic", r = d)$iso
      F_PP[, j] <- spatstat.explore::Fest(spatstat.geom::unmark(PP_xy), correction = "rs",        r = d)$rs
      G_PP[, j] <- spatstat.explore::Gest(spatstat.geom::unmark(PP_xy), correction = "rs",        r = d)$rs
      J_PP[, j] <- spatstat.explore::Jest(spatstat.geom::unmark(PP_xy), correction = "rs",        r = d)$rs - 1
      E_PP[, j] <- spatstat.explore::Emark(PP_xy, correction = "isotropic", r = d)$iso
      V_PP[, j] <- spatstat.explore::Vmark(PP_xy, correction = "isotropic", r = d)$iso

    }

  }

  if(save_sims){
    sim_list <- base::list(Ksim = K_PP,
                           Fsim = F_PP,
                           Gsim = G_PP,
                           Jsim = J_PP,
                           Esim = E_PP,
                           Vsim = V_PP,
                           n_per = n_real)
  }


  # Obtain the global envelope test for the L function
  C_ref_L <- GET::create_curve_set(base::list(r = d,
                                              obs = sqrt(K_ref$iso / pi) - d,
                                              theo = sqrt(K_ref$theo / pi) - d ,
                                              sim_m = sqrt(K_PP / pi) - d)
                                   )
  r_envL <- GET::global_envelope_test(C_ref_L, type = "rank")

  # Obtain the global envelope test for the F function
  F_val <- base::max(base::apply(F_PP, 2, function(x) base::min(base::which(x >= 1))))
  C_ref_F <- GET::create_curve_set(base::list(r = d[1:F_val],
                                              obs = spatstat.explore::Fest(spatstat.geom::unmark(reference_data),
                                                                           r = d[1:F_val])$rs,
                                              theo = spatstat.explore::Fest(spatstat.geom::unmark(reference_data),
                                                                            r = d[1:F_val])$theo,
                                              sim_m = F_PP[1:F_val,])
                                   )
  r_envF <- GET::global_envelope_test(C_ref_F, type = "rank")

  # Obtain the global envelope test for the G function
  G_val <- base::max(base::apply(G_PP, 2, function(x) base::min(base::which(x >= 1))))
  C_ref_G <- GET::create_curve_set(base::list(r = d[1:G_val],
                                              obs = spatstat.explore::Gest(spatstat.geom::unmark(reference_data),
                                                                           r = d[1:G_val])$rs,
                                              theo = spatstat.explore::Gest(spatstat.geom::unmark(reference_data),
                                                                           r = d[1:G_val])$theo,
                                              sim_m = G_PP[1:G_val,])
                                   )
  r_envG <- GET::global_envelope_test(C_ref_G, type = "rank")

  # Obtain the global envelope test for the J function
  J_val <- base::min(c(base::min(base::apply(F_PP, 2, function(x) base::sum(x < 1))),
                       base::min(base::apply(G_PP, 2, function(x) base::sum(x < 1))),
                       base::sum(!base::is.na(spatstat.explore::Jest(spatstat.geom::unmark(reference_data),
                                                                     r = d)$rs))))

  if (any(is.infinite(J_PP[1:J_val,]) | is.na(J_PP[1:J_val,]))) {
    warning("J_PP contains Inf or NA values.")
  }
  C_ref_J <- GET::create_curve_set(base::list(r = d[1:J_val],
                                              obs = spatstat.explore::Jest(spatstat.geom::unmark(reference_data),
                                                                           r = d[1:J_val])$rs - 1,
                                              theo = spatstat.explore::Jest(spatstat.geom::unmark(reference_data),
                                                                           r = d[1:J_val])$theo - 1,
                                              sim_m = J_PP[1:J_val,]),
                                   allfinite = FALSE
                                   )
  r_envJ <- GET::global_envelope_test(C_ref_J, type = "rank")

  # Obtain the global envelope test for the E function
  C_ref_E <- GET::create_curve_set(base::list(r = d,
                               obs = spatstat.explore::Emark(reference_data, correction = "isotropic", r = d)$iso,
                               theo = spatstat.explore::Emark(reference_data, correction = "isotropic", r = d)$theo,
                               sim_m = E_PP))
  r_envE <- GET::global_envelope_test(C_ref_E, type = "rank")

  # Obtain the global envelope test for the V function
  C_ref_V <- GET::create_curve_set(base::list(r = d,
                               obs = spatstat.explore::Vmark(reference_data, correction = "isotropic", r = d)$iso,
                               theo = spatstat.explore::Vmark(reference_data, correction = "isotropic", r = d)$theo,
                               sim_m = V_PP))
  r_envV <- GET::global_envelope_test(C_ref_V, type = "rank")

  # Obtain the global envelope test for the combined L, F, G, and J functions
  rComb <- c(d,
             d[1:F_val] + base::max(d),
             d[1:G_val] + base::max(d[1:F_val] + base::max(d)),
             d[1:J_val] + base::max(d[1:G_val] + base::max(d[1:F_val] + base::max(d)))
  )
  K_data <- spatstat.explore::Kest(reference_data,
                                   correction = "isotropic",
                                   r = d)$iso
  F_data <- spatstat.explore::Fest(reference_data,
                                   correction = "rs",
                                   r = d[1:F_val])$rs
  G_data <- spatstat.explore::Gest(reference_data,
                                   correction = "rs",
                                   r = d[1:G_val])$rs
  J_data <- spatstat.explore::Jest(reference_data,
                                   correction = "rs",
                                   r = d[1:J_val])$rs - 1
  J_scale <- base::max(base::apply(J_PP[1:J_val,], 1, base::max))
  Comb_data <- c(sqrt(K_data / pi) - d,
                 F_data,
                 G_data,
                 J_data / J_scale)
  Comb_ref <- GET::create_curve_set(base::list(r = rComb,
                                               obs = Comb_data,
                                               sim_m = base::rbind(sqrt(K_PP / pi) - d,
                                                                   F_PP[1:F_val,],
                                                                   G_PP[1:G_val,],
                                                                   J_PP[1:J_val,] / J_scale)
                                               )
                                   )
  r_envComb <-  GET::global_envelope_test(Comb_ref, type = "rank")

  if(save_sims){
    results_list <- base::list(global_envelope_tests = list(L_env = r_envL,
                                                            F_env = r_envF,
                                                            G_env = r_envG,
                                                            J_env = r_envJ,
                                                            E_env = r_envE,
                                                            V_env = r_envV,
                                                            Combined_env = r_envComb),
                               sim_metric_values = sim_list)
  } else{
    results_list <- base::list(L_env = r_envL,
                               F_env = r_envF,
                               G_env = r_envG,
                               J_env = r_envJ,
                               E_env = r_envE,
                               V_env = r_envV,
                               Combined_env = r_envComb)
  }

  return(results_list)

}
