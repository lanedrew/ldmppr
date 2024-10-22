#' check fit of estimated model on the reference data
#'
#' @param ref_data a ppp object for the reference data
#' @param Tmin minimum value for T (defaults to 0)
#' @param Tmax maximum value for T
#' @param params vector of parameter estimates
#' @param M_n point to condition on
#' @param raster_list a list of raster objects
#' @param mark_model a trained mark model object
#' @param xy_bounds a vector of domain bounds (2 for x, 2 for y)
#' @param include_comp_inds `TRUE` or `FALSE` indicating whether to generate and use competition indices as covariates
#' @param correction type of correction to apply ("none" or "toroidal")
#' @param n_sim number of simulated datasets to generate
#' @param verbose `TRUE` or `FALSE` indicating whether to show progress
#'
#' @return a list of model fit summaries
#' @export
#'
check_model_fit <- function(ref_data, Tmin = 0, Tmax, params,
                            M_n, raster_list, mark_model,
                            xy_bounds,
                            include_comp_inds = FALSE,
                            correction = "none",
                            n_sim = 2500,
                            verbose = TRUE){

  d <- spatstat.explore::Kest(spatstat.geom::unmark(ref_data))$r
  d_length <- base::length(d)

  K_PP <- base::matrix(0, nrow = d_length, ncol = n_sim)
  F_PP <- base::matrix(0, nrow = d_length, ncol = n_sim)
  G_PP <- base::matrix(0, nrow = d_length, ncol = n_sim)
  J_PP <- base::matrix(0, nrow = d_length, ncol = n_sim)
  E_PP <- base::matrix(0, nrow = d_length, ncol = n_sim)
  V_PP <- base::matrix(0, nrow = d_length, ncol = n_sim)
  n_real <- base::numeric()

  if(verbose) print("Beginning simulations!")
  for (j in 1:n_sim){
    sim_j <- sim_spatial_temporal_sc_cpp(Tmin = 0, Tmax = 1, params, M_n, xy_bounds)[[2]]
    pred_marks_j <- predict_marks(sim_realization = sim_j,
                                  raster_list = raster_list,
                                  size_model = mark_model,
                                  xy_bounds = xy_bounds,
                                  include_comp_inds = include_comp_inds,
                                  correction = correction)
    PP_xy <- generate_mpp(location_data = sim_j, marks = pred_marks_j, window = xy_bounds)
    n_real[j] <- base::length(sim_j[,1])
    K_PP[, j] <- spatstat.explore::Kest(spatstat.geom::unmark(PP_xy), correction = "isotropic", r = d)$iso
    F_PP[, j] <- spatstat.explore::Fest(spatstat.geom::unmark(PP_xy), correction = "rs",        r = d)$rs
    G_PP[, j] <- spatstat.explore::Gest(spatstat.geom::unmark(PP_xy), correction = "rs",        r = d)$rs
    J_PP[, j] <- spatstat.explore::Jest(spatstat.geom::unmark(PP_xy), correction = "rs",        r = d)$rs - 1
    E_PP[, j] <- spatstat.explore::Emark(PP_xy, correction = "isotropic", r = d)$iso
    V_PP[, j] <- spatstat.explore::Vmark(PP_xy, correction = "isotropic", r = d)$iso

    if(verbose == TRUE & j %% 100 == 0){
      print(paste0("Simulation iterations complete: ", j, "/", n_sim))
    }

  }

  J_scale <- base::max(base::apply(J_PP, 1, base::max))
  sim_list <- base::list(Ksim = K_PP,
                         Fsim = F_PP,
                         Gsim = G_PP,
                         Jsim = J_PP,
                         Esim = E_PP,
                         Vsim = V_PP,
                         n_per = n_real)

  if(verbose) print("Simulations complete!")


  K_ref <- spatstat.explore::Kest(spatstat.geom::unmark(ref_data))
  d <- K_ref$r
  C_ref_L <- GET::create_curve_set(base::list(r = d,
                                              obs = sqrt(K_ref$iso / pi) - d,
                                              theo = sqrt(K_ref$theo / pi) - d ,
                                              sim_m = sqrt(K_PP / pi) - d)
                                   )
  r_envL <- GET::global_envelope_test(C_ref_L, type = "rank")

  F_val <- base::max(base::apply(F_PP, 2, function(x) base::min(base::which(x >= 1))))
  C_ref_F <- GET::create_curve_set(base::list(r = d[1:F_val],
                                              obs = spatstat.explore::Fest(spatstat.geom::unmark(ref_data),
                                                                           r = d[1:F_val])$rs,
                                              theo = spatstat.explore::Fest(spatstat.geom::unmark(ref_data),
                                                                            r = d[1:F_val])$theo,
                                              sim_m = F_PP[1:F_val,])
                                   )
  r_envF <- GET::global_envelope_test(C_ref_F, type = "rank")

  G_val <- base::max(base::apply(G_PP, 2, function(x) base::min(base::which(x >= 1))))
  C_ref_G <- GET::create_curve_set(base::list(r = d[1:G_val],
                                              obs = spatstat.explore::Gest(spatstat.geom::unmark(ref_data),
                                                                           r = d[1:G_val])$rs,
                                              theo = spatstat.explore::Gest(spatstat.geom::unmark(ref_data),
                                                                           r = d[1:G_val])$theo,
                                              sim_m = G_PP[1:G_val,])
                                   )
  r_envG <- GET::global_envelope_test(C_ref_G, type = "rank")

  J_val <- base::min(c(base::min(base::apply(F_PP, 2, function(x) base::sum(x < 1))),
                       base::min(base::apply(G_PP, 2, function(x) base::sum(x < 1)))))
  # J_val <- base::min(c(F_val - 1, G_val - 1))
  C_ref_J <- GET::create_curve_set(base::list(r = d[1:J_val],
                                              obs = spatstat.explore::Jest(spatstat.geom::unmark(ref_data),
                                                                           r = d[1:J_val])$rs,
                                              theo = spatstat.explore::Jest(spatstat.geom::unmark(ref_data),
                                                                           r = d[1:J_val])$theo,
                                              sim_m = J_PP[1:J_val,]),
                                   allfinite = FALSE
                                   )
  r_envJ <- GET::global_envelope_test(C_ref_J, type = "rank")

  C_ref_E <- GET::create_curve_set(base::list(r = d,
                               obs = spatstat.explore::Emark(ref_data, correction = "isotropic", r = d)$iso,
                               theo = spatstat.explore::Emark(ref_data, correction = "isotropic", r = d)$theo,
                               sim_m = E_PP))
  r_envE <- GET::global_envelope_test(C_ref_E, type = "rank")

  C_ref_V <- GET::create_curve_set(base::list(r = d,
                               obs = spatstat.explore::Vmark(ref_data, correction = "isotropic", r = d)$iso,
                               theo = spatstat.explore::Vmark(ref_data, correction = "isotropic", r = d)$theo,
                               sim_m = V_PP))
  r_envV <- GET::global_envelope_test(C_ref_V, type = "rank")

  rComb <- c(d,
             d[1:F_val] + base::max(d),
             d[1:G_val] + base::max(d[1:F_val] + base::max(d)),
             d[1:J_val] + base::max(d[1:G_val] + base::max(d[1:F_val] + base::max(d)))
  )
  K_data <- spatstat.explore::Kest(ref_data,
                                   correction = "isotropic",
                                   r = d)$iso
  F_data <- spatstat.explore::Fest(ref_data,
                                   correction = "rs",
                                   r = d[1:F_val])$rs
  G_data <- spatstat.explore::Gest(ref_data,
                                   correction = "rs",
                                   r = d[1:G_val])$rs
  J_data <- spatstat.explore::Jest(ref_data,
                                   correction = "rs",
                                   r = d[1:J_val])$rs - 1
  Comb_data <- c(sqrt(K_data / pi) - d,
                 F_data,
                 G_data,
                 J_data / J_scale)
  Comb_ref <- GET::create_curve_set(base::list(r = rComb,
                                               obs = Comb_data,
                                               sim_m = base::rbind(sqrt(K_PP / pi) - d,
                                                                   F_PP,
                                                                   G_PP,
                                                                   J_PP / J_scale)
                                               )
                                   )
  r_envComb <-  GET::global_envelope_test(Comb_ref, type = "rank")


  return(base::list(L_env = r_envL,
                    F_env = r_envF,
                    G_env = r_envG,
                    J_env = r_envJ,
                    E_env = r_envE,
                    V_env = r_envV,
                    Comb_env = r_envComb,
                    Sims = sim_list
                    )
         )

}
