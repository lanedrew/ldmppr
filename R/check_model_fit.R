#' check fit of estimated model on the reference data
#'
#' @param ref_data a ppp object for the reference data
#' @param Tmin minimum value for T (defaults to 0)
#' @param Tmax maximum value for T
#' @param params vector of parameter estimates
#' @param M_n point to condition on
#' @param n.sim number of simulated datasets to generate
#' @param verbose `TRUE` or `FALSE` indicating whether to show progress
#' @param mark_model a trained mark model object
#' @param raster_list a list of raster objects
#'
#' @return a list of model fit summaries
#' @export
#'
check_model_fit <- function(ref_data, Tmin = 0, Tmax, params, mark_model, raster_list, M_n, n.sim = 2500, verbose = TRUE){

  d <- spatstat.explore::Kest(spatstat.geom::unmark(ref_data))$r
  dF <- d[1:200]
  dG <- d[1:200]
  dJ <- d[1:200]
  dE <- d
  dV <- d
  sim_window <- spatstat.geom::as.owin(ref_data)
  K_PP <- base::matrix(0, nrow = base::length(d), ncol = n.sim)
  F_PP <- base::matrix(0, nrow = base::length(dF), ncol = n.sim)
  G_PP <- base::matrix(0, nrow = base::length(dG), ncol = n.sim)
  J_PP <- base::matrix(0, nrow = base::length(dJ), ncol = n.sim)
  E_PP <- base::matrix(0, nrow = base::length(dE), ncol = n.sim)
  V_PP <- base::matrix(0, nrow = base::length(dV), ncol = n.sim)
  n_real <- base::numeric()

  if(verbose) print("Beginning simulations!")
  for (j in 1:n.sim){
    sim_j <- sim_spatial_temporal_sc(Tmin, Tmax, params, M_n)[[2]]
    pred_marks_j <- predict_marks(sim_j, raster_list, mark_model)
    PP_xy <- generate_mpp(location_data = sim_j, marks = pred_marks_j, window = sim_window)
    n_real[j] <- base::length(sim_j[,1])
    K_PP[, j] <- spatstat.explore::Kest(spatstat.geom::unmark(PP_xy), correction = "isotropic", r = d)$iso
    F_PP[, j] <- spatstat.explore::Fest(spatstat.geom::unmark(PP_xy), correction = "rs",        r = dF)$rs
    G_PP[, j] <- spatstat.explore::Gest(spatstat.geom::unmark(PP_xy), correction = "rs",        r = dG)$rs
    J_PP[, j] <- spatstat.explore::Jest(spatstat.geom::unmark(PP_xy), correction = "rs",        r = dJ)$rs - 1
    E_PP[, j] <- spatstat.explore::Emark(PP_xy, correction = "isotropic", r = dE)$iso
    V_PP[, j] <- spatstat.explore::Vmark(PP_xy, correction = "isotropic", r = dV)$iso

    if(verbose == TRUE & j %% 100 == 0){
      print(paste0("Simulation iterations complete: ", j, "/", n.sim))
    }

  }
  sim_list <- base::list(Ksim = K_PP,
                         Fsim = F_PP,
                         Gsim = G_PP,
                         Jsim = J_PP,
                         Esim = E_PP,
                         Vsim = V_PP,
                         n_per = n_real)

  if(verbose) print("Simulations complete!")

  K_min <- base::apply(K_PP, 1, base::min)
  K_max <- base::apply(K_PP, 1, base::max)
  K_mean <- base::apply(K_PP, 1, base::mean)
  F_min <- base::apply(F_PP, 1, base::min)
  F_max <- base::apply(F_PP, 1, base::max)
  F_mean <- base::apply(F_PP, 1, base::mean)
  G_min <- base::apply(G_PP, 1, base::min)
  G_max <- base::apply(G_PP, 1, base::max)
  G_mean <- base::apply(G_PP, 1, base::mean)
  J_min <- base::apply(J_PP, 1, base::min)
  J_max <- base::apply(J_PP, 1, base::max)
  J_mean <- base::apply(J_PP, 1, base::mean)
  Jscale <- base::max(J_max)
  E_min <- base::apply(E_PP, 1, base::min)
  E_max <- base::apply(E_PP, 1, base::max)
  E_mean <- base::apply(E_PP, 1, base::mean)
  V_min <- base::apply(V_PP, 1, base::min)
  V_max <- base::apply(V_PP, 1, base::max)
  V_mean <- base::apply(V_PP, 1, base::mean)


  K_ref <- spatstat.explore::Kest(spatstat.geom::unmark(ref_data))
  d <- K_ref$r
  C_ref_L <- GET::create_curve_set(base::list(r = d,
                                              obs = sqrt(K_ref$iso / pi) - d,
                                              theo = sqrt(K_ref$theo / pi) - d ,
                                              sim_m = sqrt(K_PP / pi) - d)
                                   )
  r_envL <- GET::global_envelope_test(C_ref_L, type = "rank")

  C_ref_F <- GET::create_curve_set(base::list(r = dF,
                                              obs = spatstat.explore::Fest(spatstat.geom::unmark(ref_data),
                                                                           r = dF)$rs,
                                              theo = spatstat.explore::Fest(spatstat.geom::unmark(ref_data),
                                                                            r = dF)$theo,
                                              sim_m = F_PP)
                                   )
  r_envF <- GET::global_envelope_test(C_ref_F, type = "rank")

  C_ref_G <- GET::create_curve_set(base::list(r = dG,
                                              obs = spatstat.explore::Gest(spatstat.geom::unmark(ref_data),
                                                                           r = dG)$rs,
                                              theo = spatstat.explore::Gest(spatstat.geom::unmark(ref_data),
                                                                           r = dG)$theo,
                                              sim_m = G_PP)
                                   )

  r_envG <- GET::global_envelope_test(C_ref_G, type = "rank")

  C_ref_J <- GET::create_curve_set(base::list(r = dJ,
                                              obs = spatstat.explore::Jest(spatstat.geom::unmark(ref_data),
                                                                           r = dJ)$rs,
                                              theo = spatstat.explore::Jest(spatstat.geom::unmark(ref_data),
                                                                           r = dJ)$theo,
                                              sim_m = J_PP)
                                   )

  r_envJ <- GET::global_envelope_test(C_ref_J, type = "rank")

  C_ref_E <- GET::create_curve_set(base::list(r = dE,
                               obs = spatstat.explore::Emark(ref_data, correction = "isotropic", r = dE)$iso,
                               theo = spatstat.explore::Emark(ref_data, correction = "isotropic", r = dE)$theo,
                               sim_m = E_PP))
  r_envE <- GET::global_envelope_test(C_ref_E, type = "rank")

  C_ref_V <- GET::create_curve_set(base::list(r = dV,
                               obs = spatstat.explore::Vmark(ref_data, correction = "isotropic", r = dV)$iso,
                               theo = spatstat.explore::Vmark(ref_data, correction = "isotropic", r = dV)$theo,
                               sim_m = V_PP))
  r_envV <- GET::global_envelope_test(C_ref_V, type = "rank")

  rComb <- c(d,
             dF + base::max(d),
             dG + base::max(dF + base::max(d)),
             dJ + base::max(dG + base::max(dF + base::max(d)))
  )
  K_data <- spatstat.explore::Kest(ref_data,
                                   correction = "isotropic",
                                   r = d)$iso
  F_data <- spatstat.explore::Fest(ref_data,
                                   correction = "rs",
                                   r = dF)$rs
  G_data <- spatstat.explore::Gest(ref_data,
                                   correction = "rs",
                                   r = dG)$rs
  J_data <- spatstat.explore::Jest(ref_data,
                                   correction = "rs",
                                   r = dJ)$rs - 1
  Comb_data <- c(sqrt(K_data / pi) - d,
                 F_data,
                 G_data,
                 J_data / Jscale)
  Comb_ref <- GET::create_curve_set(base::list(r = rComb,
                                               obs = Comb_data,
                                               sim_m = base::rbind(sqrt(K_PP / pi) - d,
                                                                   F_PP,
                                                                   G_PP,
                                                                   J_PP / Jscale)
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
