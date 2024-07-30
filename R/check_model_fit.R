#' check fit of estimated model on the reference data
#'
#' @param ref_data a ppp object for the reference data
#' @param Tmin minimum value for T (defaults to 0)
#' @param Tmax maximum value for T
#' @param params vector of parameter estimates
#' @param M_n point to condition on
#' @param n.sim number of simulated datasets to generate
#' @param verbose TRUE or FALSE indicating whether to show progress
#'
#' @return a list of model fit summaries
#' @export
#'
check_model_fit <- function(ref_data, Tmin = 0, Tmax, params, M_n, n.sim = 2500, verbose = TRUE){

  d <- spatstat.explore::Kest(spatstat.geom::unmark(ref_data))$r
  dF <- d[1:200]
  dG <- d[1:200]
  dJ <- d[1:200]
  sim_window <- spatstat.geom::as.owin(ref_data)
  K_PP <- base::matrix(0, nrow = length(d), ncol = n.sim)
  F_PP <- base::matrix(0, nrow = length(dF), ncol = n.sim)
  G_PP <- base::matrix(0, nrow = length(dG), ncol = n.sim)
  J_PP <- base::matrix(0, nrow = length(dJ), ncol = n.sim)
  n.real <- numeric()

  if(verbose) print("Beginning simulations!")
  for (j in 1:n.sim){
    sim_j <- sim_spatial_temporal_sc(Tmin, Tmax, params, M_n)[[2]]
    PP_xy <- spatstat.geom::ppp(sim_j[, 2], sim_j[, 3], window = sim_window)
    n.real[j] <- length(sim_j[,1])
    K_PP[, j] <- spatstat.explore::Kest(PP_xy, correction = "isotropic", r = d)$iso
    F_PP[, j] <- spatstat.explore::Fest(PP_xy, correction = "rs",        r = dF)$rs
    G_PP[, j] <- spatstat.explore::Gest(PP_xy, correction = "rs",        r = dG)$rs
    J_PP[, j] <- spatstat.explore::Jest(PP_xy, correction = "rs",        r = dJ)$rs - 1

    if(verbose == TRUE & j %% 100 == 0){
      print(paste0("Simulation iterations complete: ", j, "/", n.sim))
    }
  }
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
}
