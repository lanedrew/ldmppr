#' Self-correcting temporal model likelihood function
#'
#' @param par a vector of parameters (alpha_1, beta_1, gamma_1)
#' @param evalt the vector of times to evaluate the likelihood at
#' @param obst the observed time to evaluate the likelihood at
#'
#' @return the evaluated likelihood for the temporal component
#' @export
#'
TempSpruce <- function(par, evalt, obst) {
  alpha1 = par[1];
  beta1 =  par[2];
  gamma1 = par[3];
  Nrow = base::length(obst);
  oneDist = rdistC(evalt, obst);
  N_t = CondSumCpp(obst, evalt, rep(1, Nrow));
  result = base::exp(alpha1 + beta1 * evalt - gamma1 * N_t);
  return(result)
}


#' Simulate the temporal component of the self-correcting model
#'
#' @param Tmin minimum time value
#' @param Tmax maximum time value
#' @param par a vector of parameters (alpha_1, beta_1, gamma_1)
#'
#' @return a list of thinned and unthinned temporal samples
#' @export
#'
Sim_Temporal <- function(Tmin, Tmax, par){
  alpha1 = par[1];
  beta1  = par[2];
  gamma1 = par[3];
  Lmax <- base::exp( alpha1 + beta1*Tmax )
  N_max <- Lmax*Tmax
  sample_size <- stats::rpois( 1, N_max )
  sample <- base::sort( stats::runif( sample_size, Tmin, Tmax ) )
  hist <- base::matrix(0, ncol=1)
  U  <- stats::runif(sample_size, 0, 1)
  t_sim <- base::numeric(sample_size)
  lambda_star <- base::numeric()
  for (i in 1:sample_size){
    lambda_star[i] <- TempSpruce(par, sample[i], hist)
    if(U[i] < (lambda_star[i]/Lmax)){
      t_sim[i] <- c( sample[i])
      hist <- base::rbind(hist, sample[i])
    }
    else {t_sim[i] <- NA
    hist <- hist
    }
  }
  return (base::list(unthin=t_sim, all=sample))
}


#' Simulate the spatial component of the self-correcting model
#'
#' @param M_n vector of (x,y)-coordinates for largest point
#' @param par a vector of parameters (alpha_2, beta_2)
#' @param nsim.t number of points to simulate
#'
#' @return a matrix of point locations in the (x,y)-plane
#' @export
#'
Sim_spatial <- function(M_n, par, nsim.t){
  Loc <- base::matrix(M_n, ncol=2)  # x and y corresponding to M_n(t_0)
  lengLoc <- 1
  while(lengLoc < nsim.t){
    newp <- c(stats::runif(1, 0, 100), stats::runif(1, 0, 100))
    if( stats::runif(1, 0, 1) <= interactionCpp( Loc, newp, par ) ){
      Loc <- base::rbind(Loc,newp)
      lengLoc <-  lengLoc + 1
    }
  }
  return(Loc)
}

#' Simulate from the self-correcting model
#'
#' @param Tmin minimum value for T
#' @param Tmax maximum value for T
#' @param par vector of parameter estimates
#' @param M_n point to condition on
#'
#' @return a list of simulated values both thinned and unthinned
#' @export
#'
Sim_spatio_temp <- function(Tmin, Tmax, par, M_n){
  Sim_time = stats::na.omit(Sim_Temporal(Tmin, Tmax, par[1:3])$unthin)
  sim_loc =  Sim_spatial(M_n, par[4:5], length(Sim_time))
  Sim_time[1] = 0
  txy_sim = base::cbind(Sim_time, sim_loc)
  thin_vals = (stats::runif(base::nrow(txy_sim), 0, 1) < interactionCpp_st(txy_sim, par[6:8]))
  txy_sim_thin = txy_sim[thin_vals,]
  sim_df <- base::data.frame(time = txy_sim[,1], x = txy_sim[,2], y = txy_sim[,3])
  sim_thin_df <- base::data.frame(time = txy_sim_thin[,1], x = txy_sim_thin[,2], y = txy_sim_thin[,3])
  return(base::list(unthinned = sim_df, thinned = sim_thin_df))
}

#' Predict values from the mark distribution
#'
#' @param sim_realization a realization of simulated values
#' @param raster_list list of rasters
#' @param size_model a predictive model
#'
#' @return a vector of predictions
#' @export
#'
pred_marks <- function(sim_realization, raster_list, size_model){

  s <- sim_realization[,c("x", "y")]
  covars <- base::do.call(base::cbind, base::lapply(raster_list, function(x) terra::extract(x, s, method = "bilinear")))
  X <- covars
  X$x <- s[,1]
  X$y <- s[,2]
  X$age <- sim_realization[,1]

  return(stats::predict(size_model, X))
}


#' Generate a marked process given locations and marks
#'
#' @param location_data a set of locations
#' @param marks a vector of marks
#' @param window a vector of window bounds (2 for x, 2 for y)
#'
#' @return a ppp object with marks
#' @export
#'
generate_mpp <- function(location_data, marks, window){

  location_data <- base::as.data.frame(location_data)
  marked_pp <- spatstat.geom::ppp(location_data[,c("x")],
                                  location_data[,c("y")],
                                  window[1:2], window[3:4],
                                  marks = marks)

  return(marked_pp)
}
