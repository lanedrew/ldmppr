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
Sim.s_t_int <- function(Tmin, Tmax, par, M_n){
  Sim_time = stats::na.omit(Sim_Temporal(Tmin, Tmax, par[1:3])$unthin)
  sim_loc =  Sim_spatial(M_n, par[4:5], length(Sim_time))
  Sim_time[1] = 0
  txy_sim = base::cbind(Sim_time, sim_loc)
  thin_vals = (stats::runif(base::nrow(txy_sim), 0, 1) < interactionCpp_st(txy_sim, par[6:8]))
  txy_sim_thin = txy_sim[thin_vals,]
  return(base::list(txy_sim, txy_sim_thin))
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

  s <- sim_realization[,2:3]
  covars <- base::do.call(base::cbind, base::lapply(raster_list, function(x) terra::extract(x, s, method = "bilinear")))
  X <- covars
  X$x <- s[,1]
  X$y <- s[,2]
  X$age <- sim_realization[,1]

  return(stats::predict(size_model, newdata = X))
}
