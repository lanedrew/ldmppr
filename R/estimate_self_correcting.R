#' estimate parameters of self-correcting model using log-likelihood maximization
#'
#' @param xgrid a vector of grid values for x
#' @param ygrid a vector of gridyvalues for y
#' @param tgrid a vector of grid values for t
#' @param data a matrix of times and locations
#' @param parameter_inits a vector of parameter initialization values
#' @param bounds a vector of bounds for time, x, and y
#' @param opt_algorithm the NLopt algorithm to use for maximization
#' @param verbose `TRUE` or `FALSE` indicating whether to show progress of optimization
#'
#' @return an nloptr object with details of the optimization including solution
#' @export
#'
estimate_parameters_sc <- function(xgrid, ygrid, tgrid, data, parameter_inits,
                                   bounds, opt_algorithm = "NLOPT_LN_NELDERMEAD", verbose = TRUE)  {

  opt_likeli <- function(parameters) {
    likeli <- full_sc_lhood(xgrid = xgrid, ygrid = ygrid, tgrid = tgrid, tobs = data[,1],
                            data = data, params = parameters, bounds = bounds)
    return(-likeli)
  }

  if(verbose){
    print_opts <- 2
  } else{
    print_opts <- 0
  }

  ptm <- base::proc.time()
  parameter_estimates <- nloptr::nloptr(x0 = parameter_inits,
                                        eval_f = opt_likeli,
                                        lb = c(-Inf, 0, 0, 0, 0, 0, 0, 0),
                                        ub = c(Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf),
                                        opts = list("algorithm" = opt_algorithm,
                                                    "maxeval" = 300,
                                                    "xtol_rel" = 1e-2,
                                                    "print_level" = print_opts))
  if(verbose){
    print("Total time to run:")
    print(base::proc.time() - ptm)
  }

  return(parameter_estimates)
}


