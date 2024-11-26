#' Estimate parameters of the self-correcting model using log-likelihood optimization
#'
#' @description
#' Estimate the parameters of the self-correcting model using the [nloptr::nloptr()] function given a formatted dataset.
#'
#' @param data a matrix or data frame of times and locations in the form (time, x, y).
#' @param x_grid a vector of grid values for x.
#' @param y_grid a vector of grid values for y.
#' @param t_grid a vector of grid values for t.
#' @param parameter_inits a vector of parameter initialization values.
#' @param upper_bounds a vector of upper bounds for time, x, and y.
#' @param opt_algorithm the NLopt algorithm to use for optimization.
#' @param nloptr_options a list of named options for [nloptr::nloptr()] including "maxeval", "xtol_rel", and "maxtime".
#' @param verbose `TRUE` or `FALSE` indicating whether to show progress of optimization.
#'
#' @return an [nloptr::nloptr()] object with details of the optimization including solution.
#' @export
#'
#' @details
#' This function estimates the parameters of the self-correcting model presented in Møller et al. (2016) using the full likelihood.
#' Details regarding the self-correcting model and the estimation procedure can be found in the references.
#'
#' @references
#' Møller, J., Ghorbani, M., & Rubak, E. (2016). Mechanistic spatio-temporal point process models
#' for marked point processes, with a view to forest stand data. \emph{Biometrics}, 72(3), 687–696.
#' \doi{10.1111/biom.12466}.
#'
#' @examples
#' # Load the small example data
#' data(small_example_data)
#' small_example_data <- small_example_data %>%
#'   dplyr::mutate(time = power_law_mapping(size, .5)) %>%
#'   dplyr::select(time, x, y)
#'
#' # Define the grid values
#' x_grid <- seq(0, 25, length.out = 5)
#' y_grid <- seq(0, 25, length.out = 5)
#' t_grid <- seq(0, 1, length.out = 5)
#'
#' # Define the parameter initialization values
#' parameter_inits <- c(1.5, 8.5, .015, 1.5, 3.2, .75, 3, .08)
#'
#' # Define the upper bounds
#' upper_bounds <- c(1, 25, 25)
#'
#' # Estimate the parameters
#' estimate_parameters_sc(
#'   data = small_example_data,
#'   x_grid = x_grid,
#'   y_grid = y_grid,
#'   t_grid = t_grid,
#'   parameter_inits = parameter_inits,
#'   upper_bounds = upper_bounds,
#'   opt_algorithm = "NLOPT_LN_SBPLX",
#'   nloptr_options = list(
#'     maxeval = 25,
#'     xtol_rel = 1e-2
#'   ),
#'   verbose = TRUE
#' )
#'
estimate_parameters_sc <- function(data,
                                   x_grid = NULL,
                                   y_grid = NULL,
                                   t_grid = NULL,
                                   parameter_inits = NULL,
                                   upper_bounds = NULL,
                                   opt_algorithm = "NLOPT_LN_SBPLX",
                                   nloptr_options = list(
                                     maxeval = 400,
                                     xtol_rel = 1e-5,
                                     maxtime = NULL
                                   ),
                                   verbose = TRUE) {
  # Check the arguments
  if (!is.data.frame(data) & !is.matrix(data)) stop("Provide a matrix or data frame of times and locations in the form (time, x, y) for the data argument.", .call = FALSE)
  if (is.data.frame(data)) {
    data <- as.matrix(data)
  }
  if (is.null(x_grid) | is.null(y_grid) | is.null(t_grid)) stop("Provide grid values for the x_grid, y_grid, and t_grid arguments.", .call = FALSE)
  if (length(parameter_inits) != 8 | anyNA(parameter_inits) | any(parameter_inits[2:8] < 0)) stop("Provide valid initialization values for the parameter_inits argument.", .call = FALSE)
  if (is.null(upper_bounds) | !(length(upper_bounds) == 3)) stop("Provide upper bounds for t, x, and y in the form of (b_t, b_x, b_y) for the upper_bounds argument.", .call = FALSE)
  if (upper_bounds[1] < max(t_grid) | upper_bounds[2] < max(x_grid) | upper_bounds[3] < max(y_grid)) stop("Grid values for t, x, or y exceed upper bounds.", .call = FALSE)


  # Define a function to evaluate the self-correcting model full likelihood given a set of parameter values
  opt_likeli <- function(parameters) {
    likeli <- full_sc_lhood(
      xgrid = x_grid,
      ygrid = y_grid,
      tgrid = t_grid,
      tobs = data[, 1],
      data = data,
      params = parameters,
      bounds = upper_bounds
    )
    return(-likeli)
  }

  if (verbose) {
    # Prints values per iteration
    print_opts <- 2
  } else {
    # Prints nothing
    print_opts <- 0
  }

  # Track the time to run
  ptm <- base::proc.time()

  # Perform the optimization of the likelihood
  parameter_estimates <- nloptr::nloptr(
    x0 = parameter_inits,
    eval_f = opt_likeli,
    lb = c(-Inf, 0, 0, 0, 0, 0, 0, 0),
    ub = c(Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf),
    opts = list(
      "algorithm" = opt_algorithm,
      "maxeval" = nloptr_options[["maxeval"]],
      "xtol_rel" = nloptr_options[["xtol_rel"]],
      "maxtime" = nloptr_options[["maxtime"]],
      "print_level" = print_opts
    )
  )
  # Print the time to run if verbose = TRUE
  if (verbose) {
    base::message("Total time to run:")
    print(base::proc.time() - ptm)
  }

  # Return the nloptr object generated by the optimization
  return(parameter_estimates)
}
