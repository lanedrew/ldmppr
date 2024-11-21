#' Estimate parameters of the self-correcting model using log-likelihood maximization in parallel
#'
#' @description
#' Estimate the parameters of the self-correcting model using nloptr given a set of delta values.
#' Makes use of parallel computation to allow the user to identify the optimal delta value more quickly
#' given available computational resources.
#'
#' @param data a data frame of locations and sizes labelled ("x", "y", "size").
#' @param x_grid a vector of grid values for x.
#' @param y_grid a vector of grid values for y.
#' @param t_grid a vector of grid values for t.
#' @param delta_values a vector of delta values.
#' @param parameter_inits a vector of parameter initialization values.
#' @param upper_bounds a vector of upper bounds for time, x, and y.
#' @param opt_algorithm the NLopt algorithm to use for maximization.
#' @param nloptr_options a list of named options for nloptr including "maxeval", "xtol_rel", and "maxtime".
#' @param verbose `TRUE` or `FALSE` indicating whether to show progress of optimization.
#' @param num_cores number of cores to use in parallel estimation.
#' @param set_future_plan `TRUE` or `FALSE` indicating whether to change the parallelization plan for use in the function.
#'
#' @return a list containing the optimal obtained values and the full results.
#' @export
#'
#' @references
#' Møller, J., Ghorbani, M., & Rubak, E. (2016). Mechanistic spatio-temporal point process models
#' for marked point processes, with a view to forest stand data. \emph{Biometrics}, 72(3), 687–696.
#' \doi{10.1111/biom.12466}.
#'
#' @examples
#'
#' # Note: This function is designed to be run in parallel and may be computationally intensive.
#'
#' \dontrun{
#' # Load the small example data
#' data(small_example_data)
#'
#' # Define the grid values
#' x_grid <- seq(0, 25, out.length = 3)
#' y_grid <- seq(0, 25, out.length = 3)
#' t_grid <- seq(0, 1, out.length = 3)
#'
#' # Define the delta values
#' delta_values <- seq(0.75, 1, by = 0.25)
#'
#' # Define the parameter initialization values
#' parameter_inits <- c(1.5, 8.5, .015, 1.5, 3.2, .75, 3, .08)
#'
#' # Define the upper bounds
#' upper_bounds <- c(1, 25, 25)
#'
#' # Estimate the parameters in parallel
#' estimate_parameters_sc_parallel(data = small_example_data,
#'                                 x_grid = x_grid,
#'                                 y_grid = y_grid,
#'                                 t_grid = t_grid,
#'                                 delta_values = delta_values,
#'                                 parameter_inits = parameter_inits,
#'                                 upper_bounds = upper_bounds,
#'                                 opt_algorithm = "NLOPT_LN_SBPLX",
#'                                 nloptr_options = list(maxeval = 50,
#'                                                       xtol_rel = 1e-2),
#'                                 verbose = TRUE,
#'                                 set_future_plan = TRUE)
#' }
#'
estimate_parameters_sc_parallel <- function(data,
                                            x_grid = NULL,
                                            y_grid = NULL,
                                            t_grid = NULL,
                                            delta_values = NULL,
                                            parameter_inits = NULL,
                                            upper_bounds = NULL,
                                            opt_algorithm = "NLOPT_LN_SBPLX",
                                            nloptr_options = list(maxeval = 400,
                                                                  xtol_rel = 1e-5,
                                                                  maxtime = NULL),
                                            verbose = TRUE,
                                            num_cores = parallel::detectCores() - 1,
                                            set_future_plan = FALSE) {

  # Check the arguments
  if(!is.data.frame(data)) stop("Provide a data frame of locations and sizes in the form (x, y, size) for the data argument.", .call = FALSE)
  if(is.null(x_grid) | is.null(y_grid) | is.null(t_grid)) stop("Provide grid values for the x_grid, y_grid, and t_grid arguments.", .call = FALSE)
  if(is.null(delta_values) | any(delta_values < 0)) stop("Provide valid delta values for the size time mapping for the delta_values argument.", .call = FALSE)
  if(length(parameter_inits) != 8 | anyNA(parameter_inits) | any(parameter_inits[2:8] < 0)) stop("Provide valid initialization values for the parameter_inits argument.", .call = FALSE)
  if(is.null(upper_bounds) | !(length(upper_bounds) == 3)) stop("Provide upper bounds for t, x, and y in the form of (b_t, b_x, b_y) for the upper_bounds argument.", .call = FALSE)
  if(upper_bounds[1] < max(t_grid) | upper_bounds[2] < max(x_grid) | upper_bounds[3] < max(y_grid)) stop("Grid values for t, x, or y exceed upper bounds.", .call = FALSE)


  # Save the original plan and restore it on exit
  if (set_future_plan) {
    original_plan <- future::plan()
    base::on.exit(future::plan(original_plan), add = TRUE)

    # Adjust the number of cores if fewer are necessary than the total available
    if(num_cores > length(delta_values)){
      num_cores <- length(delta_values)
    }

    # Set up a new plan for parallel execution
    future::plan("future::multisession", workers = num_cores)
    if(verbose){
      message("Number of workers: ", num_cores)
    }
  }

  # Create storage for possible datasets
  generated_datasets <- list()
  for(i in 1:base::length(delta_values)){

    # Create generated datasets for different mappings using the power-law function
    generated_datasets[[i]] <- data %>%
      dplyr::arrange(dplyr::desc(size)) %>%
      dplyr::mutate(time = power_law_mapping(sizes = size, delta = delta_values[i])) %>%
      dplyr::select(-size) %>%
      dplyr::relocate(time) %>%
      base::as.matrix()
  }

  # Define a wrapper function to call in parallel
  parallel_function <- function(data_vals) {
    estimate_parameters_sc(data = data_vals,
                           x_grid = x_grid,
                           y_grid = y_grid,
                           t_grid = t_grid,
                           parameter_inits = parameter_inits,
                           upper_bounds = upper_bounds,
                           opt_algorithm = opt_algorithm,
                           nloptr_options = nloptr_options,
                           verbose = FALSE)
  }

  # Start timing
  if(verbose) {
    message("Starting parallel optimization.")
    start_time <- base::proc.time()
  }

  # Apply the parallel_function across the different datasets
  results <- furrr::future_map(generated_datasets, parallel_function, .options = furrr::furrr_options(seed = TRUE))

  if (verbose) {
    end_time <- base::proc.time()
    duration <- end_time - start_time
    message("Parallel optimization complete.")
    message("Total time for parallel computation: ", round(duration[3], 2), " seconds")
  }

  # Determine the optimal mapping and generate a list to return
  min_objective <- base::which.min(base::sapply(results, function(x) x$objective))
  optimal_delta <- delta_values[min_objective]
  optimal_params <- results[[min_objective]]$solution
  optimal_exit_status <- results[[min_objective]]$status

  results_list <- list(optimal_results = list(optimal_delta = optimal_delta,
                                              optimal_params = optimal_params,
                                              optimal_status = optimal_exit_status),
                          full_results = results)

  return(results_list)
}
