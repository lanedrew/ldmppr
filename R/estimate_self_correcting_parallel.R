#' estimate parameters of self-correcting model using log-likelihood maximization
#'
#' @param xgrid a vector of grid values for x
#' @param ygrid a vector of grid values for y
#' @param tgrid a vector of grid values for t
#' @param data a dataframe of locations and sizes labelled (x, y, size)
#' @param delta_vals a vector of delta values
#' @param parameter_inits a vector of parameter initialization values
#' @param bounds a vector of bounds for time, x, and y
#' @param opt_algorithm the NLopt algorithm to use for maximization
#' @param verbose `TRUE` or `FALSE` indicating whether to show progress of optimization
#' @param num_cores number of cores to use in parallel estimation
#' @param set_future_plan `TRUE` or `FALSE` indicating whether to change the parallelization plan for use in the function
#'
#' @return a list of the optimal found values and the full results
#' @export
#'
estimate_parameters_parallel <- function(xgrid, ygrid, tgrid, data,
                                         delta_vals = c(.5, 1, 1.5),
                                         parameter_inits,
                                         bounds,
                                         opt_algorithm = "NLOPT_LN_NELDERMEAD",
                                         verbose = TRUE,
                                         num_cores = parallel::detectCores() - 1,
                                         set_future_plan = FALSE) {



  # Save the original plan and restore it on exit
  if (set_future_plan) {
    original_plan <- future::plan()
    base::on.exit(future::plan(original_plan), add = TRUE)

    # Set up a new plan for parallel execution
    future::plan("future::multisession", workers = num_cores)
  }

  # Create storage for possible datasets
  generated_datasets <- list()
  for(i in 1:base::length(delta_vals)){
    # Create generated datasets for different mappings using the power-law function
    generated_datasets[i] <- data %>%
      dplyr::arrange(dplyr::desc(size)) %>%
      dplyr::mutate(time = power_law_mapping(sizes = size, delta = delta_vals[i])) %>%
      dplyr::select(-size) %>%
      base::as.matrix()
  }


  # Define a wrapper function to call in parallel
  parallel_function <- function(data) {
    estimate_parameters_sc(xgrid = xgrid, ygrid = ygrid, tgrid = tgrid,
                           data = data, parameter_inits = parameter_inits,
                           bounds = bounds, opt_algorithm = opt_algorithm, verbose = verbose)
  }

  # Apply the parallel_function across the different datasets
  results <- furrr::future_map(generated_datasets, parallel_function)

  # Determine the optimal mapping and generate a list to return
  min_objective <- base::which.min(base::sapply(results, function(x) x$objective))
  optimal_delta <- delta_vals[min_objective]
  optimal_params <- results[[min_objective]]$solution
  optimal_exit_status <- results[[min_objective]]$status

  optimal_results <- list(optimal_delta = optimal_delta,
                          optimal_params = optimal_params,
                          optimal_status = optimal_exit_status,
                          full_results = results)

  return(optimal_results)
}
