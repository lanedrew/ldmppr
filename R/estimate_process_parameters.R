#' Estimate point process parameters using log-likelihood maximization
#'
#' Estimate spatio-temporal point process parameters by maximizing the (approximate) full log-likelihood
#' using \code{\link[nloptr:nloptr]{nloptr}}. For the self-correcting process, the arrival times must be on \eqn{(0,1)}
#' and can either be supplied directly in \code{data} as \code{time}, or constructed from \code{size} via the
#' gentle-decay (power-law) mapping \code{\link{power_law_mapping}} using \code{delta} (single fit) or
#' \code{delta_values} (delta search).
#'
#' @param data a data.frame or matrix. Must contain either columns \code{(time, x, y)} or \code{(x, y, size)}.
#'   If a matrix is provided for delta search, it must have column names \code{c("x","y","size")}.
#' @param process character string specifying the process model. Currently supports \code{"self_correcting"}.
#' @param x_grid,y_grid,t_grid numeric vectors defining the integration grid for \eqn{(x,y,t)}.
#' @param upper_bounds numeric vector of length 3 giving \code{c(b_t, b_x, b_y)}.
#' @param parameter_inits numeric vector of length 8 giving initialization values for the model parameters.
#' @param delta (optional) numeric scalar used only when \code{data} contains \code{(x,y,size)} but not \code{time}.
#' @param delta_values (optional) numeric vector. If supplied, the function fits the model for each value of \code{delta_values}
#'   (mapping \code{size -> time} via \code{\link{power_law_mapping}}) and returns the best fit (lowest objective).
#' @param parallel logical. If \code{TRUE}, uses furrr/future to parallelize either
#'   (a) over \code{delta_values} (when provided) or (b) over multi-start initializations
#'   (when \code{delta_values} is \code{NULL} and \code{n_starts > 1}).
#' @param num_cores Integer number of workers to use when \code{set_future_plan = TRUE}.
#' @param set_future_plan \code{TRUE} or \code{FALSE}, if \code{TRUE}, temporarily sets \code{future::plan(multisession, workers = num_cores)}
#'   and restores the original plan on exit.
#' @param strategy Character string specifying the estimation strategy:
#'  - \code{"local"}: single-level local optimization from \code{parameter_inits}.
#'  - \code{"global_local"}: single-level global optimization (from \code{parameter_inits}) followed by local polish.
#'  - \code{"multires_global_local"}: multi-resolution fitting over \code{grid_levels} (coarsest level uses global + local; finer levels use local polish only).
#' @param grid_levels (optional) list defining the multi-resolution grid schedule when \code{strategy = "multires_global_local"}.
#' Each entry can be a numeric vector \code{c(nx, ny, nt)} or a list with named entries \code{list(nx=..., ny=..., nt=...)}.
#' If \code{NULL}, uses the supplied \code{(x_grid, y_grid, t_grid)} as a single level.
#' @param refine_best_delta \code{TRUE} or \code{FALSE}, if \code{TRUE} and \code{delta_values} is supplied, performs a final refinement
#'  fit at the best delta found using the full multi-resolution strategy.
#' @param global_algorithm,local_algorithm character strings specifying the NLopt algorithms to use for
#' the global and local optimization stages, respectively.
#' @param global_options,local_options named lists of options to pass to \code{nloptr::nloptr()} for
#' the global and local optimization stages, respectively.
#' @param n_starts integer number of multi-start initializations to use for the local optimization stage.
#' @param jitter_sd numeric standard deviation used to jitter the multi-start initializations.
#' @param seed integer random seed used for multi-start initialization jittering.
#' @param finite_bounds (optional) list with components \code{lb} and \code{ub} giving finite lower and upper bounds
#' for all 8 parameters. Used only when the selected optimization algorithms require finite bounds.
#' @param verbose \code{TRUE} or \code{FALSE}, if \code{TRUE}, prints progress messages during fitting.
#'
#'
#' @details
#' For the self-correcting process, the log-likelihood integral is approximated using the supplied grid
#' \code{(x_grid, y_grid, t_grid)} over the bounded domain \code{upper_bounds}.
#' When \code{delta_values} is supplied, this function performs a grid search over \code{delta} values, fitting the model
#' separately for each mapped dataset and selecting the best objective value.
#'
#' @return an object of class \code{"ldmppr_fit"} containing the best \code{nloptr} fit and (optionally) all fits from a delta search.
#'
#' @references
#' Møller, J., Ghorbani, M., & Rubak, E. (2016). Mechanistic spatio-temporal point process models
#' for marked point processes, with a view to forest stand data. \emph{Biometrics}, 72(3), 687–696.
#' \doi{10.1111/biom.12466}.
#'
#' @examples
#' data(small_example_data)
#'
#' x_grid <- seq(0, 25, length.out = 5)
#' y_grid <- seq(0, 25, length.out = 5)
#' t_grid <- seq(0, 1,  length.out = 5)
#'
#' parameter_inits <- c(1.5, 8.5, .015, 1.5, 3.2, .75, 3, .08)
#' upper_bounds <- c(1, 25, 25)
#'
#' fit <- estimate_process_parameters(
#'   data = small_example_data,
#'   process = "self_correcting",
#'   x_grid = x_grid,
#'   y_grid = y_grid,
#'   t_grid = t_grid,
#'   upper_bounds = upper_bounds,
#'   parameter_inits = parameter_inits,
#'   delta = 1,
#'   strategy = "global_local",
#'   global_algorithm = "NLOPT_GN_CRS2_LM",
#'   local_algorithm = "NLOPT_LN_BOBYQA",
#'   global_options = list(maxeval = 150),
#'   local_options = list(maxeval = 25, xtol_rel = 1e-2),
#'   verbose = TRUE
#' )
#'
#' coef(fit)
#' logLik(fit)
#'
#' \donttest{
#' # Delta-search example (data has x,y,size; time is derived internally for each delta)
#' fit_delta <- estimate_process_parameters(
#'   data = small_example_data, # x,y,size
#'   process = "self_correcting",
#'   x_grid = x_grid,
#'   y_grid = y_grid,
#'   t_grid = t_grid,
#'   upper_bounds = upper_bounds,
#'   parameter_inits = parameter_inits,
#'   delta_values = c(0.35, 0.5, 0.65, 0.9, 1.0),
#'   parallel = TRUE,
#'   set_future_plan = TRUE,
#'   num_cores = 2,
#'   strategy = "multires_global_local",
#'   grid_levels = list(
#'   list(nx = 5, ny = 5, nt = 5),
#'   list(nx = 8, ny = 8, nt = 8),
#'   list(nx = 10, ny = 10, nt = 10)
#'   ),
#'   global_options = list(maxeval = 100),
#'   local_options  = list(maxeval = 100, xtol_rel = 1e-3),
#'   n_starts = 3,
#'   refine_best_delta = TRUE,
#'   verbose = TRUE
#' )
#' plot(fit_delta)
#' }
#'
#' @export
estimate_process_parameters <- function(data,
                                        process = c("self_correcting"),
                                        x_grid = NULL,
                                        y_grid = NULL,
                                        t_grid = NULL,
                                        upper_bounds = NULL,
                                        parameter_inits = NULL,
                                        delta = NULL,
                                        delta_values = NULL,
                                        parallel = FALSE,
                                        num_cores = max(1L, parallel::detectCores() - 1L),
                                        set_future_plan = FALSE,
                                        strategy = c("local", "global_local", "multires_global_local"),
                                        grid_levels = NULL,
                                        refine_best_delta = TRUE,
                                        global_algorithm = "NLOPT_GN_CRS2_LM",
                                        local_algorithm = "NLOPT_LN_BOBYQA",
                                        global_options = list(maxeval = 150),
                                        local_options = list(maxeval = 300, xtol_rel = 1e-5, maxtime = NULL),
                                        n_starts = 1L,
                                        jitter_sd = 0.35,
                                        seed = 1L,
                                        finite_bounds = NULL,
                                        verbose = TRUE) {

  process <- match.arg(process)
  strategy <- match.arg(strategy)

  if (process != "self_correcting") {
    stop("Only process='self_correcting' is currently implemented.", call. = FALSE)
  }

  .validate_common_inputs(x_grid, y_grid, t_grid, upper_bounds, parameter_inits)

  n_starts <- as.integer(n_starts)
  if (n_starts < 1L) stop("n_starts must be >= 1.", call. = FALSE)

  # Will we actually use furrr?
  will_parallelize <- isTRUE(parallel) && (
    (!is.null(delta_values)) ||
      (is.null(delta_values) && n_starts > 1L)
  )

  max_workers <- if (!is.null(delta_values)) length(delta_values) else if (n_starts > 1L) n_starts else NULL
  original_plan <- .maybe_set_future_plan(
    set_future_plan = set_future_plan,
    will_parallelize = will_parallelize,
    num_cores = num_cores,
    max_workers = max_workers,
    verbose = verbose
  )
  if (!is.null(original_plan)) on.exit(future::plan(original_plan), add = TRUE)

  grid_schedule <- .make_grid_schedule(
    x_grid = x_grid, y_grid = y_grid, t_grid = t_grid,
    upper_bounds = upper_bounds,
    grid_levels = grid_levels
  )

  if (is.null(finite_bounds)) {
    finite_bounds <- .derive_finite_bounds(parameter_inits)
  } else {
    .validate_finite_bounds(finite_bounds)
  }

  t_start <- proc.time()

  # -------------------------------------------------------------------
  # MODE A: single fit
  # -------------------------------------------------------------------
  if (is.null(delta_values)) {

    data_mat  <- .build_sc_matrix(data, delta = delta)
    data_orig <- attr(data_mat, "ldmppr_original")
    if (is.null(data_orig)) {
      # Fallback: keep *something* reproducible rather than NULL
      # (train_mark_model may still error later if `size` is absent)
      data_orig <- if (is.data.frame(data)) data else as.data.frame(data_mat)
    }

    do_global_first <- (strategy %in% c("global_local", "multires_global_local"))
    do_multires     <- (strategy == "multires_global_local")

    current_params <- parameter_inits
    stage_store <- list()

    n_levels <- length(grid_schedule)
    for (lvl in seq_len(n_levels)) {
      grids_lvl <- grid_schedule[[lvl]]
      do_global <- do_global_first && (lvl == 1L)
      do_multistart <- (n_starts > 1L)

      lvl_fit <- .fit_one_level_sc(
        level_grids = grids_lvl,
        data_mat = data_mat,
        start_params = current_params,
        upper_bounds = upper_bounds,
        global_algorithm = global_algorithm,
        local_algorithm = local_algorithm,
        global_options = global_options,
        local_options = local_options,
        finite_bounds = finite_bounds,
        do_global = do_global,
        do_multistart = do_multistart,
        n_starts = n_starts,
        jitter_sd = jitter_sd,
        seed = seed,
        worker_parallel = isTRUE(parallel) && (n_starts > 1L),
        worker_verbose = verbose && (lvl == n_levels)
      )

      stage_store[[lvl]] <- lvl_fit
      current_params <- lvl_fit$best$solution
      if (!do_multires) break
    }

    secs <- (proc.time() - t_start)[3]

    final_grid <- list(
      x_grid = grid_schedule[[length(stage_store)]]$x,
      y_grid = grid_schedule[[length(stage_store)]]$y,
      t_grid = grid_schedule[[length(stage_store)]]$t,
      upper_bounds = upper_bounds
    )

    # IMPORTANT: pass data + data_original through the constructor (not after)
    fit_obj <- new_ldmppr_fit(
      process = "self_correcting",
      fit = stage_store[[length(stage_store)]]$best,
      fits = stage_store,
      mapping = list(delta = if (!is.null(delta)) delta else NA_real_),
      grid = final_grid,
      data_summary = list(n = nrow(data_mat)),
      data = data_mat,
      data_original = data_orig,
      engine = "nloptr",
      call = match.call(),
      timing = list(seconds = secs)
    )

    return(fit_obj)
  }

  # -------------------------------------------------------------------
  # MODE B: delta search
  # -------------------------------------------------------------------
  if (is.data.frame(data) && all(c("time", "x", "y") %in% names(data))) {
    stop("`delta_values` is only valid when `data` does NOT contain time (i.e., has x,y,size).", call. = FALSE)
  }
  if (is.matrix(data)) {
    cn <- colnames(data)
    ok <- !is.null(cn) && all(c("x", "y", "size") %in% cn)
    if (!ok) stop("For `delta_values`, matrix `data` must have colnames including 'x', 'y', 'size'.", call. = FALSE)
  }

  if (isTRUE(verbose)) message("Starting delta search with ", length(delta_values), " delta values.")

  coarse_grids <- grid_schedule[[1]]
  do_global_coarse <- (strategy %in% c("global_local", "multires_global_local"))

  fit_delta_coarse <- function(d) {
    .fit_delta_coarse_sc(
      delta = d,
      data = data,
      parameter_inits = parameter_inits,
      coarse_grids = coarse_grids,
      upper_bounds = upper_bounds,
      global_algorithm = global_algorithm,
      local_algorithm = local_algorithm,
      global_options = global_options,
      local_options = local_options,
      finite_bounds = finite_bounds,
      do_global_coarse = do_global_coarse
    )
  }

  if (isTRUE(parallel)) {
    if (!requireNamespace("furrr", quietly = TRUE) || !requireNamespace("future", quietly = TRUE)) {
      stop("Parallel delta search requires packages 'future' and 'furrr'.", call. = FALSE)
    }
    coarse_results <- furrr::future_map(
      as.list(delta_values),
      fit_delta_coarse,
      .options = furrr::furrr_options(seed = TRUE)
    )
  } else {
    coarse_results <- lapply(delta_values, fit_delta_coarse)
  }

  coarse_obj <- vapply(coarse_results, function(r) r$best$objective, numeric(1))
  best_idx <- which.min(coarse_obj)
  best_delta <- delta_values[best_idx]

  if (isTRUE(verbose)) {
    message("Coarse delta search complete. Best delta = ", best_delta,
            " (objective=", signif(coarse_obj[best_idx], 8), ").")
  }

  refined_store <- NULL
  final_fit <- coarse_results[[best_idx]]$best
  final_params <- final_fit$solution

  # Build best transformed data now (and preserve original)
  mat_best  <- .build_sc_matrix(data, delta = best_delta)
  orig_best <- attr(mat_best, "ldmppr_original")
  if (is.null(orig_best)) {
    orig_best <- if (is.data.frame(data)) data else as.data.frame(mat_best)
  }

  if (isTRUE(refine_best_delta) && length(grid_schedule) > 1L) {
    if (isTRUE(verbose)) message("Refining best delta across finer grids...")

    refined_store <- list()
    for (lvl in 2:length(grid_schedule)) {
      grids_lvl <- grid_schedule[[lvl]]
      lvl_fit <- .fit_one_level_sc(
        level_grids = grids_lvl,
        data_mat = mat_best,
        start_params = final_params,
        upper_bounds = upper_bounds,
        global_algorithm = global_algorithm,
        local_algorithm = local_algorithm,
        global_options = global_options,
        local_options = local_options,
        finite_bounds = finite_bounds,
        do_global = FALSE,
        do_multistart = (n_starts > 1L),
        n_starts = n_starts,
        jitter_sd = jitter_sd,
        seed = seed,
        worker_parallel = FALSE,
        worker_verbose = verbose && (lvl == length(grid_schedule))
      )
      refined_store[[lvl]] <- lvl_fit
      final_fit <- lvl_fit$best
      final_params <- lvl_fit$best$solution
    }

    if (isTRUE(verbose)) message("Refinement complete (best delta = ", best_delta, ").")
  }

  secs <- (proc.time() - t_start)[3]

  final_grid <- list(
    x_grid = grid_schedule[[length(grid_schedule)]]$x,
    y_grid = grid_schedule[[length(grid_schedule)]]$y,
    t_grid = grid_schedule[[length(grid_schedule)]]$t,
    upper_bounds = upper_bounds
  )

  # IMPORTANT: pass data + data_original through constructor
  fit_obj <- new_ldmppr_fit(
    process = "self_correcting",
    fit = final_fit,
    fits = list(coarse = coarse_results, refined_best = refined_store),
    mapping = list(
      delta = best_delta,
      delta_values = delta_values,
      chosen_index = best_idx,
      objectives = coarse_obj,
      refined = isTRUE(refine_best_delta) && length(grid_schedule) > 1L
    ),
    grid = final_grid,
    data_summary = list(n = nrow(mat_best)),
    data = mat_best,
    data_original = orig_best,
    engine = "nloptr",
    call = match.call(),
    timing = list(seconds = secs)
  )

  fit_obj
}

