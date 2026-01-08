# ---------------------------------------------------------------------
# estimate_helpers.R  (internal helpers for estimate_process_parameters)
# ---------------------------------------------------------------------

#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x

# ---------------------------------------------------------------------
# Data preparation (self-correcting)
# ---------------------------------------------------------------------

#' @noRd
.build_sc_matrix <- function(data, delta = NULL) {
  # Accept data.frame or matrix.
  # Cases:
  #  - already has time,x,y  -> use it
  #  - has x,y,size (+ delta) -> build time via power_law_mapping after sorting desc(size)

  if (is.matrix(data)) {
    if (ncol(data) != 3) stop("`data` matrix must have 3 columns.", call. = FALSE)

    cn <- colnames(data)

    # If matrix has named x,y,size -> map to time
    if (!is.null(cn) && all(c("x", "y", "size") %in% cn)) {
      if (is.null(delta)) stop("`delta` must be provided when `data` has (x,y,size) but no time.", call. = FALSE)
      o <- order(-data[, "size"])
      size <- data[o, "size"]
      x <- data[o, "x"]
      y <- data[o, "y"]
      time <- power_law_mapping(size, delta)
      return(cbind(time = time, x = x, y = y))
    }

    # Otherwise assume it's already (time,x,y)
    return(data)
  }

  if (!is.data.frame(data)) stop("`data` must be a data.frame or matrix.", call. = FALSE)

  nms <- names(data)

  if (all(c("time", "x", "y") %in% nms)) {
    return(as.matrix(data[, c("time", "x", "y")]))
  }

  if (all(c("x", "y", "size") %in% nms)) {
    if (is.null(delta)) stop("`delta` must be provided when `data` contains (x,y,size) but no time.", call. = FALSE)
    o <- order(-data$size)
    size <- data$size[o]
    x <- data$x[o]
    y <- data$y[o]
    time <- power_law_mapping(size, delta)
    return(cbind(time = time, x = x, y = y))
  }

  stop("`data` must contain either (time,x,y) or (x,y,size).", call. = FALSE)
}

# ---------------------------------------------------------------------
# Validation helpers
# ---------------------------------------------------------------------

#' @noRd
.validate_common_inputs <- function(x_grid, y_grid, t_grid, upper_bounds, parameter_inits) {
  if (is.null(x_grid) || is.null(y_grid) || is.null(t_grid)) {
    stop("Provide grid values for x_grid, y_grid, and t_grid.", call. = FALSE)
  }
  if (is.null(upper_bounds) || length(upper_bounds) != 3) {
    stop("Provide `upper_bounds = c(b_t, b_x, b_y)`.", call. = FALSE)
  }
  if (upper_bounds[1] < max(t_grid) || upper_bounds[2] < max(x_grid) || upper_bounds[3] < max(y_grid)) {
    stop("Grid values for t, x, or y exceed upper_bounds.", call. = FALSE)
  }
  if (is.null(parameter_inits) || length(parameter_inits) != 8 ||
      anyNA(parameter_inits) || any(parameter_inits[2:8] < 0)) {
    stop("Provide valid `parameter_inits` (length 8; entries 2:8 >= 0).", call. = FALSE)
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------
# Grid schedule (multi-resolution)
# ---------------------------------------------------------------------

#' @noRd
.make_grid_level <- function(level, upper_bounds) {
  if (is.numeric(level) && length(level) == 3) {
    nx <- as.integer(level[1]); ny <- as.integer(level[2]); nt <- as.integer(level[3])
    return(list(
      x = seq(0, upper_bounds[2], length.out = nx),
      y = seq(0, upper_bounds[3], length.out = ny),
      t = seq(0, upper_bounds[1], length.out = nt)
    ))
  }
  if (is.list(level) && all(c("nx", "ny", "nt") %in% names(level))) {
    nx <- as.integer(level$nx); ny <- as.integer(level$ny); nt <- as.integer(level$nt)
    return(list(
      x = seq(0, upper_bounds[2], length.out = nx),
      y = seq(0, upper_bounds[3], length.out = ny),
      t = seq(0, upper_bounds[1], length.out = nt)
    ))
  }
  stop("grid_levels must be NULL or a list of c(nx,ny,nt) or list(nx=,ny=,nt=) entries.", call. = FALSE)
}

#' @noRd
.make_grid_schedule <- function(x_grid, y_grid, t_grid, upper_bounds, grid_levels = NULL) {
  if (is.null(grid_levels)) {
    return(list(list(x = x_grid, y = y_grid, t = t_grid)))
  }
  lapply(grid_levels, .make_grid_level, upper_bounds = upper_bounds)
}

# ---------------------------------------------------------------------
# Bounds + algorithm helpers
# ---------------------------------------------------------------------

#' @noRd
.needs_finite_bounds <- function(alg) {
  alg %in% c("NLOPT_LN_BOBYQA", "NLOPT_LN_NEWUOA", "NLOPT_LN_PRAXIS")
}

#' @noRd
.derive_finite_bounds <- function(init) {
  # Heuristic "safe-ish" bounds for bound-requiring local optimizers.
  # Parameter 1 can be negative; others constrained >= 0 in your model checks.
  mult <- 25
  min_ub <- c(NA_real_, 50, 5, 50, 50, 50, 50, 5)
  ub_pos <- pmax(min_ub[-1], abs(init[-1]) * mult, 1e-3)
  a1_span <- max(10, abs(init[1]) * mult)

  lb <- c(-a1_span, rep(0, 7))
  ub <- c(a1_span, ub_pos)
  list(lb = lb, ub = ub)
}

#' @noRd
.validate_finite_bounds <- function(finite_bounds) {
  if (!is.list(finite_bounds) || !all(c("lb", "ub") %in% names(finite_bounds))) {
    stop("finite_bounds must be a list with components $lb and $ub.", call. = FALSE)
  }
  if (length(finite_bounds$lb) != 8 || length(finite_bounds$ub) != 8) {
    stop("finite_bounds$lb and finite_bounds$ub must be length 8.", call. = FALSE)
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------
# Objective + optimizer wrappers
# ---------------------------------------------------------------------

#' @noRd
.make_objective_sc <- function(xg, yg, tg, data_mat, upper_bounds) {
  force(xg); force(yg); force(tg); force(data_mat); force(upper_bounds)
  function(parameters) {
    likeli <- full_sc_lhood_fast(
      xgrid  = xg,
      ygrid  = yg,
      tgrid  = tg,
      tobs   = data_mat[, 1],
      data   = data_mat,
      params = parameters,
      bounds = upper_bounds
    )
    -likeli
  }
}

#' @noRd
.run_nloptr <- function(x0, obj_fun, algorithm, opts, print_level = 0L, lb = NULL, ub = NULL) {
  if (is.null(lb)) lb <- c(-Inf, rep(0, 7))
  if (is.null(ub)) ub <- rep(Inf, 8)

  opts_full <- utils::modifyList(
    list(algorithm = algorithm, print_level = as.integer(print_level)),
    opts
  )

  nloptr::nloptr(
    x0 = x0,
    eval_f = obj_fun,
    lb = lb,
    ub = ub,
    opts = opts_full
  )
}

# ---------------------------------------------------------------------
# Multi-start jittering
# ---------------------------------------------------------------------

#' @noRd
.jitter_inits <- function(base_init, n, sd, seed) {
  # Deterministic per-call when seed is fixed (good for reproducibility).
  if (n <= 1L) return(list(base_init))

  if (!is.null(seed)) set.seed(seed)

  out <- vector("list", n)
  out[[1]] <- base_init

  for (k in 2:n) {
    p <- base_init

    # jitter alpha_1 additively (can be negative)
    p[1] <- p[1] + stats::rnorm(1, sd = sd * max(1, abs(p[1])))

    # jitter positive params multiplicatively on log scale
    pos0 <- pmax(p[-1], 1e-6)
    pos_j <- pos0 * exp(stats::rnorm(length(pos0), sd = sd))
    p[-1] <- pmax(0, pos_j)

    out[[k]] <- p
  }

  out
}

# ---------------------------------------------------------------------
# One grid-level fit: optional global stage + local stage (optional multistart)
# ---------------------------------------------------------------------

#' @noRd
.fit_one_level_sc <- function(level_grids,
                              data_mat,
                              start_params,
                              upper_bounds,
                              global_algorithm,
                              local_algorithm,
                              global_options,
                              local_options,
                              finite_bounds,
                              do_global,
                              do_multistart,
                              n_starts,
                              jitter_sd,
                              seed,
                              worker_parallel = FALSE,
                              worker_verbose = FALSE) {
  obj <- .make_objective_sc(level_grids$x, level_grids$y, level_grids$t, data_mat, upper_bounds)
  print_level <- if (isTRUE(worker_verbose)) 2L else 0L

  p0 <- start_params
  global_fit <- NULL

  if (isTRUE(do_global)) {
    global_fit <- .run_nloptr(
      x0 = start_params,
      obj_fun = obj,
      algorithm = global_algorithm,
      opts = global_options,
      print_level = print_level,
      lb = finite_bounds$lb,
      ub = finite_bounds$ub
    )
    p0 <- global_fit$solution
  }

  inits <- if (isTRUE(do_multistart) && n_starts > 1L) {
    .jitter_inits(p0, n = n_starts, sd = jitter_sd, seed = seed)
  } else {
    list(p0)
  }

  if (.needs_finite_bounds(local_algorithm)) {
    lb_loc <- finite_bounds$lb
    ub_loc <- finite_bounds$ub
  } else {
    lb_loc <- c(-Inf, rep(0, 7))
    ub_loc <- rep(Inf, 8)
  }

  run_local_from_init <- function(init) {
    .run_nloptr(
      x0 = init,
      obj_fun = obj,
      algorithm = local_algorithm,
      opts = local_options,
      print_level = 0L,
      lb = lb_loc,
      ub = ub_loc
    )
  }

  if (isTRUE(worker_parallel) && length(inits) > 1L) {
    if (!requireNamespace("furrr", quietly = TRUE) || !requireNamespace("future", quietly = TRUE)) {
      stop("Parallel multi-start requires packages 'future' and 'furrr'.", call. = FALSE)
    }
    local_fits <- furrr::future_map(
      inits,
      run_local_from_init,
      .options = furrr::furrr_options(seed = TRUE)
    )
  } else {
    local_fits <- lapply(inits, run_local_from_init)
  }

  objectives <- vapply(local_fits, function(f) f$objective, numeric(1))
  best_idx <- which.min(objectives)
  best_local <- local_fits[[best_idx]]

  list(best = best_local, local_fits = local_fits, global_fit = global_fit)
}

# ---------------------------------------------------------------------
# Parallel plan helper
# ---------------------------------------------------------------------

#' @noRd
.maybe_set_future_plan <- function(set_future_plan,
                                   will_parallelize,
                                   num_cores,
                                   max_workers,
                                   verbose = FALSE) {
  if (!isTRUE(set_future_plan) || !isTRUE(will_parallelize)) return(invisible(NULL))
  if (!requireNamespace("future", quietly = TRUE)) stop("Package 'future' is required.", call. = FALSE)

  original_plan <- future::plan()

  num_cores <- max(1L, as.integer(num_cores))
  if (!is.null(max_workers)) num_cores <- min(num_cores, as.integer(max_workers))

  future::plan(future::multisession, workers = num_cores)
  if (isTRUE(verbose)) message("future plan set to multisession with workers = ", num_cores)

  original_plan
}
