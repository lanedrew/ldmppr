# estimate_helpers.R
# Internal helpers for estimate_process_parameters()
# ---------------------------------------------------------------------

#' @noRd
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------

#' @noRd
.validate_common_inputs <- function(x_grid, y_grid, t_grid, upper_bounds, parameter_inits) {
  if (is.null(x_grid) || is.null(y_grid) || is.null(t_grid)) {
    stop("Provide grid values for x_grid, y_grid, and t_grid.", call. = FALSE)
  }
  if (!is.numeric(x_grid) || !is.numeric(y_grid) || !is.numeric(t_grid)) {
    stop("x_grid, y_grid, and t_grid must be numeric vectors.", call. = FALSE)
  }
  if (anyNA(x_grid) || anyNA(y_grid) || anyNA(t_grid)) {
    stop("x_grid, y_grid, and t_grid may not contain NA.", call. = FALSE)
  }

  if (is.null(upper_bounds) || length(upper_bounds) != 3) {
    stop("Provide `upper_bounds = c(b_t, b_x, b_y)`.", call. = FALSE)
  }
  if (anyNA(upper_bounds) || any(!is.finite(upper_bounds))) {
    stop("upper_bounds must be finite numeric values.", call. = FALSE)
  }
  if (upper_bounds[1] < max(t_grid) || upper_bounds[2] < max(x_grid) || upper_bounds[3] < max(y_grid)) {
    stop("Grid values for t, x, or y exceed upper_bounds.", call. = FALSE)
  }

  if (is.null(parameter_inits) || length(parameter_inits) != 8 ||
      anyNA(parameter_inits) || any(!is.finite(parameter_inits))) {
    stop("Provide valid `parameter_inits` (numeric length 8, finite).", call. = FALSE)
  }
  if (any(parameter_inits[2:8] < 0)) {
    stop("Provide valid `parameter_inits`: entries 2:8 must be >= 0.", call. = FALSE)
  }

  invisible(TRUE)
}

# ---------------------------------------------------------------------
# Future plan helpers
# ---------------------------------------------------------------------

#' @noRd
.maybe_set_future_plan <- function(set_future_plan,
                                   will_parallelize,
                                   num_cores,
                                   max_workers = NULL,
                                   verbose = TRUE) {
  if (!isTRUE(set_future_plan) || !isTRUE(will_parallelize)) return(NULL)

  if (!requireNamespace("future", quietly = TRUE)) {
    stop("Package 'future' is required when set_future_plan=TRUE and parallelization is used.", call. = FALSE)
  }

  original_plan <- future::plan()

  num_cores <- max(1L, as.integer(num_cores))
  if (!is.null(max_workers)) num_cores <- min(num_cores, as.integer(max_workers))

  future::plan(future::multisession, workers = num_cores)

  if (isTRUE(verbose)) message("future plan set to multisession with workers = ", num_cores)

  original_plan
}

# ---------------------------------------------------------------------
# Grid schedule helpers
# ---------------------------------------------------------------------

#' @noRd
.make_grid_level <- function(level, upper_bounds) {
  # level can be c(nx,ny,nt) or list(nx=,ny=,nt=)
  if (is.numeric(level) && length(level) == 3) {
    nx <- as.integer(level[1]); ny <- as.integer(level[2]); nt <- as.integer(level[3])
  } else if (is.list(level) && all(c("nx", "ny", "nt") %in% names(level))) {
    nx <- as.integer(level$nx); ny <- as.integer(level$ny); nt <- as.integer(level$nt)
  } else {
    stop("grid_levels entries must be c(nx,ny,nt) or list(nx=,ny=,nt=).", call. = FALSE)
  }

  if (nx < 2L || ny < 2L || nt < 2L) {
    stop("Each grid level must have nx, ny, nt >= 2.", call. = FALSE)
  }

  list(
    x = seq(0, upper_bounds[2], length.out = nx),
    y = seq(0, upper_bounds[3], length.out = ny),
    t = seq(0, upper_bounds[1], length.out = nt)
  )
}

#' @noRd
.make_grid_schedule <- function(x_grid, y_grid, t_grid, upper_bounds, grid_levels = NULL) {
  if (is.null(grid_levels)) {
    return(list(list(x = x_grid, y = y_grid, t = t_grid)))
  }
  if (!is.list(grid_levels) || length(grid_levels) < 1) {
    stop("grid_levels must be NULL or a non-empty list.", call. = FALSE)
  }
  lapply(grid_levels, .make_grid_level, upper_bounds = upper_bounds)
}

# ---------------------------------------------------------------------
# Bounds helpers (finite bounds used by BOBYQA / NEWUOA / PRAXIS)
# ---------------------------------------------------------------------

#' @noRd
.needs_finite_bounds <- function(alg) {
  alg %in% c("NLOPT_LN_BOBYQA", "NLOPT_LN_NEWUOA", "NLOPT_LN_PRAXIS")
}

#' @noRd
.derive_finite_bounds <- function(init,
                                  mult = 25,
                                  min_ub = c(NA_real_, 50, 5, 50, 50, 50, 50, 5)) {
  # init is length 8; parameter 1 may be signed, others nonnegative
  if (length(init) != 8) stop("init must be length 8.", call. = FALSE)

  ub_pos <- pmax(min_ub[-1], abs(init[-1]) * mult, 1e-3)

  # alpha_1 can be negative; give symmetric bounds around 0
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
  if (any(!is.finite(finite_bounds$lb)) || any(!is.finite(finite_bounds$ub))) {
    stop("finite_bounds must be finite numeric values.", call. = FALSE)
  }
  if (any(finite_bounds$ub <= finite_bounds$lb)) {
    stop("finite_bounds must satisfy ub > lb elementwise.", call. = FALSE)
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------
# Data construction for SC likelihood
# ---------------------------------------------------------------------

#' @noRd
.build_sc_matrix <- function(data, delta = NULL) {
  # Returns a (time,x,y) matrix.
  # Attaches attributes for downstream convenience:
  #  - attr(mat,"ldmppr_original") : data.frame (time,x,y,size?) in the order used
  #  - attr(mat,"ldmppr_delta")    : delta used (NA if not applicable)
  #
  # Supported inputs:
  #  - data.frame with (time,x,y) optionally plus size
  #  - data.frame with (x,y,size)   -> requires delta
  #  - matrix with 3 cols:
  #     * if has colnames time,x,y  -> treated as txy
  #     * if has colnames x,y,size  -> requires delta and maps to time
  #     * if no colnames            -> treated as txy (time,x,y) in column order

  # ---- matrix input ----
  if (is.matrix(data)) {
    if (ncol(data) != 3) stop("`data` matrix must have 3 columns.", call. = FALSE)

    cn <- colnames(data)

    if (!is.null(cn) && all(c("x", "y", "size") %in% cn)) {
      if (is.null(delta)) stop("`delta` must be provided when `data` has (x,y,size) but no time.", call. = FALSE)
      o <- order(-data[, "size"])
      size <- as.numeric(data[o, "size"])
      x <- as.numeric(data[o, "x"])
      y <- as.numeric(data[o, "y"])
      time <- power_law_mapping(size, delta)
      mat <- cbind(time = time, x = x, y = y)
      attr(mat, "ldmppr_original") <- data.frame(time = time, x = x, y = y, size = size)
      attr(mat, "ldmppr_delta") <- delta
      return(mat)
    }

    if (!is.null(cn) && all(c("time", "x", "y") %in% cn)) {
      mat <- cbind(
        time = as.numeric(data[, "time"]),
        x = as.numeric(data[, "x"]),
        y = as.numeric(data[, "y"])
      )
      attr(mat, "ldmppr_original") <- as.data.frame(mat)
      attr(mat, "ldmppr_delta") <- NA_real_
      return(mat)
    }

    # no usable colnames => assume already (time,x,y) in order
    mat <- cbind(time = as.numeric(data[, 1]), x = as.numeric(data[, 2]), y = as.numeric(data[, 3]))
    attr(mat, "ldmppr_original") <- as.data.frame(mat)
    attr(mat, "ldmppr_delta") <- NA_real_
    return(mat)
  }

  # ---- data.frame input ----
  if (!is.data.frame(data)) stop("`data` must be a data.frame or matrix.", call. = FALSE)

  nms <- names(data)

  # Already has time,x,y (+ maybe size)
  if (all(c("time", "x", "y") %in% nms)) {
    mat <- as.matrix(data[, c("time", "x", "y")])
    storage.mode(mat) <- "double"

    # Preserve size if present
    orig <- data[, intersect(c("time", "x", "y", "size"), nms), drop = FALSE]
    attr(mat, "ldmppr_original") <- orig
    attr(mat, "ldmppr_delta") <- NA_real_
    return(mat)
  }

  # Has x,y,size => map to time via delta
  if (all(c("x", "y", "size") %in% nms)) {
    if (is.null(delta)) stop("`delta` must be provided when `data` contains (x,y,size) but no time.", call. = FALSE)

    o <- order(-data$size)
    size <- as.numeric(data$size[o])
    x <- as.numeric(data$x[o])
    y <- as.numeric(data$y[o])
    time <- power_law_mapping(size, delta)

    mat <- cbind(time = time, x = x, y = y)
    attr(mat, "ldmppr_original") <- data.frame(time = time, x = x, y = y, size = size)
    attr(mat, "ldmppr_delta") <- delta
    return(mat)
  }

  stop("`data` must contain either (time,x,y) or (x,y,size).", call. = FALSE)
}

# ---------------------------------------------------------------------
# Likelihood objective + optimizer wrappers
# ---------------------------------------------------------------------

#' @noRd
.make_objective_sc <- function(xg, yg, tg, data_mat, upper_bounds) {
  force(xg); force(yg); force(tg); force(data_mat); force(upper_bounds)

  function(parameters) {
    # full_sc_lhood_fast returns log-likelihood
    likeli <- full_sc_lhood_fast(
      xgrid = xg,
      ygrid = yg,
      tgrid = tg,
      tobs  = data_mat[, 1],
      data  = data_mat,
      params = parameters,
      bounds = upper_bounds
    )
    -likeli
  }
}

#' @noRd
.run_nloptr <- function(x0, obj_fun, algorithm, opts, print_level = 0L, lb = NULL, ub = NULL) {
  if (!requireNamespace("nloptr", quietly = TRUE)) {
    stop("Package 'nloptr' is required.", call. = FALSE)
  }
  if (is.null(lb)) lb <- c(-Inf, rep(0, 7))
  if (is.null(ub)) ub <- rep(Inf, 8)

  opts_full <- utils::modifyList(
    list(algorithm = algorithm, print_level = as.integer(print_level)),
    opts %||% list()
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
# Multi-start jitter helper
# ---------------------------------------------------------------------

#' @noRd
.jitter_inits <- function(base_init, n, sd, seed = NULL) {
  # NOTE:
  # - If seed is non-NULL, jittering is deterministic across calls.
  # - If seed is NULL, jittering uses the current RNG state and will vary across calls.
  n <- as.integer(n)
  if (n <= 1L) return(list(base_init))

  if (!is.null(seed)) set.seed(as.integer(seed))

  out <- vector("list", n)
  out[[1]] <- base_init

  for (k in 2:n) {
    p <- base_init

    # jitter alpha_1 additively (signed)
    p[1] <- p[1] + stats::rnorm(1, sd = sd * max(1, abs(p[1])))

    # jitter positive parameters multiplicatively on log-scale
    pos0 <- pmax(p[-1], 1e-6)
    posj <- pos0 * exp(stats::rnorm(length(pos0), sd = sd))
    p[-1] <- pmax(0, posj)

    out[[k]] <- p
  }

  out
}

# ---------------------------------------------------------------------
# One-level fitting routine (global -> local; optional multistart; optional parallel)
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
                              do_global = FALSE,
                              do_multistart = FALSE,
                              n_starts = 1L,
                              jitter_sd = 0.35,
                              seed = 1L,
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

  inits <- if (isTRUE(do_multistart) && as.integer(n_starts) > 1L) {
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
# Delta search coarse fit helper
# ---------------------------------------------------------------------

#' @noRd
.fit_delta_coarse_sc <- function(delta,
                                 data,
                                 parameter_inits,
                                 coarse_grids,
                                 upper_bounds,
                                 global_algorithm,
                                 local_algorithm,
                                 global_options,
                                 local_options,
                                 finite_bounds,
                                 do_global_coarse = TRUE) {
  mat <- .build_sc_matrix(data, delta = delta)

  res <- .fit_one_level_sc(
    level_grids = coarse_grids,
    data_mat = mat,
    start_params = parameter_inits,
    upper_bounds = upper_bounds,
    global_algorithm = global_algorithm,
    local_algorithm = local_algorithm,
    global_options = global_options,
    local_options = local_options,
    finite_bounds = finite_bounds,
    do_global = isTRUE(do_global_coarse),
    do_multistart = FALSE,
    n_starts = 1L,
    jitter_sd = 0,
    seed = 1L,
    worker_parallel = FALSE,
    worker_verbose = FALSE
  )

  res$best$ldmppr_delta <- delta
  res
}
