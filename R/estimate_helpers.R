# estimate_helpers.R
# Internal helpers for estimate_process_parameters()
# ---------------------------------------------------------------------

#' @noRd
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------

#' @noRd
.validate_epp_inputs <- function(data,
                                 grids,
                                 budgets,
                                 parameter_inits,
                                 delta,
                                 strategy,
                                 global_algorithm,
                                 local_algorithm,
                                 starts) {

  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("`data` must be a data.frame or matrix.", call. = FALSE)
  }
  if (!is_ldmppr_grids(grids)) stop("`grids` must be an ldmppr_grids object.", call. = FALSE)
  if (!is_ldmppr_budgets(budgets)) stop("`budgets` must be an ldmppr_budgets object.", call. = FALSE)

  ub <- grids$upper_bounds
  if (is.null(ub) || length(ub) != 3L || anyNA(ub) || any(!is.finite(ub))) {
    stop("`grids$upper_bounds` must be finite numeric length 3: c(b_t, b_x, b_y).", call. = FALSE)
  }

  # data shape checks
  has_time <- FALSE
  has_size <- FALSE

  if (is.data.frame(data)) {
    has_time <- all(c("time", "x", "y") %in% names(data))
    has_size <- all(c("x", "y", "size") %in% names(data))
  } else {
    cn <- colnames(data)
    has_time <- !is.null(cn) && all(c("time", "x", "y") %in% cn)
    has_size <- !is.null(cn) && all(c("x", "y", "size") %in% cn)
    if (is.null(cn) && ncol(data) == 3L) has_time <- TRUE # legacy
  }

  if (!has_time && !has_size) {
    stop("`data` must contain either (time,x,y) or (x,y,size).", call. = FALSE)
  }

  # delta rules
  if (has_time && !is.null(delta) && length(delta) > 1L) {
    stop("Delta search (length(delta)>1) is not valid when `data` already contains time.", call. = FALSE)
  }
  if (!has_time && has_size && is.null(delta)) {
    stop("`data` has (x,y,size) but no time. Provide `delta`.", call. = FALSE)
  }

  # parameter_inits if provided
  if (!is.null(parameter_inits)) {
    if (!is.numeric(parameter_inits) || length(parameter_inits) != 8L ||
        anyNA(parameter_inits) || any(!is.finite(parameter_inits))) {
      stop("`parameter_inits` must be numeric length 8, finite, non-NA.", call. = FALSE)
    }
    if (any(parameter_inits[2:8] < 0)) stop("`parameter_inits[2:8]` must be >= 0.", call. = FALSE)
  }

  # strategy vs grids
  if (strategy == "global_local" && length(grids$levels) != 1L) {
    stop("strategy='global_local' requires a single grid level in `grids`.", call. = FALSE)
  }

  # starts normalization
  .normalize_epp_starts(starts)

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


# -------------------------------------------------------------------------------------
# Bounds helpers (finite bounds used by global algorithms and BOBYQA / NEWUOA / PRAXIS)
# -------------------------------------------------------------------------------------

#' @noRd
.needs_finite_bounds <- function(alg) {
  # NLopt docs: all global algorithms require bound constraints
  # (and several local methods behave best / require bounds too).
  is_global <- is.character(alg) && length(alg) == 1L && grepl("^NLOPT_G", alg)

  is_local_bound_req <- alg %in% c(
    "NLOPT_LN_BOBYQA",
    "NLOPT_LN_NEWUOA",
    "NLOPT_LN_PRAXIS"
  )

  is_global || is_local_bound_req
}

#' @noRd
.derive_finite_bounds <- function(init,
                                  mult = 15,
                                  min_ub = c(NA_real_, 25, 2, 25, 10, 25, 25, 2),
                                  max_ub = c(NA_real_, 500, 50, 500, 50, 500, 500, 50)) {
  if (length(init) != 8) stop("init must be length 8.", call. = FALSE)

  ub_pos <- pmax(min_ub[-1], abs(init[-1]) * mult, 1e-3)
  ub_pos <- pmin(ub_pos, max_ub[-1])

  a1_span <- max(10, abs(init[1]) * mult)
  lb <- c(-a1_span, rep(0, 7))
  ub <- c( a1_span, ub_pos)

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


#' @noRd
.clamp_to_bounds <- function(x, lb, ub, eps = 1e-8) {
  stopifnot(length(x) == length(lb), length(x) == length(ub))

  # compute a per-parameter epsilon that won't invert bounds on tiny spans
  span <- ub - lb
  eps_i <- rep(eps, length(x))
  ok_span <- is.finite(span) & span > 0
  eps_i[ok_span] <- pmin(eps, span[ok_span] * 1e-10)

  x2 <- x

  # clamp only where bounds are finite
  fin_lb <- is.finite(lb)
  fin_ub <- is.finite(ub)

  x2[fin_lb] <- pmax(x2[fin_lb], lb[fin_lb] + eps_i[fin_lb])
  x2[fin_ub] <- pmin(x2[fin_ub], ub[fin_ub] - eps_i[fin_ub])

  x2
}

#' @noRd
.ensure_x0_in_bounds <- function(x0, lb, ub, verbose = FALSE, context = NULL) {
  bad <- (is.finite(lb) & (x0 < lb)) | (is.finite(ub) & (x0 > ub))
  if (!any(bad)) return(x0)

  x_new <- .clamp_to_bounds(x0, lb, ub)

  if (isTRUE(verbose)) {
    msg <- "Clamped starting parameters to satisfy bounds"
    if (!is.null(context) && nzchar(context)) msg <- paste0(msg, " (", context, ")")
    msg <- paste0(msg, ". Indices: ", paste(which(bad), collapse = ", "))
    message(msg)
  }

  attr(x_new, "ldmppr_was_clamped") <- TRUE
  x_new
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
# Default initialization (SC)
# ---------------------------------------------------------------------

#' @noRd
.default_parameter_inits_sc <- function(data,
                                        upper_bounds,
                                        delta = NULL,
                                        min_gamma1 = 0.01,
                                        min_beta1  = 1e-3,
                                        min_alpha2 = 1e-3,
                                        min_beta2  = 1e-3,
                                        min_alpha3 = 0,
                                        min_beta3  = 0,
                                        min_gamma3 = 0) {
  if (is.null(upper_bounds) || length(upper_bounds) != 3L ||
      anyNA(upper_bounds) || any(!is.finite(upper_bounds))) {
    stop("default_parameter_inits_sc(): `upper_bounds` must be finite numeric length 3: c(b_t, b_x, b_y).",
         call. = FALSE)
  }

  # Use the same construction logic as the likelihood pipeline
  mat <- .build_sc_matrix(data, delta = delta)
  t <- as.numeric(mat[, 1])
  x <- as.numeric(mat[, 2])
  y <- as.numeric(mat[, 3])
  n <- length(t)

  # --- temporal component: exp(alpha1 + beta1 t - gamma1 N(t)) ---
  t_span <- max(t, na.rm = TRUE) - min(t, na.rm = TRUE)
  if (!is.finite(t_span) || t_span <= 0) t_span <- 1

  tbar <- mean(t)
  alpha_base <- log(n / t_span)

  gamma1 <- max(min_gamma1, alpha_base / max(1, n))  # ~ log(n)/n scale
  beta1  <- max(min_beta1, 1e-2)

  # Make log-intensity around mid-run roughly alpha_base
  alpha1 <- alpha_base - beta1 * tbar + gamma1 * (n/2)

  # --- spatial inhibition scale: alpha2 ---
  w <- upper_bounds[2]
  h <- upper_bounds[3]
  diag_len <- sqrt(w^2 + h^2)

  alpha2 <- 0.05 * min(w, h)  # fallback

  if (n >= 2L) {
    coords <- cbind(x, y)

    nn <- NULL
    if (requireNamespace("spatstat.geom", quietly = TRUE)) {
      nn <- spatstat.geom::nndist(coords)
      nn <- nn[is.finite(nn)]
    }

    # Fallback only if nndist isn't available for some reason
    if (is.null(nn) || !length(nn)) {
      dmat <- as.matrix(stats::dist(coords))
      diag(dmat) <- Inf
      nn <- apply(dmat, 1L, min, na.rm = TRUE)
      nn <- nn[is.finite(nn)]
    }

    if (length(nn)) {
      nn_med <- stats::median(nn)
      alpha2 <- max(1.5 * nn_med, 0.02 * diag_len)
    }
  }

  alpha2 <- max(min_alpha2, min(alpha2, 0.5 * diag_len))

  # --- spatial interaction shape: beta2 ---
  beta2 <- max(min_beta2, 2)

  # --- spatio-temporal thinning interaction: (alpha3, beta3, gamma3) ---
  # Conservative defaults so thinning doesn't annihilate the process
  alpha3 <- max(min_alpha3, 0.5)
  beta3  <- max(min_beta3, alpha2)
  gamma3 <- max(min_gamma3, 0.05)

  c(alpha1, beta1, gamma1, alpha2, beta2, alpha3, beta3, gamma3)
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

  # Clamp x0 if needed (prevents nloptr from erroring on out-of-bounds starts)
  x0 <- .ensure_x0_in_bounds(
    x0 = x0, lb = lb, ub = ub,
    verbose = (as.integer(print_level) > 0L),
    context = paste0("algorithm=", algorithm)
  )

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
.jitter_inits <- function(base_init, n, sd, seed = NULL, lb = NULL, ub = NULL, clamp = FALSE) {
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

    if (isTRUE(clamp) && !is.null(lb) && !is.null(ub)) {
      p <- .ensure_x0_in_bounds(p, lb = lb, ub = ub, verbose = FALSE)
    }

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
                              global_n_starts = 1L,
                              worker_parallel = FALSE,
                              worker_verbose = FALSE) {

  obj <- .make_objective_sc(level_grids$x, level_grids$y, level_grids$t, data_mat, upper_bounds)
  print_level <- if (isTRUE(worker_verbose)) 2L else 0L

  # ---- bounds for global/local ----
  if (.needs_finite_bounds(global_algorithm)) {
    lb_glob <- finite_bounds$lb
    ub_glob <- finite_bounds$ub
  } else {
    lb_glob <- c(-Inf, rep(0, 7))
    ub_glob <- rep(Inf, 8)
  }

  if (.needs_finite_bounds(local_algorithm)) {
    lb_loc <- finite_bounds$lb
    ub_loc <- finite_bounds$ub
  } else {
    lb_loc <- c(-Inf, rep(0, 7))
    ub_loc <- rep(Inf, 8)
  }

  # ---- ensure x0 is valid for global/local stages ----
  p0 <- .ensure_x0_in_bounds(
    x0 = start_params,
    lb = lb_glob,
    ub = ub_glob,
    verbose = worker_verbose,
    context = "global start"
  )

  global_fit <- NULL
  global_fits <- NULL

  # ---- GLOBAL stage (optional) ----
  if (isTRUE(do_global)) {

    global_n_starts <- max(1L, as.integer(global_n_starts))

    # Multiple *seeded* global restarts are usually more meaningful than different x0 for CRS2
    run_one_global <- function(k) {
      gseed <- as.integer(seed) + as.integer(k) - 1L
      gopts <- utils::modifyList(global_options %||% list(), list(seed = gseed))

      # make sure the x0 is inside bounds (and stable)
      x0g <- .ensure_x0_in_bounds(p0, lb_glob, ub_glob, verbose = FALSE)

      .run_nloptr(
        x0 = x0g,
        obj_fun = obj,
        algorithm = global_algorithm,
        opts = gopts,
        print_level = print_level,
        lb = lb_glob,
        ub = ub_glob
      )
    }

    if (isTRUE(worker_parallel) && global_n_starts > 1L) {
      if (!requireNamespace("furrr", quietly = TRUE) || !requireNamespace("future", quietly = TRUE)) {
        stop("Parallel global restarts require packages 'future' and 'furrr'.", call. = FALSE)
      }
      global_fits <- furrr::future_map(
        seq_len(global_n_starts),
        run_one_global,
        .options = furrr::furrr_options(seed = TRUE)
      )
    } else {
      global_fits <- lapply(seq_len(global_n_starts), run_one_global)
    }

    gobj <- vapply(global_fits, function(f) f$objective, numeric(1))
    best_g <- which.min(gobj)
    global_fit <- global_fits[[best_g]]
    p0 <- global_fit$solution
  }

  # ---- LOCAL multistart stage ----
  inits <- if (isTRUE(do_multistart) && as.integer(n_starts) > 1L) {
    .jitter_inits(p0, n = n_starts, sd = jitter_sd, seed = seed)
  } else {
    list(p0)
  }

  # clamp all local inits to local bounds (this fixes “start outside bounds” errors)
  inits <- lapply(inits, .ensure_x0_in_bounds, lb = lb_loc, ub = ub_loc, verbose = FALSE)

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

  list(
    best = best_local,
    local_fits = local_fits,
    global_fit = global_fit,
    global_fits = global_fits
  )
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
                                 do_global_coarse = TRUE,
                                 global_n_starts = 1L,
                                 jitter_sd = 0.35,
                                 seed = 1L) {

  data_mat <- .build_sc_matrix(data, delta = delta)

  # IMPORTANT: forward global_n_starts (and seed/jitter) down to the level fitter
  lvl_fit <- .fit_one_level_sc(
    level_grids = coarse_grids,
    data_mat = data_mat,
    start_params = parameter_inits,
    upper_bounds = upper_bounds,
    global_algorithm = global_algorithm,
    local_algorithm = local_algorithm,
    global_options = global_options,
    local_options = local_options,
    finite_bounds = finite_bounds,
    do_global = isTRUE(do_global_coarse),
    do_multistart = FALSE,          # coarse stage: keep cheap by default
    jitter_sd = jitter_sd,
    seed = as.integer(seed),
    global_n_starts = as.integer(global_n_starts),
    worker_parallel = FALSE,
    worker_verbose = FALSE
  )

  lvl_fit
}


#' @noRd
.normalize_epp_starts <- function(starts) {
  if (is.null(starts)) starts <- list()
  if (!is.list(starts)) stop("`starts` must be a list.", call. = FALSE)

  out <- list(
    global    = as.integer(starts$global %||% 1L),
    local     = as.integer(starts$local  %||% 1L),
    jitter_sd = as.numeric(starts$jitter_sd %||% 0.35),
    seed      = as.integer(starts$seed %||% 1L)
  )

  if (is.na(out$global) || out$global < 1L) stop("starts$global must be >= 1.", call. = FALSE)
  if (is.na(out$local)  || out$local  < 1L) stop("starts$local must be >= 1.", call. = FALSE)
  if (!is.finite(out$jitter_sd) || out$jitter_sd < 0) stop("starts$jitter_sd must be >= 0.", call. = FALSE)
  if (is.na(out$seed)) stop("starts$seed must be a finite integer.", call. = FALSE)

  out
}
