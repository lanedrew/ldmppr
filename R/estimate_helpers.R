# estimate_helpers.R
# Internal helpers for estimate_process_parameters()
# ---------------------------------------------------------------------

# NOTE:
# Shared internal helpers (e.g., `%||%`, `.build_sc_matrix`) are defined in
# internal_helpers.R to avoid duplicate definitions across helper files.

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
                                  data = NULL,
                                  upper_bounds = NULL,
                                  delta = NULL,
                                  mult = 15,
                                  # original-style fallback caps
                                  min_ub = c(NA_real_, 25, 2, 25, 10, 25, 25, 2),
                                  max_ub = c(NA_real_, 500, 50, 500, 50, 500, 500, 50)) {
  if (length(init) != 8) stop("init must be length 8.", call. = FALSE)

  # ---- defaults: old behavior ----
  ub_pos <- pmax(min_ub[-1], abs(init[-1]) * mult, 1e-3)
  ub_pos <- pmin(ub_pos, max_ub[-1])

  # ---- data-derived caps (if possible) ----
  t_span <- NA_real_
  diag_len <- NA_real_

  if (!is.null(data)) {
    mat <- .build_sc_matrix(data, delta = delta)
    tt <- as.numeric(mat[, 1])
    if (length(tt) >= 2L && all(is.finite(tt))) {
      t_span <- max(tt) - min(tt)
      if (!is.finite(t_span) || t_span <= 0) t_span <- NA_real_
    }
  }

  if (!is.null(upper_bounds) && length(upper_bounds) >= 3L &&
      all(is.finite(upper_bounds[2:3])) && all(upper_bounds[2:3] > 0)) {
    w <- upper_bounds[2]
    h <- upper_bounds[3]
    diag_len <- sqrt(w^2 + h^2)
    if (!is.finite(diag_len) || diag_len <= 0) diag_len <- NA_real_
  }

  # If times were mapped to [0,1], t_span will usually be 1.
  # gamma3 cap = max temporal separation
  if (is.finite(t_span)) {
    ub_pos[7] <- min(ub_pos[7], t_span)      # gamma3 is the 8th param -> ub_pos[7]
  }

  # beta3 cap = max spatial separation
  if (is.finite(diag_len)) {
    ub_pos[6] <- min(ub_pos[6], diag_len)    # beta3 is the 7th param -> ub_pos[6]
  }

  # alpha2 cap = max spatial separation (same as beta3)
  if (is.finite(diag_len)) ub_pos[3] <- min(ub_pos[3], diag_len)

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
                                        min_gamma3 = 0,
                                        beta1_factor_target = 5,   # exp(beta1 * t_span) ~ beta1_factor_target
                                        beta1_cap = 50) {          # hard cap on beta1 for stability
  if (is.null(upper_bounds) || length(upper_bounds) != 3L ||
      anyNA(upper_bounds) || any(!is.finite(upper_bounds))) {
    stop("default_parameter_inits_sc(): `upper_bounds` must be finite numeric length 3: c(b_t, b_x, b_y).",
         call. = FALSE)
  }

  # Build (t,x,y) with the same construction logic used elsewhere
  mat <- .build_sc_matrix(data, delta = delta)
  t <- as.numeric(mat[, 1])
  x <- as.numeric(mat[, 2])
  y <- as.numeric(mat[, 3])
  n <- length(t)

  # ---- basic spans / geometry ----
  t0 <- min(t, na.rm = TRUE)
  t1 <- max(t, na.rm = TRUE)
  t_span <- t1 - t0
  if (!is.finite(t_span) || t_span <= 0) t_span <- 1

  w <- upper_bounds[2]
  h <- upper_bounds[3]
  diag_len <- sqrt(w^2 + h^2)
  if (!is.finite(diag_len) || diag_len <= 0) diag_len <- 1

  tbar <- mean(t, na.rm = TRUE)

  # ---- temporal component: exp(alpha1 + beta1 t - gamma1 N(t)) ----
  # Baseline scale from average rate over the observed time span
  alpha_base <- log(max(n, 1) / t_span)

  # gamma1: keep in a conservative ~log(n)/n regime (avoid explosive growth)
  gamma1 <- max(min_gamma1, alpha_base / max(1, n))

  # initialize beta1 using a crude Poisson regression on binned event counts
  # (fallback to exp(beta1 * t_span) ~ beta1_factor_target if regression is unstable)
  beta1_hat <- NA_real_
  if (n >= 20 && all(is.finite(t))) {
    nbins <- max(4L, min(10L, floor(n / 25L)))
    brks <- unique(as.numeric(stats::quantile(t, probs = seq(0, 1, length.out = nbins + 1L),
                                              na.rm = TRUE, names = FALSE)))
    if (length(brks) >= 3L) {
      bin_id <- cut(t, breaks = brks, include.lowest = TRUE, labels = FALSE)
      # counts + bin midpoints + bin widths
      counts <- as.numeric(tabulate(bin_id, nbins))
      mids <- 0.5 * (brks[-1] + brks[-length(brks)])
      dts  <- pmax(brks[-1] - brks[-length(brks)], 1e-8)

      dat <- data.frame(counts = counts, mid = mids, dt = dts)
      ok <- is.finite(dat$counts) & is.finite(dat$mid) & is.finite(dat$dt) & dat$dt > 0
      dat <- dat[ok, , drop = FALSE]

      if (nrow(dat) >= 4 && sum(dat$counts) > 0) {
        dat$log_dt <- log(dat$dt)
        fit <- try(stats::glm(counts ~ mid + offset(log_dt),
                              family = stats::poisson(), data = dat),
                    silent = TRUE)
        if (!inherits(fit, "try-error")) {
          b <- stats::coef(fit)[["mid"]]
          if (is.finite(b)) beta1_hat <- as.numeric(b)
        }
      }
    }
  }

  # fallback heuristic: modest multiplicative change over the interval
  if (!is.finite(beta1_hat)) {
    beta1_hat <- log(max(beta1_factor_target, 1.5)) / t_span
  }

  # enforce minimum + cap
  beta1 <- max(min_beta1, min(abs(beta1_hat), beta1_cap))

  # choose alpha1 so log-intensity around mid-run is ~ alpha_base
  alpha1 <- alpha_base - beta1 * tbar + gamma1 * (n / 2)

  # ---- spatial inhibition scale: alpha2 ----
  alpha2 <- 0.05 * min(w, h)  # fallback

  if (n >= 2L) {
    coords <- cbind(x, y)

    nn <- NULL
    if (requireNamespace("spatstat.geom", quietly = TRUE)) {
      nn <- spatstat.geom::nndist(coords)
      nn <- nn[is.finite(nn)]
    }

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

  # ---- spatial interaction shape: beta2 ----
  beta2 <- max(min_beta2, 2)

  # ---- spatio-temporal thinning interaction: (alpha3, beta3, gamma3) ----
  # Natural caps:
  #  - beta3 is a spatial distance scale, so cap by max possible distance ~ diag_len
  #  - gamma3 is a temporal distance scale, so cap by max time gap ~ t_span (often 1)
  alpha3 <- max(min_alpha3, 0.5)

  beta3  <- max(min_beta3, min(alpha2, diag_len))
  gamma3 <- max(min_gamma3, min(0.10 * t_span, t_span))

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
                              worker_verbose = FALSE,
                              global_to_local = 1L,
                              bound_tol = 1e-6,
                              bound_penalty = 1e-6) {

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

  # ---- helper: soft-penalize boundary hits (tie-breaker only) ----
  .on_bound <- function(x, lb, ub, tol = bound_tol) {
    any(is.finite(lb) & abs(x - lb) <= tol) || any(is.finite(ub) & abs(x - ub) <= tol)
  }
  .score_with_penalty <- function(fit, lb, ub) {
    val <- fit$objective
    if (isTRUE(.on_bound(fit$solution, lb, ub))) val <- val + bound_penalty
    val
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
    global_to_local <- max(1L, as.integer(global_to_local))

    run_one_global <- function(k) {
      gseed <- as.integer(seed) + as.integer(k) - 1L
      gopts <- utils::modifyList(global_options %||% list(), list(seed = gseed))
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

    gobj <- vapply(global_fits, .score_with_penalty, numeric(1), lb = lb_glob, ub = ub_glob)
    best_g <- which.min(gobj)
    global_fit <- global_fits[[best_g]]
    p0 <- global_fit$solution
  }

  # ---- LOCAL stage inits ----
  inits <- list()

  # (A) include top-k global solutions as local seeds (if requested)
  if (isTRUE(do_global) && length(global_fits)) {
    gobj_raw <- vapply(global_fits, function(f) f$objective, numeric(1))
    ord <- order(gobj_raw)
    k <- min(as.integer(global_to_local), length(ord))
    g_solutions <- lapply(ord[seq_len(k)], function(i) global_fits[[i]]$solution)
    inits <- c(inits, g_solutions)
  } else {
    inits <- c(inits, list(p0))
  }

  # (B) optional multistart jitter around best seed (first element)
  if (isTRUE(do_multistart) && as.integer(n_starts) > 1L) {
    jit <- .jitter_inits(inits[[1]], n = n_starts, sd = jitter_sd, seed = seed)
    inits <- c(inits, jit)
  }

  # clamp all local inits to local bounds
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

  lobj <- vapply(local_fits, .score_with_penalty, numeric(1), lb = lb_loc, ub = ub_loc)
  best_idx <- which.min(lobj)
  best_local <- local_fits[[best_idx]]

  # candidates for rescoring at higher resolution
  candidates <- list(
    local = local_fits,
    global = global_fits
  )

  list(
    best = best_local,
    local_fits = local_fits,
    global_fit = global_fit,
    global_fits = global_fits,
    candidates = candidates
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

  # Forward global restart controls to the one-level fitter.
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


# Rescoring helpers
#' @noRd
.safe_unique_params <- function(par_list, tol = 1e-12) {
  # de-dup by exact-ish numeric equality
  keep <- logical(length(par_list))
  kept <- list()
  for (i in seq_along(par_list)) {
    p <- as.numeric(par_list[[i]])
    if (!length(kept)) {
      keep[i] <- TRUE
      kept[[1]] <- p
    } else {
      d <- vapply(kept, function(q) sqrt(sum((p - q)^2)), numeric(1))
      if (all(d > tol)) {
        keep[i] <- TRUE
        kept[[length(kept) + 1]] <- p
      }
    }
  }
  par_list[keep]
}

#' @noRd
.rel_l2 <- function(a, b) {
  a <- as.numeric(a); b <- as.numeric(b)
  den <- max(1e-12, sqrt(sum(b^2)))
  sqrt(sum((a - b)^2)) / den
}

#' @noRd
.is_on_bound <- function(sol, lb, ub, eps = 1e-8) {
  sol <- as.numeric(sol)
  near_lb <- is.finite(lb) & (sol <= lb + eps)
  near_ub <- is.finite(ub) & (sol >= ub - eps)
  any(near_lb | near_ub)
}

#' @noRd
.nudge_interior <- function(sol, lb, ub, eps = 1e-8) {
  # move any bound-huggers slightly into interior
  sol <- as.numeric(sol)
  out <- sol
  for (k in seq_along(sol)) {
    if (is.finite(lb[k]) && out[k] <= lb[k] + eps) out[k] <- lb[k] + 10 * eps
    if (is.finite(ub[k]) && out[k] >= ub[k] - eps) out[k] <- ub[k] - 10 * eps
  }
  out
}

#' @noRd
.rescore_on_grid <- function(param_list, grids_lvl, data_mat, upper_bounds) {
  obj <- .make_objective_sc(grids_lvl$x, grids_lvl$y, grids_lvl$t, data_mat, upper_bounds)
  vals <- vapply(param_list, function(p) obj(as.numeric(p)), numeric(1))
  ord <- order(vals)
  list(params = param_list[ord], obj = vals[ord])
}

#' @noRd
.run_local_from_starts <- function(starts_list, grids_lvl, data_mat, upper_bounds,
                                   local_algorithm, local_options, finite_bounds) {

  obj <- .make_objective_sc(grids_lvl$x, grids_lvl$y, grids_lvl$t, data_mat, upper_bounds)

  if (.needs_finite_bounds(local_algorithm)) {
    lb <- finite_bounds$lb
    ub <- finite_bounds$ub
  } else {
    lb <- c(-Inf, rep(0, 7))
    ub <- rep(Inf, 8)
  }

  starts_list <- lapply(starts_list, .ensure_x0_in_bounds, lb = lb, ub = ub, verbose = FALSE)

  fits <- lapply(starts_list, function(x0) {
    .run_nloptr(
      x0 = x0,
      obj_fun = obj,
      algorithm = local_algorithm,
      opts = local_options,
      print_level = 0L,
      lb = lb,
      ub = ub
    )
  })

  objs <- vapply(fits, function(f) f$objective, numeric(1))
  best_i <- which.min(objs)
  list(best = fits[[best_i]], fits = fits)
}
