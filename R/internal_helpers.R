#' Internal helpers (not part of the public API)
#'
#' These functions are used internally by ldmppr and are not intended to be
#' called directly by users.
#'
#' @name ldmppr-internal
#' @keywords internal
NULL


# ---- Constructors for internal S3 classes ----

# Constructor for ldmppr_model_check
#' @rdname ldmppr-internal
#' @keywords internal
new_ldmppr_model_check <- function(combined_env,
                                   envs,
                                   curve_sets,
                                   sim_metrics = NULL,
                                   settings = list(),
                                   call = NULL) {
  structure(
    list(
      combined_env = combined_env,
      envs = envs,               # named list: L, F, G, J, E, V
      curve_sets = curve_sets,   # named list: L, F, G, J, E, V
      sim_metrics = sim_metrics, # optional
      settings = settings,
      call = call
    ),
    class = "ldmppr_model_check"
  )
}


# Constructor for ldmppr_sim
#' @rdname ldmppr-internal
#' @keywords internal
new_ldmppr_sim <- function(process,
                           mpp,
                           realization,
                           params,
                           bounds,
                           anchor_point,
                           thinning,
                           edge_correction,
                           include_comp_inds,
                           competition_radius,
                           call = NULL,
                           meta = list()) {
  structure(
    list(
      process = process,
      mpp = mpp,
      realization = realization,
      params = params,
      bounds = bounds,
      anchor_point = anchor_point,
      thinning = thinning,
      edge_correction = edge_correction,
      include_comp_inds = include_comp_inds,
      competition_radius = competition_radius,
      call = call,
      meta = meta
    ),
    class = c(paste0("ldmppr_sim_", process), "ldmppr_sim")
  )
}


# Constructor for ldmppr_mark_model
#' @rdname ldmppr-internal
#' @keywords internal
new_ldmppr_mark_model <- function(engine,
                                  fit_engine = NULL,
                                  xgb_raw = NULL,
                                  recipe = NULL,
                                  outcome = "size",
                                  feature_names = NULL,
                                  info = list()) {
  structure(
    list(
      engine = engine,            # "xgboost" or "ranger"
      fit_engine = fit_engine,    # xgb.Booster or ranger fit; may be NULL after save
      xgb_raw = xgb_raw,          # raw UBJ for xgboost
      recipe = recipe,            # prepped recipe
      outcome = outcome,          # outcome column name
      feature_names = feature_names,
      info = info,
      cache = new.env(parent = emptyenv())
    ),
    class = "ldmppr_mark_model"
  )
}


# Constructor for ldmppr_fit
#' @rdname ldmppr-internal
#' @keywords internal
new_ldmppr_fit <- function(process,
                           fit,
                           fits = NULL,
                           mapping = NULL,
                           grid = NULL,
                           data_summary = NULL,
                           data = NULL,
                           data_original = NULL,
                           engine = "nloptr",
                           call = NULL,
                           timing = NULL) {

  if (!is.null(data_original) && !is.data.frame(data_original)) {
    stop("`data_original` must be a data.frame when provided.", call. = FALSE)
  }

  structure(
    list(
      process = process,
      engine = engine,
      fit = fit,
      fits = fits,
      mapping = mapping,
      grid = grid,
      data = data,
      data_original = data_original,
      data_summary = data_summary,
      call = call,
      timing = timing
    ),
    class = c(paste0("ldmppr_fit_", process), "ldmppr_fit")
  )
}


# ---- internal helpers ----

# Preprocess new data using the stored recipe and feature names
#' @rdname ldmppr-internal
#' @keywords internal
preprocess_new_data <- function(object, new_data) {
  baked <- if (is.null(object$recipe)) {
    new_data
  } else {
    recipes::bake(object$recipe, new_data = new_data)
  }

  if (!is.null(object$outcome) && object$outcome %in% names(baked)) {
    baked[[object$outcome]] <- NULL
  }

  if (!is.null(object$feature_names)) {
    missing <- setdiff(object$feature_names, names(baked))
    if (length(missing)) {
      stop("New data is missing required predictors: ", paste(missing, collapse = ", "))
    }
    baked <- baked[, object$feature_names, drop = FALSE]
  }

  baked
}


# Rehydrate an xgboost booster from raw data if needed
#' @rdname ldmppr-internal
#' @keywords internal
rehydrate_xgb <- function(object) {
  if (!is.null(object$fit_engine)) return(object$fit_engine)

  if (exists("booster", envir = object$cache, inherits = FALSE)) {
    return(get("booster", envir = object$cache, inherits = FALSE))
  }

  if (is.null(object$xgb_raw)) stop("No xgb_raw found to rehydrate xgboost booster.")
  booster <- xgboost::xgb.load.raw(object$xgb_raw)
  assign("booster", booster, envir = object$cache)
  booster
}


#' Ensure an object is a mark model
#' @rdname ldmppr-internal
#' @keywords internal
as_mark_model <- function(mark_model) {
  if (is.null(mark_model)) stop("Provide a mark model for `mark_model`.", call. = FALSE)

  # Allow passing a path
  if (is.character(mark_model) && length(mark_model) == 1 && file.exists(mark_model)) {
    if (exists("load_mark_model", mode = "function")) {
      mark_model <- load_mark_model(mark_model)
    } else {
      mark_model <- readRDS(mark_model)
    }
  }

  # Legacy bundle
  if (inherits(mark_model, "bundle")) {
    if (!requireNamespace("bundle", quietly = TRUE)) {
      stop("mark_model is a bundled object but the 'bundle' package is not available.", call. = FALSE)
    }
    mark_model <- bundle::unbundle(mark_model)
  }

  mark_model
}


# ---------------------------------------------------------------------
# Helper: build (time, x, y) matrix used by the likelihood
# ---------------------------------------------------------------------
#' @rdname ldmppr-internal
#' @keywords internal
.build_sc_matrix <- function(data, delta = NULL) {
  # Returns numeric matrix with columns: time, x, y
  # Ensures nondecreasing time (required by full_sc_lhood_fast()).

  as_num_mat <- function(df) {
    m <- as.matrix(df)
    storage.mode(m) <- "double"
    m
  }

  # ---- matrix input ----
  if (is.matrix(data)) {
    if (ncol(data) != 3) stop("`data` matrix must have exactly 3 columns.", call. = FALSE)
    cn <- colnames(data)

    # Case A: matrix has x,y,size by name
    if (!is.null(cn) && all(c("x", "y", "size") %in% cn)) {
      if (is.null(delta)) stop("`delta` must be provided when `data` has (x,y,size) but no time.", call. = FALSE)
      o <- order(-data[, "size"], data[, "x"], data[, "y"])
      size <- data[o, "size"]
      x <- data[o, "x"]
      y <- data[o, "y"]
      time <- power_law_mapping(size, delta)
      out <- cbind(time = time, x = x, y = y)
      storage.mode(out) <- "double"
      return(out)
    }

    # Case B: matrix has time,x,y by name but maybe wrong order
    if (!is.null(cn) && all(c("time", "x", "y") %in% cn)) {
      out <- data[, c("time", "x", "y"), drop = FALSE]
      storage.mode(out) <- "double"
      o <- order(out[, "time"], out[, "x"], out[, "y"])
      return(out[o, , drop = FALSE])
    }

    # Case C: unnamed matrix: assume already (time,x,y)
    out <- data
    colnames(out) <- c("time", "x", "y")
    storage.mode(out) <- "double"
    o <- order(out[, "time"], out[, "x"], out[, "y"])
    return(out[o, , drop = FALSE])
  }

  # ---- data.frame input ----
  if (!is.data.frame(data)) stop("`data` must be a data.frame or matrix.", call. = FALSE)

  nms <- names(data)

  # Case 1: already has time,x,y
  if (all(c("time", "x", "y") %in% nms)) {
    out <- as_num_mat(data[, c("time", "x", "y")])
    o <- order(out[, "time"], out[, "x"], out[, "y"])
    return(out[o, , drop = FALSE])
  }

  # Case 2: has x,y,size -> build time via mapping
  if (all(c("x", "y", "size") %in% nms)) {
    if (is.null(delta)) stop("`delta` must be provided when `data` contains (x,y,size) but no time.", call. = FALSE)

    o <- order(-data$size, data$x, data$y)
    size <- data$size[o]
    x <- data$x[o]
    y <- data$y[o]
    time <- power_law_mapping(size, delta)

    out <- cbind(time = time, x = x, y = y)
    storage.mode(out) <- "double"
    return(out)
  }

  stop("`data` must contain either (time,x,y) or (x,y,size).", call. = FALSE)
}

# ---------------------------------------------------------------------
# Helper: default parameter bounds for self-correcting process
# ---------------------------------------------------------------------
#' @rdname ldmppr-internal
#' @keywords internal
.default_sc_param_bounds <- function(txy, upper_bounds) {
  stopifnot(is.matrix(txy), ncol(txy) == 3, length(upper_bounds) == 3)

  bt <- upper_bounds[1]
  bx <- upper_bounds[2]
  by <- upper_bounds[3]

  # domain diagonal as a natural upper scale for distance thresholds
  diag_dom <- sqrt(bx^2 + by^2)

  # also consider observed max distance (often tighter than full diagonal)
  dx <- outer(txy[,2], txy[,2], "-")
  dy <- outer(txy[,3], txy[,3], "-")
  obs_maxdist <- sqrt(max(dx*dx + dy*dy, na.rm = TRUE))
  dmax <- max(1e-6, min(diag_dom, obs_maxdist))

  lb <- c(
    alpha1 = -30,   # baseline log-intensity can be negative
    beta1  = 0,     # your code constrains >=0
    gamma1 = 0,
    alpha2 = 1e-6,  # must be >0 to avoid log issues; allow ~0 for "no spatial interaction"
    beta2  = 0,
    alpha3 = 0,
    beta3  = 0,
    gamma3 = 0
  )

  ub <- c(
    alpha1 = 30,
    beta1  = 50,        # wide; can tighten later if you see extremes
    gamma1 = 50,        # wide; remember N(t) can be up to n
    alpha2 = dmax,      # interaction radius can’t exceed domain/observed scale sensibly
    beta2  = 10,        # keeps pow() stable-ish; widen if needed
    alpha3 = 50,
    beta3  = dmax,      # distance threshold
    gamma3 = bt         # time-lag threshold in [0, bt]
  )

  list(lb = unname(lb), ub = unname(ub))
}


#' Null coalescing operator
#' @rdname ldmppr-internal
#' @keywords internal
`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}


# ----------------------------
# helpers for train_mark_model
# ----------------------------

#' @rdname ldmppr-internal
#' @keywords internal
.require_pkgs <- function(pkgs) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop("Package '", p, "' is required for train_mark_model().", call. = FALSE)
    }
  }
}


#' @rdname ldmppr-internal
#' @keywords internal
.coerce_training_df <- function(x, delta = NULL, xy_bounds = NULL) {
  # x can be data.frame or ldmppr_fit
  fit <- NULL

  delta_provided <- !is.null(delta) && !is.na(delta)

  if (inherits(x, "ldmppr_fit")) {
    fit <- x

    # prefer data_original; fallback to NULL (we can't invent size)
    x <- fit$data_original %||% NULL
    if (is.null(x)) {
      stop("`data` is an ldmppr_fit but `data_original` is NULL; cannot train mark model without size.",
           call. = FALSE)
    }

    # default xy_bounds from fit if not supplied
    if (is.null(xy_bounds)) {
      ub <- fit$grid$upper_bounds %||% NULL
      if (!is.null(ub) && length(ub) == 3) {
        xy_bounds <- c(0, ub[2], 0, ub[3])
      }
    }

    # If user supplied delta, ALWAYS overwrite time using that delta.
    # Otherwise, only derive time if missing, using fit$mapping$delta if available.
    if (!("size" %in% names(x))) {
      stop("Training data must contain column `size`.", call. = FALSE)
    }

    if (delta_provided) {
      x$time <- power_law_mapping(x$size, delta)
    } else if (!("time" %in% names(x))) {
      d_use <- fit$mapping$delta %||% NULL
      if (is.null(d_use) || is.na(d_use)) {
        stop("Training data has no `time`. Provide `delta`, or ensure fit$mapping$delta is present.",
             call. = FALSE)
      }
      x$time <- power_law_mapping(x$size, d_use)
    }
  }

  if (!is.data.frame(x)) {
    stop("Provide a data.frame (x,y,size,time) or an ldmppr_fit for `data`.", call. = FALSE)
  }

  needed1 <- c("x", "y", "size")
  if (!all(needed1 %in% names(x))) {
    stop("Training data must contain columns: x, y, size (and optionally time).", call. = FALSE)
  }

  # If user supplied delta, ALWAYS overwrite time, even if time exists.
  if (delta_provided) {
    x$time <- power_law_mapping(x$size, delta)
  } else if (!("time" %in% names(x))) {
    stop("Training data has no `time`. Provide `delta`.", call. = FALSE)
  }

  list(df = x, xy_bounds = xy_bounds, fit = fit)
}


# ---------------------------
# helpers for check_model_fit
# ---------------------------

#' @rdname ldmppr-internal
#' @keywords internal
infer_xy_bounds_from_ppp <- function(ppp) {
  w <- spatstat.geom::as.owin(ppp)
  xr <- w$xrange
  yr <- w$yrange
  c(xr[1], xr[2], yr[1], yr[2])
}


#' @rdname ldmppr-internal
#' @keywords internal
infer_anchor_from_ppp <- function(ppp) {
  m <- spatstat.geom::marks(ppp)
  if (is.null(m)) {
    return(c(ppp$x[1], ppp$y[1]))
  }
  if (is.numeric(m)) {
    i <- which.max(m)
    return(c(ppp$x[i], ppp$y[i]))
  }
  c(ppp$x[1], ppp$y[1])
}


#' @rdname ldmppr-internal
#' @keywords internal
infer_anchor_from_df <- function(df) {
  if (!all(c("x", "y") %in% names(df))) return(NULL)
  if ("size" %in% names(df) && is.numeric(df$size)) {
    i <- which.max(df$size)
    return(c(df$x[i], df$y[i]))
  }
  c(df$x[1], df$y[1])
}


#' @rdname ldmppr-internal
#' @keywords internal
resolve_sc_params <- function(process_fit) {
  if (inherits(process_fit, "ldmppr_fit")) {
    p <- coef(process_fit)
    return(as.numeric(p))
  }
  if (is.numeric(process_fit) && length(process_fit) == 8 && !anyNA(process_fit) && all(process_fit[2:8] >= 0)) {
    return(as.numeric(process_fit))
  }
  stop("`process_fit` must be an ldmppr_fit object or a numeric vector of length 8.", call. = FALSE)
}


#' @rdname ldmppr-internal
#' @keywords internal
resolve_reference_ppp <- function(reference_data, process_fit, xy_bounds) {
  if (!is.null(reference_data)) {
    if (!spatstat.geom::is.ppp(reference_data)) {
      stop("If provided, `reference_data` must be a spatstat.geom::ppp object.", call. = FALSE)
    }
    return(reference_data)
  }

  if (!inherits(process_fit, "ldmppr_fit")) {
    stop("If `reference_data` is NULL, `process_fit` must be an ldmppr_fit containing `data_original` or `data` with (x,y,size).", call. = FALSE)
  }

  df <- NULL
  if (!is.null(process_fit$data_original)) df <- process_fit$data_original
  if (is.null(df) && !is.null(process_fit$data)) df <- process_fit$data

  if (is.null(df) || !is.data.frame(df)) {
    stop("Could not derive reference data: `process_fit$data_original` (preferred) or `process_fit$data` must be present as a data.frame.", call. = FALSE)
  }
  if (!all(c("x", "y") %in% names(df))) {
    stop("Derived reference data must contain columns `x` and `y`.", call. = FALSE)
  }
  if (!("size" %in% names(df))) {
    stop("Derived reference data must contain a `size` column to construct a marked reference pattern.", call. = FALSE)
  }
  if (is.null(xy_bounds) || length(xy_bounds) != 4) {
    stop("When deriving `reference_data` internally, please supply `xy_bounds = c(a_x,b_x,a_y,b_y)` (or provide a `reference_data` ppp so bounds can be inferred).", call. = FALSE)
  }

  generate_mpp(
    locations = df[, c("x", "y")],
    marks = df$size,
    xy_bounds = xy_bounds
  )
}


# ------------------------
# helpers for simulate_mpp
# ------------------------

#' @rdname ldmppr-internal
#' @keywords internal
.as_sc_params <- function(process_fit) {
  if (inherits(process_fit, "ldmppr_fit")) {
    p <- tryCatch(stats::coef(process_fit), error = function(e) NULL)
    if (is.null(p)) p <- process_fit$fit$solution
    return(as.numeric(p))
  }
  as.numeric(process_fit)
}


#' @rdname ldmppr-internal
#' @keywords internal
.infer_xy_bounds <- function(process_fit) {
  if (!inherits(process_fit, "ldmppr_fit")) return(NULL)
  ub <- process_fit$grid$upper_bounds %||% NULL
  if (is.null(ub) || length(ub) != 3) return(NULL)
  # upper_bounds = c(b_t, b_x, b_y) with implicit 0 lower bounds
  c(0, ub[2], 0, ub[3])
}


#' @rdname ldmppr-internal
#' @keywords internal
.infer_anchor_point <- function(process_fit) {
  if (!inherits(process_fit, "ldmppr_fit")) return(NULL)

  # Prefer original data if it exists and has x,y
  d0 <- process_fit$data_original %||% NULL
  if (is.data.frame(d0) && all(c("x", "y") %in% names(d0)) && nrow(d0) >= 1) {
    return(c(as.numeric(d0$x[1]), as.numeric(d0$y[1])))
  }

  # Fallback: transformed data matrix (time,x,y)
  d <- process_fit$data %||% NULL
  if (is.matrix(d) && ncol(d) >= 3 && nrow(d) >= 1) {
    return(c(as.numeric(d[1, 2]), as.numeric(d[1, 3])))
  }

  NULL
}


## ------------------------ ldmppr_budgets helpers ------------------------
#' @rdname ldmppr-internal
#' @keywords internal
new_ldmppr_budgets <- function(global_options,
                               local_budget_first_level,
                               local_budget_refinement_levels = NULL) {

  x <- list(
    global_options = global_options %||% list(),
    local_budget_first_level = local_budget_first_level %||% list(),
    local_budget_refinement_levels = local_budget_refinement_levels
  )
  class(x) <- "ldmppr_budgets"
  .validate_ldmppr_budgets(x)
  x
}

#' @rdname ldmppr-internal
#' @keywords internal
is_ldmppr_budgets <- function(x) inherits(x, "ldmppr_budgets")

#' @rdname ldmppr-internal
#' @keywords internal
as_ldmppr_budgets <- function(x, ...) {
  if (is_ldmppr_budgets(x)) return(x)

  # allow coercion from a plain list with the expected names
  if (is.list(x) &&
      all(c("global_options", "local_budget_first_level") %in% names(x))) {
    return(new_ldmppr_budgets(
      global_options = x$global_options %||% list(),
      local_budget_first_level = x$local_budget_first_level %||% list(),
      local_budget_refinement_levels = x$local_budget_refinement_levels
    ))
  }

  stop("Cannot coerce to ldmppr_budgets. Provide an ldmppr_budgets object or a list with ",
       "global_options and local_budget_first_level.", call. = FALSE)
}

#' @rdname ldmppr-internal
#' @keywords internal
.validate_ldmppr_budgets <- function(b) {
  if (!is_ldmppr_budgets(b)) stop("Invalid ldmppr_budgets object.", call. = FALSE)

  if (!is.list(b$global_options)) {
    stop("ldmppr_budgets: global_options must be a list.", call. = FALSE)
  }
  if (!is.list(b$local_budget_first_level)) {
    stop("ldmppr_budgets: local_budget_first_level must be a list.", call. = FALSE)
  }
  if (!is.null(b$local_budget_refinement_levels) && !is.list(b$local_budget_refinement_levels)) {
    stop("ldmppr_budgets: local_budget_refinement_levels must be NULL or a list.", call. = FALSE)
  }

  invisible(TRUE)
}


## ------------------------ ldmppr_grids helpers ------------------------
#' @rdname ldmppr-internal
#' @keywords internal
new_ldmppr_grids <- function(levels, upper_bounds, labels = NULL, include_endpoints = TRUE) {
  x <- list(
    levels = levels,
    upper_bounds = upper_bounds,
    labels = labels %||% rep("", length(levels)),
    include_endpoints = isTRUE(include_endpoints)
  )
  class(x) <- "ldmppr_grids"
  .validate_ldmppr_grids(x)
  x
}

#' @rdname ldmppr-internal
#' @keywords internal
is_ldmppr_grids <- function(x) inherits(x, "ldmppr_grids")

#' @rdname ldmppr-internal
#' @keywords internal
as_ldmppr_grids <- function(x, ...) {
  if (is_ldmppr_grids(x)) return(x)

  if (is.list(x) && !is.null(x$levels) && !is.null(x$upper_bounds)) {
    return(new_ldmppr_grids(
      levels = x$levels,
      upper_bounds = x$upper_bounds,
      labels = x$labels %||% NULL,
      include_endpoints = x$include_endpoints %||% TRUE
    ))
  }

  stop("Cannot coerce to ldmppr_grids.", call. = FALSE)
}

#' @rdname ldmppr-internal
#' @keywords internal
.validate_ldmppr_grids <- function(g) {
  if (!is.list(g) || !is_ldmppr_grids(g)) stop("Invalid ldmppr_grids object.", call. = FALSE)

  ub <- g$upper_bounds
  if (is.null(ub) || length(ub) != 3L || anyNA(ub) || any(!is.finite(ub))) {
    stop("ldmppr_grids: upper_bounds must be finite numeric length 3: c(b_t, b_x, b_y).", call. = FALSE)
  }
  if (ub[1] <= 0 || ub[2] <= 0 || ub[3] <= 0) {
    stop("ldmppr_grids: upper_bounds must be positive.", call. = FALSE)
  }

  lvls <- g$levels
  if (!is.list(lvls) || length(lvls) < 1L) {
    stop("ldmppr_grids: must contain at least one grid level.", call. = FALSE)
  }

  for (i in seq_along(lvls)) {
    L <- lvls[[i]]
    if (!is.list(L) || !all(c("x", "y", "t") %in% names(L))) {
      stop("ldmppr_grids: each level must be a list with x, y, t.", call. = FALSE)
    }
    for (nm in c("x", "y", "t")) {
      v <- L[[nm]]
      if (!is.numeric(v) || length(v) < 2L || anyNA(v) || any(!is.finite(v))) {
        stop("ldmppr_grids: each of x, y, t must be numeric vectors (len>=2), finite, non-NA.", call. = FALSE)
      }
      if (is.unsorted(v, strictly = FALSE)) {
        stop("ldmppr_grids: each of x, y, t must be non-decreasing (sorted).", call. = FALSE)
      }
    }

    if (min(L$x) < 0 || max(L$x) > ub[2]) stop("ldmppr_grids: x grid exceeds upper_bounds.", call. = FALSE)
    if (min(L$y) < 0 || max(L$y) > ub[3]) stop("ldmppr_grids: y grid exceeds upper_bounds.", call. = FALSE)
    if (min(L$t) < 0 || max(L$t) > ub[1]) stop("ldmppr_grids: t grid exceeds upper_bounds.", call. = FALSE)
  }

  invisible(TRUE)
}

#' @rdname ldmppr-internal
#' @keywords internal
.ldmppr_make_grid_schedule <- function(upper_bounds,
                                       levels,
                                       labels = NULL,
                                       include_endpoints = TRUE) {
  if (is.null(upper_bounds) || length(upper_bounds) != 3L ||
      anyNA(upper_bounds) || any(!is.finite(upper_bounds))) {
    stop("upper_bounds must be finite numeric length 3: c(b_t, b_x, b_y).", call. = FALSE)
  }
  if (!is.list(levels) || length(levels) < 1L) {
    stop("levels must be a non-empty list.", call. = FALSE)
  }

  mk_level <- function(level) {
    if (is.numeric(level) && length(level) == 3L) {
      nx <- as.integer(level[1]); ny <- as.integer(level[2]); nt <- as.integer(level[3])
      if (anyNA(c(nx, ny, nt)) || any(c(nx, ny, nt) < 2L)) {
        stop("Grid level counts must be >= 2.", call. = FALSE)
      }
      return(list(
        x = seq(0, upper_bounds[2], length.out = nx),
        y = seq(0, upper_bounds[3], length.out = ny),
        t = seq(0, upper_bounds[1], length.out = nt)
      ))
    }

    if (is.list(level) && all(c("nx", "ny", "nt") %in% names(level))) {
      nx <- as.integer(level$nx); ny <- as.integer(level$ny); nt <- as.integer(level$nt)
      if (anyNA(c(nx, ny, nt)) || any(c(nx, ny, nt) < 2L)) {
        stop("Grid level nx/ny/nt must be >= 2.", call. = FALSE)
      }
      return(list(
        x = seq(0, upper_bounds[2], length.out = nx),
        y = seq(0, upper_bounds[3], length.out = ny),
        t = seq(0, upper_bounds[1], length.out = nt)
      ))
    }

    if (is.list(level) && all(c("x", "y", "t") %in% names(level))) {
      return(list(
        x = as.numeric(level$x),
        y = as.numeric(level$y),
        t = as.numeric(level$t)
      ))
    }

    stop("Each grid level must be c(nx,ny,nt), list(nx,ny,nt), or list(x,y,t).", call. = FALSE)
  }

  lev_out <- lapply(levels, mk_level)

  if (!is.null(labels)) {
    if (!is.character(labels) || length(labels) != length(lev_out)) {
      stop("labels must be NULL or a character vector with length equal to levels.", call. = FALSE)
    }
  }

  list(
    levels = lev_out,
    upper_bounds = as.numeric(upper_bounds),
    labels = labels %||% rep("", length(lev_out)),
    include_endpoints = isTRUE(include_endpoints)
  )
}
