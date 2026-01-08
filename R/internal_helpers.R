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
                           correction,
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
      correction = correction,
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
                           engine = "nloptr",
                           call = NULL,
                           timing = NULL) {
  structure(
    list(
      process = process,
      engine = engine,
      fit = fit,          # best fit (nloptr object)
      fits = fits,        # list of nloptr fits (e.g. per-delta) or NULL
      mapping = mapping,  # list(delta, delta_values, chosen_index, objectives)
      grid = grid,        # list(x_grid, y_grid, t_grid, upper_bounds)
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
