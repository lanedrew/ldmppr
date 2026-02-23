#' Mark model object
#'
#' \code{ldmppr_mark_model} objects store a fitted mark model and preprocessing
#' information used to predict marks at new locations and times.
#' These objects are typically returned by \code{\link{train_mark_model}} and can be
#' saved/loaded with \code{\link{save_mark_model}} and \code{\link{load_mark_model}}.
#'
#' @details
#' The model may be backed by different engines (currently \code{"xgboost"} and
#' \code{"ranger"}). For \code{"xgboost"}, the object can store a serialized booster payload
#' to make saving/loading robust across R sessions.
#'
#' @return
#' \describe{
#'    \item{\code{print()}}{prints a brief summary.}
#'    \item{\code{predict()}}{returns numeric predictions for new data.}
#' }
#'
#' @name ldmppr_mark_model
#' @rdname ldmppr_mark_model
#' @docType class
NULL


#' Create a mark model object
#'
#' @param engine character string (currently \code{"xgboost"} and \code{"ranger"}).
#' @param fit_engine fitted engine object (e.g. \code{xgb.Booster} or a ranger fit).
#' @param xgb_raw raw xgboost payload (e.g. UBJ) used for rehydration.
#' @param recipe a prepped recipes object used for preprocessing new data.
#' @param outcome outcome column name (default \code{"size"}).
#' @param feature_names (optional) vector of predictor names required at prediction time.
#' @param rasters (optional) list of rasters used for prediction (e.g. for spatial covariates).
#' @param info (optional) list of metadata.
#'
#' @describeIn ldmppr_mark_model Create a mark model container.
#' @export
ldmppr_mark_model <- function(engine,
                              fit_engine = NULL,
                              xgb_raw = NULL,
                              recipe = NULL,
                              outcome = "size",
                              feature_names = NULL,
                              rasters = NULL,
                              info = list()) {
  stopifnot(is.character(engine), length(engine) == 1)

  if (!engine %in% c("xgboost", "ranger")) {
    stop("engine must be 'xgboost' or 'ranger'.")
  }

  if (engine == "xgboost" && is.null(fit_engine) && is.null(xgb_raw)) {
    stop("For engine='xgboost', provide `fit_engine` (xgb.Booster) or `xgb_raw`.")
  }
  if (engine == "ranger" && is.null(fit_engine)) {
    stop("For engine='ranger', provide `fit_engine` (ranger fit).")
  }

  new_ldmppr_mark_model(engine, fit_engine, xgb_raw, recipe, outcome, feature_names, rasters, info)
}


#' @describeIn ldmppr_mark_model Print a brief summary of the mark model.
#' @param x a \code{ldmppr_mark_model} object.
#' @param ... additional arguments (not used).
#'
#' @export
print.ldmppr_mark_model <- function(x, ...) {
  .raster_label <- function(r) {
    nm <- names(r)
    if (!length(nm)) return("<unnamed>")
    if (length(nm) == 1L) return(nm)
    paste0(nm[1], " (+", length(nm) - 1L, " more)")
  }

  cat("ldmppr Mark Model\n")
  .cat_wrapped_field("  engine:           ", x$engine %||% NA_character_)
  cat("  has_fit_engine:   ", !is.null(x$fit_engine), "\n", sep = "")
  cat("  has_xgb_raw:      ", !is.null(x$xgb_raw), "\n", sep = "")
  if (!is.null(x$feature_names)) cat("  n_features:       ", length(x$feature_names), "\n", sep = "")
  if (!is.null(x$rasters)) {
    cat("  n_rasters:        ", length(x$rasters), "\n", sep = "")
    labels <- vapply(x$rasters, .raster_label, character(1))
    .cat_wrapped_field(
      "  raster_names:     ",
      paste0(paste(utils::head(labels, 4L), collapse = ", "), if (length(labels) > 4L) " ..." else "")
    )
  }
  if (!is.null(x$info$scaled_rasters)) cat("  scaled_rasters:   ", x$info$scaled_rasters, "\n", sep = "")
  if (!is.null(x$info$include_comp_inds)) cat("  comp_indices:     ", x$info$include_comp_inds, "\n", sep = "")
  invisible(x)
}

#' @describeIn ldmppr_mark_model Summarize a mark model.
#' @param object a \code{ldmppr_mark_model} object.
#' @param ... additional arguments (not used).
#' @export
summary.ldmppr_mark_model <- function(object, ...) {
  out <- list(
    engine = object$engine,
    n_features = length(object$feature_names %||% character(0)),
    n_rasters = length(object$rasters %||% list()),
    scaled_rasters = object$info$scaled_rasters %||% NA,
    include_comp_inds = object$info$include_comp_inds %||% NA,
    competition_radius = object$info$competition_radius %||% NA_real_,
    edge_correction = object$info$edge_correction %||% NA_character_,
    cv_folds = object$info$cv_folds %||% NA_integer_,
    tuning_grid_size = object$info$tuning_grid_size %||% NA_integer_,
    selection_metric = object$info$selection_metric %||% NA_character_
  )
  class(out) <- "summary.ldmppr_mark_model"
  out
}

#' @describeIn ldmppr_mark_model Print a summary produced by \code{\link{summary.ldmppr_mark_model}}.
#' @param x an object of class \code{summary.ldmppr_mark_model}.
#' @param ... additional arguments (not used).
#' @export
print.summary.ldmppr_mark_model <- function(x, ...) {
  cat("Summary: ldmppr Mark Model\n")
  .cat_wrapped_field("  engine:           ", x$engine %||% NA_character_)
  cat("  n_features:       ", x$n_features, "\n", sep = "")
  cat("  n_rasters:        ", x$n_rasters, "\n", sep = "")
  cat("  scaled_rasters:   ", x$scaled_rasters, "\n", sep = "")
  cat("  comp_indices:     ", x$include_comp_inds, "\n", sep = "")
  cat("  comp_radius:      ", x$competition_radius, "\n", sep = "")
  .cat_wrapped_field("  edge_correction:  ", x$edge_correction)
  cat("  cv_folds:         ", x$cv_folds, "\n", sep = "")
  cat("  tuning_grid_size: ", x$tuning_grid_size, "\n", sep = "")
  .cat_wrapped_field("  selection_metric: ", x$selection_metric)
  invisible(x)
}


#' @describeIn ldmppr_mark_model Predict marks for new data.
#' @param object a \code{ldmppr_mark_model} object.
#' @param new_data a data frame of predictors (and possibly outcome columns).
#'   Ignored when \code{sim_realization} is supplied.
#' @param sim_realization optional simulation realization containing \code{x}, \code{y}, and \code{time}.
#'   When supplied, predictors are built from rasters and optional competition indices.
#' @param raster_list optional list of rasters used when \code{sim_realization} is supplied.
#'   If omitted, uses rasters stored in \code{object} when available.
#' @param scaled_rasters \code{TRUE} or \code{FALSE}; whether supplied rasters are pre-scaled.
#' @param xy_bounds domain bounds \code{c(a_x, b_x, a_y, b_y)} used for competition indices.
#' @param include_comp_inds \code{TRUE} or \code{FALSE}; include competition-index features.
#' @param competition_radius positive numeric distance used when \code{include_comp_inds = TRUE}.
#' @param edge_correction edge correction for competition indices (\code{"none"} or \code{"toroidal"}).
#' @param seed optional nonnegative integer seed.
#' @param ... additional arguments.
#'
#' @importFrom stats predict
#' @export
predict.ldmppr_mark_model <- function(object,
                                      new_data = NULL,
                                      sim_realization = NULL,
                                      raster_list = NULL,
                                      scaled_rasters = FALSE,
                                      xy_bounds = NULL,
                                      include_comp_inds = FALSE,
                                      competition_radius = 15,
                                      edge_correction = "none",
                                      seed = NULL,
                                      ...) {
  if (!is.null(sim_realization)) {
    if (is.null(raster_list)) {
      raster_list <- object$rasters
      if (!is.null(object$info) && !is.null(object$info$scaled_rasters) && isTRUE(object$info$scaled_rasters)) {
        scaled_rasters <- TRUE
      }
    }

    new_data <- .build_mark_predictors(
      sim_realization = sim_realization,
      raster_list = raster_list,
      scaled_rasters = scaled_rasters,
      xy_bounds = xy_bounds,
      include_comp_inds = include_comp_inds,
      competition_radius = competition_radius,
      edge_correction = edge_correction
    )

    if (!is.null(seed)) {
      if (is.na(seed) || seed < 0 || seed != as.integer(seed)) {
        stop("Provide a nonnegative integer `seed`.", call. = FALSE)
      }
      set.seed(as.integer(seed))
    }
  }

  if (is.null(new_data)) {
    stop("Provide `new_data`, or provide `sim_realization` with raster/context arguments.", call. = FALSE)
  }

  baked <- preprocess_new_data(object, new_data)

  if (object$engine == "xgboost") {
    booster <- rehydrate_xgb(object)
    xmat <- as.matrix(baked)
    return(stats::predict(booster, xmat))
  }

  if (object$engine == "ranger") {
    fit <- object$fit_engine
    return(stats::predict(fit, data = baked, ...)$predictions)
  }

  stop("Unknown engine: ", object$engine)
}

#' Save a mark model to disk
#'
#' @param object a \code{ldmppr_mark_model} object.
#' @param path file path to write an \code{.rds}.
#' @param ... passed to methods.
#'
#' @describeIn ldmppr_mark_model Save a mark model to disk.
#' @export
save_mark_model <- function(object, path, ...) UseMethod("save_mark_model")


#' @describeIn ldmppr_mark_model Save method for \code{ldmppr_mark_model}.
#' @export
save_mark_model.ldmppr_mark_model <- function(object, path, ...) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  if (object$engine == "xgboost") {
    # Ensure raw bytes exist; store UBJ for robustness
    if (is.null(object$xgb_raw)) {
      booster <- object$fit_engine %||% rehydrate_xgb(object)
      object$xgb_raw <- xgboost::xgb.save.raw(booster, raw_format = "ubj")
    }
    # Drop live booster from serialization (optional but recommended)
    object$fit_engine <- NULL
  }

  saveRDS(object, file = path)
  invisible(path)
}


#' Load a saved mark model
#'
#' @param path path to an \code{.rds} created by \code{\link{save_mark_model}} (or legacy objects).
#' @return an object of class \code{"ldmppr_mark_model"}.
#'
#' @describeIn ldmppr_mark_model Load a saved mark model from disk.
#' @export
load_mark_model <- function(path) {
  if (!is.character(path) || length(path) != 1 || !file.exists(path)) {
    stop("`path` must be an existing .rds file.", call. = FALSE)
  }

  obj <- readRDS(path)

  # If already in new format, return
  if (inherits(obj, "ldmppr_mark_model")) return(obj)

  # Legacy formats: bundle/parsnip/workflow/etc.
  if (exists("as_mark_model", mode = "function")) {
    return(as_mark_model(obj))
  }

  stop("File does not contain an ldmppr_mark_model object.", call. = FALSE)
}
