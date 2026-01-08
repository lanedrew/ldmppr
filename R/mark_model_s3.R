#' Mark model object
#'
#' `ldmppr_mark_model` objects store a fitted mark model and preprocessing
#' information used to predict marks at new locations and times.
#'
#' These objects are typically returned by [train_mark_model()] and can be
#' saved/loaded with [save_mark_model()] and [load_mark_model()].
#'
#' @details
#' The model may be backed by different engines (currently `"xgboost"` and
#' `"ranger"`). For xgboost, the object can store a serialized booster payload
#' to make saving/loading robust across R sessions.
#'
#' @return
#' * `print()` prints a brief summary.
#' * `predict()` returns numeric predictions for new data.
#'
#' @name ldmppr_mark_model
#' @rdname ldmppr_mark_model
#' @docType class
NULL


#' Create a mark model object
#'
#' @param engine Character scalar. One of `"xgboost"` or `"ranger"`.
#' @param fit_engine Fitted engine object (e.g. `xgb.Booster` or a ranger fit).
#' @param xgb_raw Raw xgboost payload (e.g. UBJ) used for rehydration.
#' @param recipe A prepped recipes object used for preprocessing new data.
#' @param outcome Outcome column name (default `"size"`).
#' @param feature_names Optional vector of predictor names required at prediction time.
#' @param info Optional list of metadata.
#'
#' @describeIn ldmppr_mark_model Create a mark model container.
#' @export
ldmppr_mark_model <- function(engine,
                              fit_engine = NULL,
                              xgb_raw = NULL,
                              recipe = NULL,
                              outcome = "size",
                              feature_names = NULL,
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

  new_ldmppr_mark_model(engine, fit_engine, xgb_raw, recipe, outcome, feature_names, info)
}


#' @describeIn ldmppr_mark_model Print a brief summary of the mark model.
#' @param x a `ldmppr_mark_model` object.
#' @param ... additional arguments (not used).
#'
#' @export
print.ldmppr_mark_model <- function(x, ...) {
  cat("<ldmppr_mark_model>\n")
  cat("  engine: ", x$engine, "\n", sep = "")
  cat("  has fit_engine: ", !is.null(x$fit_engine), "\n", sep = "")
  cat("  has xgb_raw: ", !is.null(x$xgb_raw), "\n", sep = "")
  if (!is.null(x$feature_names)) cat("  n_features: ", length(x$feature_names), "\n", sep = "")
  invisible(x)
}


#' @describeIn ldmppr_mark_model Predict marks for new data.
#' @param object a `ldmppr_mark_model` object.
#' @param new_data a data frame of predictors (and possibly outcome columns).
#' @param ... additional arguments.
#'
#' @export
predict.ldmppr_mark_model <- function(object, new_data, ...) {
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
#' @param object a `ldmppr_mark_model` object.
#' @param path file path to write an `.rds`.
#' @param ... passed to methods.
#'
#' @describeIn ldmppr_mark_model Save a mark model to disk.
#' @export
save_mark_model <- function(object, path, ...) UseMethod("save_mark_model")


#' @describeIn ldmppr_mark_model Save method for `ldmppr_mark_model`.
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
#' @param path path to an `.rds` created by [save_mark_model()] (or legacy objects).
#' @return an object of class `"ldmppr_mark_model"`.
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
