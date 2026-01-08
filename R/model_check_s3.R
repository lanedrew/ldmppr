#' Model fit diagnostic object
#'
#' Objects of class `ldmppr_model_check` are returned by [check_model_fit()].
#' They contain global envelope test results and curve sets for multiple summary
#' functions/statistics.
#'
#' @details
#' An `ldmppr_model_check` is a list with components such as:
#' \itemize{
#'   \item `combined_env`: a global envelope test object (typically from **GET**)
#'   \item `envs`: named list of envelope test objects (e.g., `L`, `F`, `G`, `J`, `E`, `V`)
#'   \item `curve_sets`: named list of curve set objects
#'   \item `settings`: list of settings used when generating envelopes (e.g., `n_sim`, `thinning`)
#' }
#'
#' @return
#' * `print()` prints a brief summary of the diagnostic object.
#' * `summary()` returns a `summary.ldmppr_model_check` object.
#' * `plot()` plots the combined envelope or a selected statistic envelope.
#'
#' @name ldmppr_model_check
#' @rdname ldmppr_model_check
#' @docType class
NULL


#' @describeIn ldmppr_model_check Print a brief summary of the diagnostic results.
#' @param x an object of class `ldmppr_model_check`.
#' @param ... additional arguments (not used).
#'
#' @export
print.ldmppr_model_check <- function(x, ...) {
  cat("<ldmppr_model_check>\n")
  if (!is.null(x$settings$n_sim)) cat("  n_sim: ", x$settings$n_sim, "\n", sep = "")
  if (!is.null(x$settings$thinning)) cat("  thinning: ", x$settings$thinning, "\n", sep = "")
  if (!is.null(x$settings$correction)) cat("  correction: ", x$settings$correction, "\n", sep = "")
  cat("  envelopes: ", paste(names(x$envs), collapse = ", "), "\n", sep = "")
  invisible(x)
}


#' @describeIn ldmppr_model_check Summarize p-values from the combined and individual envelopes.
#' @param object an object of class `ldmppr_model_check`.
#' @param ... additional arguments (not used).
#'
#' @export
summary.ldmppr_model_check <- function(object, ...) {
  # GET global_envelope_test objects usually have p-values; be conservative if field names differ
  p_comb <- object$combined_env$p %||% object$combined_env$p_value %||% NA_real_
  per_stat <- vapply(object$envs, function(e) e$p %||% e$p_value %||% NA_real_, numeric(1))

  out <- list(
    p_combined = p_comb,
    p_individual = per_stat,
    settings = object$settings
  )
  class(out) <- "summary.ldmppr_model_check"
  out
}


#' @describeIn ldmppr_model_check Print a summary produced by [summary.ldmppr_model_check()].
#' @param x an object of class `summary.ldmppr_model_check`.
#' @param ... additional arguments (not used).
#'
#' @export
print.summary.ldmppr_model_check <- function(x, ...) {
  cat("<summary: ldmppr_model_check>\n")
  cat("  combined p-value: ", x$p_combined, "\n", sep = "")
  cat("  individual p-values:\n")
  print(x$p_individual)
  invisible(x)
}


#' @describeIn ldmppr_model_check Plot the combined envelope or a selected statistic.
#' @param x an object of class `ldmppr_model_check`.
#' @param which which envelope to plot. `"combined"` plots the global envelope; otherwise one of `"L"`, `"F"`, `"G"`, `"J"`, `"E"`, `"V"`.
#' @param ... additional arguments passed to the underlying `plot()` method (e.g., from **GET**).
#'
#' @export
plot.ldmppr_model_check <- function(x,
                                    which = c("combined", "L", "F", "G", "J", "E", "V"),
                                    ...) {
  which <- match.arg(which)
  if (which == "combined") return(plot(x$combined_env, ...))
  plot(x$envs[[which]], ...)
}
