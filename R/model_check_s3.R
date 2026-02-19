#' Model fit diagnostic object
#'
#' Objects of class \code{ldmppr_model_check} are returned by \code{\link{check_model_fit}}.
#' They contain global envelope test results and curve sets for multiple summary
#' functions/statistics.
#'
#' @details
#' An \code{ldmppr_model_check} is a list with components such as:
#' \itemize{
#'   \item \code{combined_env}: a global envelope test object (typically from **GET**)
#'   \item \code{envs}: named list of envelope test objects (e.g., \code{L}, \code{F}, \code{G}, \code{J}, \code{E}, \code{V})
#'   \item \code{curve_sets}: named list of curve set objects
#'   \item \code{settings}: list of settings used when generating envelopes (e.g., \code{n_sim}, \code{thinning})
#' }
#'
#' @return
#' \describe{
#'    \item{\code{print()}}{prints a brief summary of the diagnostic object.}
#'    \item{\code{summary()}}{returns a \code{summary.ldmppr_model_check} object.}
#'    \item{\code{plot()}}{plots the combined envelope or a selected statistic envelope.}
#' }
#'
#' @name ldmppr_model_check
#' @rdname ldmppr_model_check
#' @docType class
NULL


#' @describeIn ldmppr_model_check Print a brief summary of the diagnostic results.
#' @param x an object of class \code{ldmppr_model_check}.
#' @param ... additional arguments (not used).
#'
#' @export
print.ldmppr_model_check <- function(x, ...) {
  p_comb <- attributes(x$combined_env)$p %||% x$combined_env$p_value %||% NA_real_
  cat("ldmppr Model Check\n")
  if (!is.null(x$settings$n_sim)) cat("  n_sim:            ", x$settings$n_sim, "\n", sep = "")
  if (!is.null(x$settings$n_sim_used)) cat("  n_sim_used:       ", x$settings$n_sim_used, "\n", sep = "")
  if (!is.null(x$settings$thinning)) cat("  thinning:         ", x$settings$thinning, "\n", sep = "")
  if (!is.null(x$settings$edge_correction)) cat("  edge_correction:  ", x$settings$edge_correction, "\n", sep = "")
  if (!is.na(p_comb)) cat("  p_combined:       ", signif(p_comb, 6), "\n", sep = "")
  cat("  envelopes:        ", paste(names(x$envs), collapse = ", "), "\n", sep = "")
  invisible(x)
}


#' @describeIn ldmppr_model_check Summarize p-values from the combined and individual envelopes.
#' @param object an object of class \code{ldmppr_model_check}.
#' @param ... additional arguments (not used).
#'
#' @export
summary.ldmppr_model_check <- function(object, ...) {
  # GET global_envelope_test objects usually have p-values; be conservative if field names differ
  p_comb <- attributes(object$combined_env)$p %||% object$combined_env$p_value %||% NA_real_
  per_stat <- vapply(object$envs, function(e) attributes(e)$p %||% e$p_value %||% NA_real_, numeric(1))
  alpha <- attributes(object$combined_env)$alpha %||% 0.05

  out <- list(
    p_combined = p_comb,
    alpha = alpha,
    p_individual = per_stat,
    settings = object$settings
  )
  class(out) <- "summary.ldmppr_model_check"
  out
}


#' @describeIn ldmppr_model_check Print a summary produced by \code{\link{summary.ldmppr_model_check}}.
#' @param x an object of class \code{summary.ldmppr_model_check}.
#' @param ... additional arguments (not used).
#'
#' @export
print.summary.ldmppr_model_check <- function(x, ...) {
  cat("Summary: ldmppr Model Check\n")
  cat("  p_combined:       ", signif(x$p_combined, 6), "\n", sep = "")
  cat("  alpha:            ", x$alpha, "\n", sep = "")
  cat("  p_individual:\n")
  print(signif(x$p_individual, 6))
  invisible(x)
}


#' @describeIn ldmppr_model_check Plot the combined envelope or a selected statistic.
#' @param x an object of class \code{ldmppr_model_check}.
#' @param which which envelope to plot. \code{"combined"} plots the global envelope;
#' otherwise one of \code{"L"}, \code{"F"}, \code{"G"}, \code{"J"}, \code{"E"}, \code{"V"}.
#' @param ... additional arguments passed to the underlying \code{plot()} method (e.g., from **GET**).
#'
#' @export
plot.ldmppr_model_check <- function(x,
                                    which = c("combined", "L", "F", "G", "J", "E", "V"),
                                    ...) {
  which <- match.arg(which)
  dots <- list(...)

  .render_env_plot <- function(env_obj, default_main) {
    out <- do.call(plot, c(list(env_obj), dots))

    if (inherits(out, "ggplot")) {
      if (is.null(dots$main)) out <- out + ggplot2::labs(title = default_main)
      print(out)
      return(invisible(out))
    }

    if (is.null(dots$main)) {
      try(graphics::title(main = default_main), silent = TRUE)
    }
    invisible(out)
  }

  if (which == "combined") {
    return(.render_env_plot(x$combined_env, paste0("ldmppr Model Check: Combined Envelope (p = ", signif(attributes(x$combined_env)$p %||% x$combined_env$p_value %||% NA_real_, 6), ")")))
  }
  if (is.null(x$envs[[which]])) {
    stop("Envelope '", which, "' is not available. Available: ", paste(names(x$envs), collapse = ", "), call. = FALSE)
  }
  .render_env_plot(x$envs[[which]], paste0("ldmppr Model Check: ", which, " Envelope (p = ", signif(attributes(x$envs[[which]])$p %||% x$envs[[which]]$p_value %||% NA_real_, 6), ")"))
}
