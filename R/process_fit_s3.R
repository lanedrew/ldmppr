#' Fitted point-process model object
#'
#' Objects of class \code{ldmppr_fit} are returned by \code{\link{estimate_process_parameters}}.
#' They contain the best-fitting optimization result (and optionally multiple fits,
#' e.g. from a delta search) along with metadata used to reproduce the fit.
#'
#' @details
#' A \code{ldmppr_fit} is a list with (at minimum):
#' \itemize{
#'   \item \code{process}: process name (e.g. \code{"self_correcting"})
#'   \item \code{fit}: best optimization result (currently an \code{nloptr} object)
#'   \item \code{mapping}: mapping information (e.g. chosen \code{delta}, objectives)
#'   \item \code{grid}: grid definitions used by likelihood approximation
#' }
#'
#' @return
#' \describe{
#'   \item{\code{print()}}{prints a brief summary of the fit.}
#'   \item{\code{coef()}}{returns the estimated parameter vector.}
#'   \item{\code{logLik()}}{returns the log-likelihood at the optimum.}
#'   \item{\code{summary()}}{returns a \code{summary.ldmppr_fit}.}
#'   \item{\code{plot()}}{plots diagnostics for multi-fit runs, if available.}
#' }
#'
#' @name ldmppr_fit
#' @docType class
NULL



#' @describeIn ldmppr_fit Print a brief summary of a fitted model.
#' @param x an object of class \code{ldmppr_fit}.
#' @param ... additional arguments (not used).
#'
#' @export
print.ldmppr_fit <- function(x, ...) {
  cat("ldmppr Fit\n")
  cat("  process:         ", x$process %||% NA_character_, "\n", sep = "")
  cat("  engine:          ", x$engine %||% NA_character_, "\n", sep = "")
  if (!is.null(x$settings$strategy)) cat("  strategy:        ", x$settings$strategy, "\n", sep = "")
  if (!is.null(x$settings$global_algorithm)) cat("  global_alg:      ", x$settings$global_algorithm, "\n", sep = "")
  if (!is.null(x$settings$local_algorithm)) cat("  local_alg:       ", x$settings$local_algorithm, "\n", sep = "")
  if (!is.null(x$settings$starts) && is.list(x$settings$starts)) {
    st <- x$settings$starts
    cat("  starts:          ",
        "global=", st$global %||% NA_integer_,
        ", local=", st$local %||% NA_integer_,
        ", jitter_sd=", signif(st$jitter_sd %||% NA_real_, 4),
        ", seed=", st$seed %||% NA_integer_, "\n", sep = "")
  }
  if (!is.null(x$data_summary$n)) cat("  n_obs:           ", x$data_summary$n, "\n", sep = "")
  if (!is.null(x$mapping$delta_values)) cat("  n_deltas:        ", length(x$mapping$delta_values), "\n", sep = "")
  if (!is.null(x$mapping$delta) && !is.na(x$mapping$delta)) cat("  selected_delta:  ", signif(x$mapping$delta, 6), "\n", sep = "")
  if (!is.null(x$fit$objective)) cat("  objective:       ", signif(x$fit$objective, 8), "\n", sep = "")
  if (!is.null(x$fit$status)) cat("  final_status:    ", x$fit$status, "\n", sep = "")
  if (!is.null(x$fit$message) && is.character(x$fit$message) && nzchar(x$fit$message)) cat("  final_outcome:   ", x$fit$message, "\n", sep = "")
  if (!is.null(x$timing$seconds)) cat("  elapsed_sec:     ", round(x$timing$seconds, 3), "\n", sep = "")
  if (!is.null(x$fit$solution)) {
    cat("  coefficients:    ", paste(signif(as.numeric(x$fit$solution), 6), collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}

#' @importFrom stats coef
#' @describeIn ldmppr_fit Extract the estimated parameter vector.
#' @param object an object of class \code{ldmppr_fit}.
#' @param ... additional arguments (not used).
#'
#' @export
coef.ldmppr_fit <- function(object, ...) {
  object$fit$solution
}


#' @describeIn ldmppr_fit Log-likelihood at the optimum.
#' @param object an object of class \code{ldmppr_fit}.
#' @param ... additional arguments (not used).
#'
#' @export
logLik.ldmppr_fit <- function(object, ...) {
  # nloptr minimizes objective; your objective is negative log-likelihood
  ll <- -object$fit$objective
  class(ll) <- "logLik"
  attr(ll, "df") <- length(object$fit$solution)
  attr(ll, "nobs") <- object$data_summary$n %||% NA_integer_
  ll
}


#' @describeIn ldmppr_fit Summarize a fitted model.
#' @param object an object of class \code{ldmppr_fit}.
#' @param ... additional arguments (not used).
#'
#' @export
summary.ldmppr_fit <- function(object, ...) {
  out <- list(
    process = object$process,
    engine = object$engine,
    settings = object$settings,
    delta = object$mapping$delta,
    solution = object$fit$solution,
    objective = object$fit$objective,
    status = object$fit$status,
    message = object$fit$message %||% NULL,
    mapping = object$mapping,
    timing = object$timing
  )
  class(out) <- "summary.ldmppr_fit"
  out
}


#' @describeIn ldmppr_fit Print a summary produced by \code{\link{summary.ldmppr_fit}}.
#' @param x an object of class \code{summary.ldmppr_fit}.
#' @param ... additional arguments (not used).
#' @export
print.summary.ldmppr_fit <- function(x, ...) {
  cat("Summary: ldmppr Fit\n")
  cat("  process:         ", x$process %||% NA_character_, "\n", sep = "")
  cat("  engine:          ", x$engine %||% NA_character_, "\n", sep = "")
  if (!is.null(x$settings$strategy)) cat("  strategy:        ", x$settings$strategy, "\n", sep = "")
  if (!is.null(x$settings$starts) && is.list(x$settings$starts)) {
    st <- x$settings$starts
    cat("  starts:          ",
        "global=", st$global %||% NA_integer_,
        ", local=", st$local %||% NA_integer_,
        ", jitter_sd=", signif(st$jitter_sd %||% NA_real_, 4),
        ", seed=", st$seed %||% NA_integer_, "\n", sep = "")
  }
  cat("  status:          ", x$status %||% NA_character_, "\n", sep = "")
  if (!is.null(x$message) && is.character(x$message) && nzchar(x$message)) cat("  outcome:         ", x$message, "\n", sep = "")
  cat("  objective:       ", signif(x$objective, 10), "\n", sep = "")
  if (!is.null(x$mapping$delta) && !is.na(x$mapping$delta)) cat("  selected_delta:  ", signif(x$mapping$delta, 6), "\n", sep = "")
  if (!is.null(x$timing$seconds)) cat("  elapsed_sec:     ", round(x$timing$seconds, 3), "\n", sep = "")
  cat("  coefficients:\n")
  print(signif(x$solution, 8))
  invisible(x)
}


#' @describeIn ldmppr_fit Plot diagnostics for a fitted model.
#' @param x an object of class \code{ldmppr_fit}.
#' @param ... additional arguments passed to \code{plot()}.
#'
#' @export
plot.ldmppr_fit <- function(x, ...) {
  dots <- list(...)

  if (is.null(x$mapping$delta_values) || is.null(x$mapping$objectives)) {
    graphics::plot.new()
    graphics::title(main = "<ldmppr_fit>", sub = "No delta-profile diagnostics stored")
    return(invisible(x))
  }

  dv <- x$mapping$delta_values
  obj <- x$mapping$objectives

  if (is.null(dots$type)) dots$type <- "b"
  if (is.null(dots$pch)) dots$pch <- 16
  if (is.null(dots$lwd)) dots$lwd <- 1.5
  if (is.null(dots$col)) dots$col <- "#2C7FB8"
  if (is.null(dots$xlab)) dots$xlab <- "delta"
  if (is.null(dots$ylab)) dots$ylab <- "objective (smaller is better)"
  if (is.null(dots$main)) dots$main <- "<ldmppr_fit> delta profile"

  do.call(graphics::plot, c(list(x = dv, y = obj), dots))

  if (!is.null(x$mapping$delta) && !is.na(x$mapping$delta)) {
    graphics::abline(v = x$mapping$delta, lty = 2, col = "#D95F02")
    graphics::legend("topright",
                     legend = c("objective", "selected delta"),
                     col = c(dots$col, "#D95F02"),
                     lty = c(1, 2),
                     pch = c(16, NA),
                     bty = "n")
  }
  invisible(x)
}


#' Convert an ldmppr_fit to the underlying nloptr result
#'
#' @param x an object of class \code{ldmppr_fit}.
#' @param ... unused.
#'
#' @describeIn ldmppr_fit Extract the underlying \code{nloptr} result.
#' @export
as_nloptr <- function(x, ...) UseMethod("as_nloptr")

#' @describeIn ldmppr_fit Extract the underlying \code{nloptr} result.
#' @param x an object of class \code{ldmppr_fit}.
#' @param ... additional arguments (not used).
#'
#' @export
as_nloptr.ldmppr_fit <- function(x, ...) x$fit

#' @describeIn ldmppr_fit Number of observations used in the fitted model.
#' @importFrom stats nobs
#' @param object an object of class \code{ldmppr_fit}.
#' @param ... additional arguments (not used).
#' @export
nobs.ldmppr_fit <- function(object, ...) {
  object$data_summary$n %||% NA_integer_
}

#' @describeIn ldmppr_fit Coerce fit summary to a one-row data frame.
#' @param x an object of class \code{ldmppr_fit}.
#' @param ... additional arguments (not used).
#' @export
as.data.frame.ldmppr_fit <- function(x, ...) {
  out <- data.frame(
    process = x$process %||% NA_character_,
    engine = x$engine %||% NA_character_,
    strategy = x$settings$strategy %||% NA_character_,
    global_algorithm = x$settings$global_algorithm %||% NA_character_,
    local_algorithm = x$settings$local_algorithm %||% NA_character_,
    n_obs = x$data_summary$n %||% NA_integer_,
    selected_delta = x$mapping$delta %||% NA_real_,
    objective = x$fit$objective %||% NA_real_,
    status = x$fit$status %||% NA_integer_,
    elapsed_sec = x$timing$seconds %||% NA_real_,
    stringsAsFactors = FALSE
  )

  sol <- as.numeric(x$fit$solution %||% numeric(0))
  if (length(sol)) {
    nm <- names(x$fit$solution)
    if (is.null(nm) || any(!nzchar(nm))) nm <- paste0("param_", seq_along(sol))
    for (i in seq_along(sol)) out[[paste0("coef_", nm[i])]] <- sol[i]
  }
  out
}
