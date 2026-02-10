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
  cat("ldmppr fit: \n")
  cat("  process: ", x$process, "\n", sep = "")
  cat("  engine:  ", x$engine, "\n", sep = "")
  if (!is.null(x$data_summary$n)) cat("  n_obs:   ", x$data_summary$n, "\n", sep = "")
  if (!is.null(x$mapping$delta_values)) cat("  n_deltas:", length(x$mapping$delta_values), "\n", sep = "")
  if (!is.null(x$mapping$delta) && !is.na(x$mapping$delta)) cat("  delta*:  ", x$mapping$delta, "\n", sep = "")
  if (!is.null(x$fit$objective)) cat("  objective (best): ", signif(x$fit$objective, 8), "\n", sep = "")
  if (!is.null(x$fit$solution)) cat(" optimal parameters: ", signif(x$fit$solution, 8), "\n", sep = " ")
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
  cat("summary: ldmppr fit \n")
  cat("  process:  ", x$process, "\n", sep = "")
  cat("  engine:   ", x$engine, "\n", sep = "")
  cat("  status:   ", x$status, "\n", sep = "")
  cat("  objective:", signif(x$objective, 10), "\n", sep = "")
  if (!is.null(x$mapping$delta) && !is.na(x$mapping$delta)) cat("  delta*:   ", x$mapping$delta, "\n", sep = "")
  if (!is.null(x$timing$seconds)) cat("  seconds:  ", round(x$timing$seconds, 3), "\n", sep = "")
  cat("  solution:\n")
  print(x$solution)
  invisible(x)
}


#' @describeIn ldmppr_fit Plot diagnostics for a fitted model.
#' @param x an object of class \code{ldmppr_fit}.
#' @param ... additional arguments passed to \code{plot()}.
#'
#' @export
plot.ldmppr_fit <- function(x, ...) {
  if (is.null(x$mapping$delta_values) || is.null(x$mapping$objectives)) {
    graphics::plot.new()
    graphics::title(main = "ldmppr_fit (no multi-fit diagnostics stored)")
    return(invisible(x))
  }
  dv <- x$mapping$delta_values
  obj <- x$mapping$objectives
  graphics::plot(dv, obj, xlab = "delta", ylab = "objective (minimized)", ...)
  if (!is.null(x$mapping$delta) && !is.na(x$mapping$delta)) {
    graphics::abline(v = x$mapping$delta, lty = 2)
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
