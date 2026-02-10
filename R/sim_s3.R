#' Simulated marked point process object
#'
#' \code{ldmppr_sim} objects are returned by \code{\link{simulate_mpp}}. They contain the simulated
#' realization, an associated marked point pattern object, and metadata used to
#' reproduce or inspect the simulation.
#'
#' @details
#' An \code{ldmppr_sim} is a list with at least:
#' \itemize{
#'   \item \code{process}: process name (e.g. \code{"self_correcting"})
#'   \item \code{mpp}: a marked point pattern object
#'   \item \code{realization}: data.frame with columns \code{time}, \code{x}, \code{y}, \code{marks}
#'   \item \code{params}, \code{bounds}, and other metadata
#' }
#'
#' @return For methods:
#' \describe{
#'   \item{\code{print()}}{prints a summary of the simulation.}
#'   \item{\code{plot()}}{returns a ggplot visualization of the marked point pattern.}
#'   \item{\code{as.data.frame()}}{returns the simulated realization as a data.frame.}
#'   \item{\code{nobs()}}{returns the number of points in the realization.}
#'   \item{\code{mpp()}}{returns the marked point pattern object.}
#' }
#'
#' @name ldmppr_sim
#' @rdname ldmppr_sim
#' @docType class
NULL


#' @describeIn ldmppr_sim Print a brief summary of the simulation.
#' @param x a \code{ldmppr_sim} object.
#' @param ... additional arguments (not used).
#' @export
print.ldmppr_sim <- function(x, ...) {
  cat("ldmppr sim: \n")
  cat("  process: ", x$process, "\n", sep = "")
  cat("  n:       ", nrow(x$realization), "\n", sep = "")
  cat("  thinning:", x$thinning, "\n", sep = "")
  cat("  correction: ", x$correction, "\n", sep = "")
  cat("  time bounds: [", x$bounds$t_min, ", ", x$bounds$t_max, "]\n", sep = "")
  cat("  xy bounds:   [", paste(x$bounds$xy_bounds, collapse = ", "), "]\n", sep = "")
  invisible(x)
}

#' @describeIn ldmppr_sim Coerce to a data.frame of the simulated realization.
#' @param x a \code{ldmppr_sim} object.
#' @param ... additional arguments (not used).
#' @export
as.data.frame.ldmppr_sim <- function(x, ...) {
  x$realization
}

#' @describeIn ldmppr_sim Number of simulated points.
#' @importFrom stats nobs
#' @param object a \code{ldmppr_sim} object.
#' @param ... additional arguments (not used).
#' @export
nobs.ldmppr_sim <- function(object, ...) {
  nrow(object$realization)
}

#' @describeIn ldmppr_sim Plot the simulated marked point pattern.
#' @param x a \code{ldmppr_sim} object.
#' @param pattern_type type of pattern to plot \code{"simulated"} (default).
#' @param ... additional arguments passed to \code{plot_mpp}.
#' @export
plot.ldmppr_sim <- function(x, pattern_type = "simulated", ...) {
  plot_mpp(x$mpp, pattern_type = pattern_type, ...)
}

#' @describeIn ldmppr_sim Extract the underlying marked point pattern object.
#' @param x a \code{ldmppr_sim} object.
#' @param ... additional arguments (not used).
#' @export
mpp.ldmppr_sim <- function(x, ...) x$mpp
