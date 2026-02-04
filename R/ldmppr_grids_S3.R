#' Grid schedule object
#'
#' Objects of class \code{ldmppr_grids} define one or more grid "levels" used by
#' \code{\link{estimate_process_parameters}}. Each level contains numeric vectors
#' \code{x}, \code{y}, and \code{t} defining the approximation grid. Levels are
#' typically ordered from coarse to fine.
#'
#' @details
#' A \code{ldmppr_grids} is a list with (at minimum):
#' \itemize{
#'   \item \code{levels}: list of levels; each level is a list with \code{x}, \code{y}, \code{t}
#'   \item \code{upper_bounds}: numeric \code{c(b_t, b_x, b_y)}
#'   \item \code{labels}: optional labels used only for printing
#'   \item \code{include_endpoints}: logical
#' }
#'
#' @return
#' \describe{
#'   \item{\code{print()}}{prints a brief description of bounds and grid levels.}
#'   \item{\code{summary()}}{returns a \code{summary.ldmppr_grids}.}
#'   \item{\code{as.data.frame()}}{returns one row per level with dimensions and ranges.}
#'   \item{\code{length()}}{returns the number of levels.}
#'   \item{\code{[ ]}}{subsets levels, preserving class.}
#'   \item{\code{as.list()}}{returns the underlying list structure.}
#' }
#'
#' @name ldmppr_grids-class
#' @docType class
NULL

#' @describeIn ldmppr_grids-class Print a brief summary of a grid schedule.
#' @param x an object of class \code{ldmppr_grids}.
#' @param ... additional arguments (not used).
#' @export
print.ldmppr_grids <- function(x, ...) {
  stopifnot(is_ldmppr_grids(x))
  ub <- x$upper_bounds
  cat("<ldmppr_grids>\n")
  cat("  upper_bounds: b_t=", ub[1], ", b_x=", ub[2], ", b_y=", ub[3], "\n", sep = "")
  cat("  levels: ", length(x$levels), "\n", sep = "")
  for (i in seq_along(x$levels)) {
    L <- x$levels[[i]]
    lab <- x$labels[[i]]
    if (nzchar(lab)) lab <- paste0(" (", lab, ")")
    cat("    - level ", i, lab, ": ",
        length(L$x), "x", length(L$y), "x", length(L$t),
        "  [x:", signif(min(L$x), 4), "..", signif(max(L$x), 4),
        ", y:", signif(min(L$y), 4), "..", signif(max(L$y), 4),
        ", t:", signif(min(L$t), 4), "..", signif(max(L$t), 4), "]\n",
        sep = "")
  }
  invisible(x)
}

#' @describeIn ldmppr_grids-class Summarize a grid schedule.
#' @param object an object of class \code{ldmppr_grids}.
#' @param ... additional arguments (not used).
#' @export
summary.ldmppr_grids <- function(object, ...) {
  stopifnot(is_ldmppr_grids(object))
  ub <- object$upper_bounds
  dims <- vapply(object$levels, function(L) c(nx = length(L$x), ny = length(L$y), nt = length(L$t)), numeric(3))
  dims <- t(dims)

  out <- list(
    upper_bounds = ub,
    n_levels = length(object$levels),
    dims = dims,
    labels = object$labels
  )
  class(out) <- "summary.ldmppr_grids"
  out
}

#' @describeIn ldmppr_grids-class Print a summary produced by \code{summary.ldmppr_grids()}.
#' @param x an object of class \code{summary.ldmppr_grids}.
#' @param ... additional arguments (not used).
#' @export
print.summary.ldmppr_grids <- function(x, ...) {
  cat("<summary: ldmppr_grids>\n")
  cat("  levels: ", x$n_levels, "\n", sep = "")
  cat("  dims:\n")
  print(x$dims)
  invisible(x)
}

#' @describeIn ldmppr_grids-class Convert a grid schedule to a data.frame.
#' @param x an object of class \code{ldmppr_grids}.
#' @param ... additional arguments (not used).
#' @export
as.data.frame.ldmppr_grids <- function(x, ...) {
  stopifnot(is_ldmppr_grids(x))
  do.call(rbind, lapply(seq_along(x$levels), function(i) {
    L <- x$levels[[i]]
    data.frame(
      level = i,
      label = x$labels[[i]],
      nx = length(L$x),
      ny = length(L$y),
      nt = length(L$t),
      x_min = min(L$x), x_max = max(L$x),
      y_min = min(L$y), y_max = max(L$y),
      t_min = min(L$t), t_max = max(L$t)
    )
  }))
}

#' @describeIn ldmppr_grids-class Number of levels in a grid schedule.
#' @param x an object of class \code{ldmppr_grids}.
#' @export
length.ldmppr_grids <- function(x) {
  stopifnot(is_ldmppr_grids(x))
  length(x$levels)
}

#' @describeIn ldmppr_grids-class Subset grid levels.
#' @param x an object of class \code{ldmppr_grids}.
#' @param i indices of levels to keep.
#' @param ... unused.
#' @export
`[.ldmppr_grids` <- function(x, i, ...) {
  stopifnot(is_ldmppr_grids(x))
  i <- as.integer(i)
  new_ldmppr_grids(
    levels = x$levels[i],
    upper_bounds = x$upper_bounds,
    labels = x$labels[i],
    include_endpoints = x$include_endpoints
  )
}

#' @describeIn ldmppr_grids-class Extract the underlying list representation.
#' @param x an object of class \code{ldmppr_grids}.
#' @param ... unused.
#' @export
as.list.ldmppr_grids <- function(x, ...) {
  stopifnot(is_ldmppr_grids(x))
  unclass(x)
}
