#' Optimization budget specification object
#'
#' Objects of class \code{ldmppr_budgets} define optimizer budget/options used by
#' \code{\link{estimate_process_parameters}}.
#'
#' @details
#' A \code{ldmppr_budgets} is a list with (at minimum):
#' \itemize{
#'   \item \code{global_options}: list of NLopt options for the global stage (e.g., \code{maxeval}, \code{maxtime}).
#'   \item \code{local_budget_first_level}: list of NLopt options for the local stage at the first/coarsest grid level.
#'   \item \code{local_budget_refinement_levels}: optional list of NLopt options for local refinement levels
#'     (used only when \code{estimate_process_parameters(strategy = "multires_global_local")}).
#' }
#'
#' @return
#' \describe{
#'   \item{\code{print()}}{prints a brief description of configured budgets.}
#'   \item{\code{summary()}}{returns a \code{summary.ldmppr_budgets}.}
#'   \item{\code{as.data.frame()}}{a compact table of the global/local budget entries (when present).}
#'   \item{\code{length()}}{number of available local budget stages (1 or 2).}
#'   \item{\code{[ ]}}{subset local budget stages (keeps global options).}
#'   \item{\code{as.list()}}{returns the underlying list structure.}
#' }
#'
#' @name ldmppr_budgets-class
#' @docType class
NULL


#' @describeIn ldmppr_budgets-class Print a brief summary of optimization budgets.
#' @param x an object of class \code{ldmppr_budgets}.
#' @param ... additional arguments (not used).
#' @export
print.ldmppr_budgets <- function(x, ...) {
  stopifnot(is_ldmppr_budgets(x))

  cat("<ldmppr_budgets>\n")

  g <- x$global_options %||% list()
  l1 <- x$local_budget_first_level %||% list()
  l2 <- x$local_budget_refinement_levels

  cat("  global_options:\n")
  cat("    - maxeval: ", g$maxeval %||% NA, "\n", sep = "")
  cat("    - maxtime: ", g$maxtime %||% NA, "\n", sep = "")

  cat("  local_budget_first_level:\n")
  cat("    - maxeval:  ", l1$maxeval  %||% NA, "\n", sep = "")
  cat("    - maxtime:  ", l1$maxtime  %||% NA, "\n", sep = "")
  cat("    - xtol_rel: ", l1$xtol_rel %||% NA, "\n", sep = "")

  if (is.null(l2)) {
    cat("  local_budget_refinement_levels: <NULL>\n")
  } else {
    cat("  local_budget_refinement_levels:\n")
    cat("    - maxeval:  ", l2$maxeval  %||% NA, "\n", sep = "")
    cat("    - maxtime:  ", l2$maxtime  %||% NA, "\n", sep = "")
    cat("    - xtol_rel: ", l2$xtol_rel %||% NA, "\n", sep = "")
  }

  invisible(x)
}


#' @describeIn ldmppr_budgets-class Summarize an optimization budget specification.
#' @param object an object of class \code{ldmppr_budgets}.
#' @param ... additional arguments (not used).
#' @export
summary.ldmppr_budgets <- function(object, ...) {
  stopifnot(is_ldmppr_budgets(object))

  out <- list(
    n_local_stages = length(object),
    has_refinement = !is.null(object$local_budget_refinement_levels),
    table = as.data.frame(object)
  )
  class(out) <- "summary.ldmppr_budgets"
  out
}

#' @describeIn ldmppr_budgets-class Print a summary produced by \code{summary.ldmppr_budgets()}.
#' @param x an object of class \code{summary.ldmppr_budgets}.
#' @param ... additional arguments (not used).
#' @export
print.summary.ldmppr_budgets <- function(x, ...) {
  cat("<summary: ldmppr_budgets>\n")
  cat("  local stages: ", x$n_local_stages, "\n", sep = "")
  cat("  has refinement budgets: ", if (isTRUE(x$has_refinement)) "yes" else "no", "\n", sep = "")
  print(x$table)
  invisible(x)
}


#' @describeIn ldmppr_budgets-class Convert budgets to a data.frame.
#' @param x an object of class \code{ldmppr_budgets}.
#' @param ... additional arguments (not used).
#' @export
as.data.frame.ldmppr_budgets <- function(x, ...) {
  stopifnot(is_ldmppr_budgets(x))

  grab <- function(lst, nm) if (!is.null(lst) && !is.null(lst[[nm]])) lst[[nm]] else NA

  rows <- list(
    data.frame(
      stage = "global",
      maxeval = grab(x$global_options, "maxeval"),
      maxtime = grab(x$global_options, "maxtime"),
      xtol_rel = grab(x$global_options, "xtol_rel"),
      stringsAsFactors = FALSE
    ),
    data.frame(
      stage = "local_first_level",
      maxeval = grab(x$local_budget_first_level, "maxeval"),
      maxtime = grab(x$local_budget_first_level, "maxtime"),
      xtol_rel = grab(x$local_budget_first_level, "xtol_rel"),
      stringsAsFactors = FALSE
    )
  )

  if (!is.null(x$local_budget_refinement_levels)) {
    rows[[length(rows) + 1L]] <- data.frame(
      stage = "local_refinement_levels",
      maxeval = grab(x$local_budget_refinement_levels, "maxeval"),
      maxtime = grab(x$local_budget_refinement_levels, "maxtime"),
      xtol_rel = grab(x$local_budget_refinement_levels, "xtol_rel"),
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}


#' @describeIn ldmppr_budgets-class Number of available local budget stages.
#' @param x an object of class \code{ldmppr_budgets}.
#' @export
length.ldmppr_budgets <- function(x) {
  stopifnot(is_ldmppr_budgets(x))
  if (is.null(x$local_budget_refinement_levels)) 1L else 2L
}


#' @describeIn ldmppr_budgets-class Subset local budget stages (keeps global options).
#' @param x an object of class \code{ldmppr_budgets}.
#' @param i indices of local stages to keep: 1 = first level, 2 = refinement levels.
#' @param ... unused.
#' @export
`[.ldmppr_budgets` <- function(x, i, ...) {
  stopifnot(is_ldmppr_budgets(x))

  i <- as.integer(i)
  i <- i[!is.na(i)]
  if (!length(i)) stop("`i` must select 1 and/or 2.", call. = FALSE)
  if (any(!i %in% c(1L, 2L))) stop("Valid indices are 1 (first level) and 2 (refinement levels).", call. = FALSE)

  keep_first <- 1L %in% i
  keep_ref   <- 2L %in% i

  if (keep_ref && is.null(x$local_budget_refinement_levels)) {
    stop("No refinement budgets are present to subset (local_budget_refinement_levels is NULL).", call. = FALSE)
  }

  new_ldmppr_budgets(
    global_options = x$global_options %||% list(),
    local_budget_first_level = if (keep_first) (x$local_budget_first_level %||% list()) else list(),
    local_budget_refinement_levels = if (keep_ref) x$local_budget_refinement_levels else NULL
  )
}


#' @describeIn ldmppr_budgets-class Extract the underlying list representation.
#' @param x an object of class \code{ldmppr_budgets}.
#' @param ... unused.
#' @export
as.list.ldmppr_budgets <- function(x, ...) {
  stopifnot(is_ldmppr_budgets(x))
  unclass(x)
}
