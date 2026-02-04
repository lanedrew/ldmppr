#' Create an optimization budget specification for estimate_process_parameters()
#'
#' \code{ldmppr_budgets()} defines per-stage optimization options (budgets) used by
#' \code{\link{estimate_process_parameters}} for NLopt via \code{\link[nloptr:nloptr]{nloptr}}.
#'
#' The returned object is an S3 class. Use \code{summary()} and \code{as.data.frame()}
#' methods (if you provide them) to inspect.
#'
#' @param global_options Optional list of NLopt options used for the global stage
#'   (only relevant when \code{strategy} uses a global optimizer).
#'   Examples: \code{list(maxeval = 2000, maxtime = 10)}.
#' @param local_budget_first_level Optional list of NLopt options used for the local stage
#'   at the first (coarsest) grid level.
#' @param local_budget_refinement_levels Optional list of NLopt options used for local refinement
#'   on subsequent (finer) grid levels in multi-resolution strategies. If \code{NULL}, the
#'   estimator will fall back to \code{local_budget_first_level}.
#'
#' @return An object of class \code{"ldmppr_budgets"}.
#' @export
#'
#' @seealso \code{\link{ldmppr_grids-class}} for methods and details.
#'
#' @examples
#' b <- ldmppr_budgets(
#'   global_options = list(maxeval = 150),
#'   local_budget_first_level = list(maxeval = 300, xtol_rel = 1e-5),
#'   local_budget_refinement_levels = list(maxeval = 150, xtol_rel = 1e-5)
#' )
#' b
ldmppr_budgets <- function(global_options = NULL,
                           local_budget_first_level = NULL,
                           local_budget_refinement_levels = NULL) {

  x <- new_ldmppr_budgets(
    global_options = global_options %||% list(),
    local_budget_first_level = local_budget_first_level %||% list(),
    local_budget_refinement_levels = local_budget_refinement_levels
  )

  x
}
