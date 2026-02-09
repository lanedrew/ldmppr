#' Estimate point process parameters using log-likelihood maximization
#'
#' Estimate spatio-temporal point process parameters by maximizing the (approximate)
#' full log-likelihood using \code{\link[nloptr:nloptr]{nloptr}}.
#'
#' For the self-correcting process, arrival times must lie on \eqn{(0,1)} and can be
#' supplied directly in \code{data} as \code{time}, or constructed from \code{size}
#' via the gentle-decay (power-law) mapping \code{\link{power_law_mapping}} using
#' \code{delta}. When \code{delta} is a vector, the model is fit for each candidate
#' value and the best objective is selected.
#'
#' This function supports multi-resolution estimation via a \code{\link{ldmppr_grids}}
#' schedule. If multiple grid levels are provided, the coarsest level may use a global
#' optimizer followed by local refinement, and subsequent levels run local refinement only.
#'
#' @param data a data.frame or matrix. Must contain either columns \code{(time, x, y)}
#'   or \code{(x, y, size)}. If a matrix is provided without time, it must have
#'   column names \code{c("x","y","size")}.
#' @param process type of process used (currently supports \code{"self_correcting"}).
#' @param grids a \code{\link{ldmppr_grids}} object specifying the integration grid schedule
#'   (single-level or multi-resolution). The integration bounds are taken from
#'   \code{grids$upper_bounds}.
#' @param budgets a \code{\link{ldmppr_budgets}} object controlling optimizer options
#'   for the global stage and local stages (first level vs refinement levels).
#' @param parameter_inits (optional) numeric vector of length 8 giving initialization values
#'   for the model parameters. If \code{NULL}, defaults are derived from \code{data} and
#'   \code{grids$upper_bounds}.
#' @param delta (optional) numeric scalar or vector. Used only when \code{data} does not contain
#'   \code{time} (i.e., data has \code{(x,y,size)}).
#'   \itemize{
#'     \item If \code{length(delta) == 1}, fits the model once using \code{power_law_mapping(size, delta)}.
#'     \item If \code{length(delta) > 1}, performs a delta-search by fitting the model for each candidate value
#'       and selecting the best objective. If \code{refine_best_delta = TRUE} and multiple grid levels are used,
#'       the best delta is refined on the remaining (finer) grid levels.
#'   }
#'   If \code{data} already contains \code{time}, \code{delta} is ignored when \code{length(delta)==1}
#'   and is an error when \code{length(delta)>1}.
#' @param parallel \code{TRUE} or \code{FALSE} specifying furrr/future to parallelize either:
#'   (a) over candidate \code{delta} values (when \code{length(delta) > 1}), and/or
#'   (b) over local multi-start initializations (when \code{starts$local > 1}), and/or
#'   (c) over global restarts (when \code{starts$global > 1}).
#' @param num_cores number of workers to use when \code{set_future_plan = TRUE}.
#' @param set_future_plan \code{TRUE} or \code{FALSE}, temporarily sets
#'   \code{future::plan(multisession, workers = num_cores)} and restores the original plan on exit.
#' @param strategy character string specifying the estimation strategy:
#'  \itemize{
#'    \item \code{"local"}: local optimization only (single-level or multi-level polish).
#'    \item \code{"global_local"}: global optimization then local polish (single grid level).
#'    \item \code{"multires_global_local"}: multi-resolution (coarsest uses global+local; refinements use local only).
#'  }
#' @param global_algorithm,local_algorithm NLopt algorithms to use for
#'   the global and local optimization stages, respectively.
#' @param starts a list controlling restart and jitter behavior:
#'   \itemize{
#'     \item \code{global}: integer, number of global restarts at the first/coarsest level (default 1).
#'     \item \code{local}: integer, number of local multi-starts per level (default 1).
#'     \item \code{jitter_sd}: numeric SD for jittering (default 0.35).
#'     \item \code{seed}: integer base seed (default 1).
#'   }
#' @param finite_bounds (optional) list with components \code{lb} and \code{ub} giving finite lower and
#'   upper bounds for all 8 parameters. If \code{NULL}, bounds are derived from \code{parameter_inits}.
#'   Global algorithms and select local algorithms in NLopt require finite bounds.
#' @param refine_best_delta \code{TRUE} or \code{FALSE}. If \code{TRUE} and \code{length(delta) > 1}, performs refinement
#'   of the best \code{delta} across additional grid levels (if available).
#' @param verbose \code{TRUE} or \code{FALSE} indicating whether to show progress of model estimation.
#'
#' @return an object of class \code{"ldmppr_fit"} containing the best \code{nloptr} fit and
#'   (optionally) stored fits from global restarts and/or a delta search.
#'
#' @references
#' Møller, J., Ghorbani, M., & Rubak, E. (2016). Mechanistic spatio-temporal point process models
#' for marked point processes, with a view to forest stand data. \emph{Biometrics}, 72(3), 687-696.
#' \doi{10.1111/biom.12466}.
#'
#' @examples
#' # Load example data
#' data(small_example_data)
#'
#' # Define grids and budgets
#' ub <- c(1, 25, 25)
#' g  <- ldmppr_grids(upper_bounds = ub, levels = list(c(10,10,10)))
#' b  <- ldmppr_budgets(
#'   global_options = list(maxeval = 150),
#'   local_budget_first_level = list(maxeval = 50, xtol_rel = 1e-2),
#'   local_budget_refinement_levels = list(maxeval = 25, xtol_rel = 1e-2)
#' )
#'
#' # Estimate parameters using a single delta value
#' fit <- estimate_process_parameters(
#'   data = small_example_data,
#'   grids = g,
#'   budgets = b,
#'   delta = 1,
#'   strategy = "global_local",
#'   global_algorithm = "NLOPT_GN_CRS2_LM",
#'   local_algorithm  = "NLOPT_LN_BOBYQA",
#'   starts = list(global = 2, local = 2, jitter_sd = 0.25, seed = 1),
#'   verbose = TRUE
#' )
#' coef(fit)
#' logLik(fit)
#'
#' \donttest{
#' # Estimate parameters using multiple delta values (delta search)
#' g2 <- ldmppr_grids(upper_bounds = ub, levels = list(c(8,8,8), c(12,12,12)))
#' fit_delta <- estimate_process_parameters(
#'   data = small_example_data, # x,y,size
#'   grids = g2,
#'   budgets = b,
#'   delta = c(0.35, 0.5, 0.65, 0.9, 1.0),
#'   parallel = TRUE,
#'   set_future_plan = TRUE,
#'   num_cores = 2,
#'   strategy = "multires_global_local",
#'   starts = list(local = 1),
#'   refine_best_delta = FALSE,
#'   verbose = FALSE
#' )
#' plot(fit_delta)
#' }
#'
#' @export
estimate_process_parameters <- function(data,
                                        process = c("self_correcting"),
                                        grids,
                                        budgets,
                                        parameter_inits = NULL,
                                        delta = NULL,
                                        parallel = FALSE,
                                        num_cores = max(1L, parallel::detectCores() - 1L),
                                        set_future_plan = FALSE,
                                        strategy = c("local", "global_local", "multires_global_local"),
                                        global_algorithm = "NLOPT_GN_CRS2_LM",
                                        local_algorithm  = "NLOPT_LN_BOBYQA",
                                        starts = list(global = 1L, local = 1L, jitter_sd = 0.35, seed = 1L),
                                        finite_bounds = NULL,
                                        refine_best_delta = TRUE,
                                        verbose = TRUE) {

  process  <- match.arg(process)
  strategy <- match.arg(strategy)

  if (process != "self_correcting") {
    stop("Only process='self_correcting' is currently implemented.", call. = FALSE)
  }

  # ---- delta validation ----
  if (!is.null(delta) && !is.numeric(delta)) {
    stop("`delta` must be NULL or numeric.", call. = FALSE)
  }
  if (!is.null(delta) && (anyNA(delta) || any(!is.finite(delta)))) {
    stop("`delta` must be finite and non-missing.", call. = FALSE)
  }

  delta_is_search <- !is.null(delta) && length(delta) > 1L
  delta_is_single <- !is.null(delta) && length(delta) == 1L

  # ---- user-friendly helpers ----
  .vmsg <- function(..., .indent = 0L) {
    if (!isTRUE(verbose)) return(invisible(NULL))
    prefix <- if (.indent > 0L) paste(rep("  ", .indent), collapse = "") else ""
    message(prefix, paste0(..., collapse = ""))
    invisible(NULL)
  }
  .fmt_sec <- function(sec) {
    sec <- as.numeric(sec)
    if (!is.finite(sec)) return("NA")
    if (sec < 60) return(sprintf("%.1fs", sec))
    if (sec < 3600) return(sprintf("%.1fmin", sec / 60))
    sprintf("%.2fh", sec / 3600)
  }
  .tic <- function() proc.time()[["elapsed"]]
  .toc <- function(t0) proc.time()[["elapsed"]] - t0

  # ---- coerce + validate user objects ----
  grids   <- as_ldmppr_grids(grids)
  budgets <- as_ldmppr_budgets(budgets)

  upper_bounds <- grids$upper_bounds

  starts <- .normalize_epp_starts(starts)

  .validate_epp_inputs(
    data = data,
    grids = grids,
    budgets = budgets,
    parameter_inits = parameter_inits,
    delta = delta,
    strategy = strategy,
    global_algorithm = global_algorithm,
    local_algorithm = local_algorithm,
    starts = starts
  )

  # ---- derive default inits if needed ----
  if (is.null(parameter_inits)) {
    delta_for_init <- if (delta_is_search) 1 else if (delta_is_single) delta else NULL
    if (delta_is_search) {
      .vmsg("Note: `delta` has multiple candidates; deriving default initialization using delta=1 for stability.", .indent = 1L)
    }
    parameter_inits <- .default_parameter_inits_sc(
      data = data,
      upper_bounds = upper_bounds,
      delta = delta_for_init
    )
    .vmsg("Using default starting values (parameter_inits) since none were provided.")
    .vmsg("Initial values: ", paste(signif(parameter_inits, 4), collapse = ", "), .indent = 1L)
  }

  # ---- bounds ----
  if (is.null(finite_bounds)) {
    finite_bounds <- .derive_finite_bounds(parameter_inits)
  } else {
    .validate_finite_bounds(finite_bounds)
  }

  # ---- will we parallelize? ----
  will_parallelize <- isTRUE(parallel) && (
    delta_is_search || starts$local > 1L || starts$global > 1L
  )
  max_workers <- if (delta_is_search) length(delta) else max(starts$local, starts$global)

  # ---- maybe set future plan ----
  original_plan <- .maybe_set_future_plan(
    set_future_plan = set_future_plan,
    will_parallelize = will_parallelize,
    num_cores = num_cores,
    max_workers = max_workers,
    verbose = FALSE
  )
  if (!is.null(original_plan)) on.exit(future::plan(original_plan), add = TRUE)

  # ---- progress summary upfront ----
  .vmsg("Estimating self-correcting process parameters")
  .vmsg("Strategy: ", strategy, if (delta_is_search) " (delta search)" else "", .indent = 1L)

  if (delta_is_search) .vmsg("Delta candidates: ", length(delta), .indent = 1L)
  if (!delta_is_search) .vmsg("Delta: ", if (delta_is_single) delta else "already in data (time provided)", .indent = 1L)

  .vmsg("Grids: ", length(grids), " level(s)", .indent = 1L)
  .vmsg("Local optimizer: ", local_algorithm, .indent = 1L)
  if (strategy %in% c("global_local", "multires_global_local")) .vmsg("Global optimizer: ", global_algorithm, .indent = 1L)
  .vmsg("Starts: global=", starts$global, ", local=", starts$local,
        ", jitter_sd=", starts$jitter_sd, ", seed=", starts$seed, .indent = 1L)

  if (will_parallelize) {
    .vmsg("Parallel: on", .indent = 1L)
    .vmsg("Requested workers: ", as.integer(num_cores), .indent = 2L)
  } else {
    .vmsg("Parallel: off", .indent = 1L)
  }

  t_all <- .tic()

  # -------------------------------------------------------------------
  # MODE A: single fit
  # -------------------------------------------------------------------
  if (!delta_is_search) {

    .vmsg("Step 1/2: Preparing data and objective function...")
    t_prep <- .tic()

    data_mat <- .build_sc_matrix(data, delta = if (delta_is_single) delta else NULL)

    # ---- safety: enforce sorted times (C++ likelihood assumes sorted) ----
    if (nrow(data_mat) > 1L) {
      tt <- data_mat[, 1]
      if (is.unsorted(tt, strictly = FALSE)) {
        o <- order(tt)
        data_mat <- data_mat[o, , drop = FALSE]
      }
    }

    data_orig <- attr(data_mat, "ldmppr_original")
    if (is.null(data_orig)) data_orig <- if (is.data.frame(data)) data else as.data.frame(data_mat)

    .vmsg("Prepared ", nrow(data_mat), " points.", .indent = 1L)
    .vmsg("Done in ", .fmt_sec(.toc(t_prep)), ".", .indent = 1L)

    do_global_first <- (strategy %in% c("global_local", "multires_global_local"))
    do_multires     <- (strategy == "multires_global_local")
    n_levels        <- length(grids$levels)

    current_params <- parameter_inits
    stage_store <- list()

    .vmsg("Step 2/2: Optimizing parameters...")

    for (lvl in seq_len(n_levels)) {
      grids_lvl <- grids$levels[[lvl]]
      do_global <- do_global_first && (lvl == 1L)
      do_multistart <- (starts$local > 1L)

      local_opts <- if (lvl == 1L) budgets$local_budget_first_level else budgets$local_budget_refinement_levels
      if (is.null(local_opts)) local_opts <- budgets$local_budget_first_level %||% list()

      nx <- length(grids_lvl$x); ny <- length(grids_lvl$y); nt <- length(grids_lvl$t)
      lvl_label <- if (do_multires) {
        paste0("Level ", lvl, " of ", n_levels, " (grid ", nx, "x", ny, "x", nt, ")")
      } else {
        paste0("Single level (grid ", nx, "x", ny, "x", nt, ")")
      }

      .vmsg(lvl_label)
      if (do_global) .vmsg("Global search: ", starts$global, " restart(s), then local refinement.", .indent = 1L)
      if (!do_global) .vmsg("Local refinement only.", .indent = 1L)
      if (do_multistart) .vmsg("Local multi-start: ", starts$local, " start(s).", .indent = 1L)

      t_lvl <- .tic()

      lvl_fit <- .fit_one_level_sc(
        level_grids      = grids_lvl,
        data_mat         = data_mat,
        start_params     = current_params,
        upper_bounds     = upper_bounds,
        global_algorithm = global_algorithm,
        local_algorithm  = local_algorithm,
        global_options   = budgets$global_options %||% list(),
        local_options    = local_opts,
        finite_bounds    = finite_bounds,
        do_global        = do_global,
        do_multistart    = do_multistart,
        n_starts         = starts$local,
        jitter_sd        = starts$jitter_sd,
        seed             = starts$seed,
        global_n_starts  = starts$global,
        worker_parallel  = isTRUE(parallel) && (starts$local > 1L || starts$global > 1L),
        worker_verbose   = FALSE
      )

      stage_store[[lvl]] <- lvl_fit
      current_params <- lvl_fit$best$solution

      .vmsg("Completed in ", .fmt_sec(.toc(t_lvl)), ".", .indent = 1L)
      .vmsg("Best objective: ", signif(lvl_fit$best$objective, 8), .indent = 1L)

      if (!do_multires) break
    }

    secs <- .toc(t_all)
    final_level <- length(stage_store)

    final_grid <- list(
      x_grid = grids$levels[[final_level]]$x,
      y_grid = grids$levels[[final_level]]$y,
      t_grid = grids$levels[[final_level]]$t,
      upper_bounds = upper_bounds
    )

    # delta used is meaningful only if time was constructed from size
    delta_used <- attr(data_mat, "ldmppr_delta")

    # fallback: if we *know* we mapped from size using a single delta, record it
    if (is.null(delta_used) || is.na(delta_used)) {
      # detect “mapped-from-size” path
      data_has_time <- is.data.frame(data) && all(c("time","x","y") %in% names(data))
      if (!data_has_time && delta_is_single) {
        delta_used <- as.numeric(delta)
      } else {
        delta_used <- NA_real_
      }
    }

    fit_obj <- new_ldmppr_fit(
      process = "self_correcting",
      fit = stage_store[[final_level]]$best,
      fits = stage_store,
      mapping = list(delta = delta_used),
      grid = final_grid,
      data_summary = list(n = nrow(data_mat)),
      data = data_mat,
      data_original = data_orig,
      engine = "nloptr",
      call = match.call(),
      timing = list(seconds = secs)
    )

    .vmsg("Finished. Total time: ", .fmt_sec(secs), ".")
    return(fit_obj)
  }

  # -------------------------------------------------------------------
  # MODE B: delta search
  # -------------------------------------------------------------------
  if (is.data.frame(data) && all(c("time", "x", "y") %in% names(data))) {
    stop("Delta search is only valid when `data` does NOT contain time (i.e., has x,y,size).", call. = FALSE)
  }
  if (is.matrix(data)) {
    cn <- colnames(data)
    ok <- !is.null(cn) && all(c("x", "y", "size") %in% cn)
    if (!ok) stop("For delta search with matrix input, `data` must have colnames including 'x', 'y', 'size'.", call. = FALSE)
  }

  .vmsg("Step 1/3: Coarse search over delta values (", length(delta), " candidates)...")
  t_coarse <- .tic()

  coarse_grids <- grids$levels[[1]]
  do_global_coarse <- (strategy %in% c("global_local", "multires_global_local"))

  fit_delta_coarse <- function(d) {
    .fit_delta_coarse_sc(
      delta            = d,
      data             = data,
      parameter_inits  = parameter_inits,
      coarse_grids     = coarse_grids,
      upper_bounds     = upper_bounds,
      global_algorithm = global_algorithm,
      local_algorithm  = local_algorithm,
      global_options   = budgets$global_options %||% list(),
      local_options    = budgets$local_budget_first_level %||% list(),
      finite_bounds    = finite_bounds,
      do_global_coarse = isTRUE(do_global_coarse),
      global_n_starts  = starts$global,
      seed             = starts$seed,
      jitter_sd        = starts$jitter_sd
    )
  }

  if (isTRUE(parallel)) {
    if (!requireNamespace("furrr", quietly = TRUE) || !requireNamespace("future", quietly = TRUE)) {
      stop("Parallel delta search requires packages 'future' and 'furrr'.", call. = FALSE)
    }
    coarse_results <- furrr::future_map(
      as.list(delta),
      fit_delta_coarse,
      .options = furrr::furrr_options(seed = TRUE)
    )
  } else {
    coarse_results <- lapply(as.list(delta), fit_delta_coarse)
  }

  coarse_obj <- vapply(coarse_results, function(r) r$best$objective, numeric(1))
  best_idx <- which.min(coarse_obj)
  best_delta <- delta[best_idx]

  .vmsg("Coarse delta search done in ", .fmt_sec(.toc(t_coarse)), ".")
  .vmsg("Best delta: ", best_delta, " (objective=", signif(coarse_obj[best_idx], 8), ").", .indent = 1L)

  .vmsg("Step 2/3: Building data under best delta and selecting best coarse fit...")
  mat_best  <- .build_sc_matrix(data, delta = best_delta)

  # ---- safety: enforce sorted times (C++ likelihood assumes sorted) ----
  if (nrow(mat_best) > 1L) {
    tt <- mat_best[, 1]
    if (is.unsorted(tt, strictly = FALSE)) {
      o <- order(tt)
      mat_best <- mat_best[o, , drop = FALSE]
    }
  }

  orig_best <- attr(mat_best, "ldmppr_original")
  if (is.null(orig_best)) orig_best <- if (is.data.frame(data)) data else as.data.frame(mat_best)

  refined_store <- NULL
  final_fit <- coarse_results[[best_idx]]$best
  final_params <- final_fit$solution

  if (isTRUE(refine_best_delta) && length(grids$levels) > 1L) {
    .vmsg("Step 3/3: Refining best delta on finer grids (", length(grids$levels) - 1L, " more level(s))...")
    t_ref <- .tic()

    refined_store <- list()
    for (lvl in 2:length(grids$levels)) {
      grids_lvl <- grids$levels[[lvl]]
      local_opts <- budgets$local_budget_refinement_levels %||% budgets$local_budget_first_level %||% list()

      nx <- length(grids_lvl$x); ny <- length(grids_lvl$y); nt <- length(grids_lvl$t)
      .vmsg("Refinement level ", lvl, " of ", length(grids$levels),
            " (grid ", nx, "x", ny, "x", nt, ")", .indent = 1L)

      t_lvl <- .tic()
      lvl_fit <- .fit_one_level_sc(
        level_grids      = grids_lvl,
        data_mat         = mat_best,
        start_params     = final_params,
        upper_bounds     = upper_bounds,
        global_algorithm = global_algorithm,
        local_algorithm  = local_algorithm,
        global_options   = budgets$global_options %||% list(),
        local_options    = local_opts,
        finite_bounds    = finite_bounds,
        do_global        = FALSE,
        do_multistart    = (starts$local > 1L),
        n_starts         = starts$local,
        jitter_sd        = starts$jitter_sd,
        seed             = starts$seed,
        worker_parallel  = FALSE,
        worker_verbose   = FALSE
      )

      refined_store[[lvl]] <- lvl_fit
      final_fit <- lvl_fit$best
      final_params <- lvl_fit$best$solution

      .vmsg("Completed in ", .fmt_sec(.toc(t_lvl)), "; objective=", signif(final_fit$objective, 8), ".", .indent = 2L)
    }

    .vmsg("Refinement done in ", .fmt_sec(.toc(t_ref)), ".")
  } else {
    .vmsg("Step 3/3: Skipping refinement (refine_best_delta=FALSE or only one grid level).")
  }

  secs <- .toc(t_all)

  final_level <- length(grids$levels)
  final_grid <- list(
    x_grid = grids$levels[[final_level]]$x,
    y_grid = grids$levels[[final_level]]$y,
    t_grid = grids$levels[[final_level]]$t,
    upper_bounds = upper_bounds
  )

  fit_obj <- new_ldmppr_fit(
    process = "self_correcting",
    fit = final_fit,
    fits = list(coarse = coarse_results, refined_best = refined_store),
    mapping = list(
      delta = best_delta,
      delta_values = delta,
      chosen_index = best_idx,
      objectives = coarse_obj,
      refined = isTRUE(refine_best_delta) && length(grids$levels) > 1L
    ),
    grid = final_grid,
    data_summary = list(n = nrow(mat_best)),
    data = mat_best,
    data_original = orig_best,
    engine = "nloptr",
    call = match.call(),
    timing = list(seconds = secs)
  )

  .vmsg("Finished. Total time: ", .fmt_sec(secs), ".")
  fit_obj
}
