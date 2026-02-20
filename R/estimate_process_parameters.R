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
#'   For small problems, parallel overhead may outweigh speed gains.
#' @param num_cores number of workers to use when \code{parallel=TRUE} and
#'   \code{set_future_plan=TRUE}.
#' @param set_future_plan \code{TRUE} or \code{FALSE}. If \code{TRUE} and
#'   \code{parallel=TRUE}, set a temporary \pkg{future} plan internally and
#'   restore the previous plan on exit.
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
#' @param rescore_control controls candidate rescoring and bound-handling behavior in
#'   multi-resolution fitting. Can be either:
#'   
#'   \itemize{
#'     \item a single logical value (toggle rescoring on/off while keeping defaults), or
#'     \item a named list with any of:
#'       \code{enabled}, \code{top}, \code{objective_tol}, \code{param_tol},
#'       \code{avoid_bound_solutions}, \code{bound_eps}.
#'   }
#'
#'   Defaults are:
#'   \code{list(enabled = TRUE, top = 5L, objective_tol = 1e-6, param_tol = 0.10,}
#'   \code{avoid_bound_solutions = TRUE, bound_eps = 1e-8)}.
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
                                        rescore_control = list(
                                          enabled = TRUE,
                                          top = 5L,
                                          objective_tol = 1e-6,
                                          param_tol = 0.10,
                                          avoid_bound_solutions = TRUE,
                                          bound_eps = 1e-8
                                        ),
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

  delta_for_init <- if (delta_is_search) 1 else if (delta_is_single) delta else NULL

  # ---- consolidated rescoring/bound controls ----
  default_rescore_control <- list(
    enabled = TRUE,
    top = 5L,
    objective_tol = 1e-6,
    param_tol = 0.10,
    avoid_bound_solutions = TRUE,
    bound_eps = 1e-8
  )

  if (is.logical(rescore_control) && length(rescore_control) == 1L && !is.na(rescore_control)) {
    rc <- default_rescore_control
    rc$enabled <- isTRUE(rescore_control)
  } else if (is.list(rescore_control)) {
    unknown <- setdiff(names(rescore_control), names(default_rescore_control))
    if (length(unknown)) {
      stop("Unknown `rescore_control` field(s): ", paste(unknown, collapse = ", "), call. = FALSE)
    }
    rc <- utils::modifyList(default_rescore_control, rescore_control)
  } else {
    stop("`rescore_control` must be a single logical value or a named list.", call. = FALSE)
  }

  if (!is.logical(rc$enabled) || length(rc$enabled) != 1L || is.na(rc$enabled)) {
    stop("`rescore_control$enabled` must be TRUE/FALSE.", call. = FALSE)
  }
  if (!is.numeric(rc$top) || length(rc$top) != 1L || is.na(rc$top) || rc$top < 1) {
    stop("`rescore_control$top` must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(rc$objective_tol) || length(rc$objective_tol) != 1L || is.na(rc$objective_tol) || rc$objective_tol < 0) {
    stop("`rescore_control$objective_tol` must be a nonnegative number.", call. = FALSE)
  }
  if (!is.numeric(rc$param_tol) || length(rc$param_tol) != 1L || is.na(rc$param_tol) || rc$param_tol < 0) {
    stop("`rescore_control$param_tol` must be a nonnegative number.", call. = FALSE)
  }
  if (!is.logical(rc$avoid_bound_solutions) || length(rc$avoid_bound_solutions) != 1L || is.na(rc$avoid_bound_solutions)) {
    stop("`rescore_control$avoid_bound_solutions` must be TRUE/FALSE.", call. = FALSE)
  }
  if (!is.numeric(rc$bound_eps) || length(rc$bound_eps) != 1L || is.na(rc$bound_eps) || rc$bound_eps <= 0) {
    stop("`rescore_control$bound_eps` must be a positive number.", call. = FALSE)
  }

  rescore_enabled <- isTRUE(rc$enabled)
  rescore_top_n <- as.integer(rc$top)
  rescore_obj_tol <- as.numeric(rc$objective_tol)
  rescore_par_tol <- as.numeric(rc$param_tol)
  avoid_bound_solutions_flag <- isTRUE(rc$avoid_bound_solutions)
  bound_eps_val <- as.numeric(rc$bound_eps)

  # ---- user-friendly helpers ----
  .vstate <- new.env(parent = emptyenv())
  .vstate$header_printed <- FALSE
  .vmsg <- function(..., .indent = 0L) {
    if (!isTRUE(verbose)) return(invisible(NULL))
    if (!isTRUE(.vstate$header_printed)) {
      message("[ldmppr::estimate_process_parameters]")
      .vstate$header_printed <- TRUE
    }
    indent <- if (.indent > 0L) paste(rep("  ", .indent), collapse = "") else ""
    message(indent, paste0(..., collapse = ""))
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

  global_opts <- budgets$global_options %||% list()
  # Backward-compatible normalization for NLopt RNG option naming.
  # NLopt expects `ranseed`; legacy user configs may still pass `seed`.
  if (!is.null(global_opts$seed) && is.null(global_opts$ranseed)) {
    global_opts$ranseed <- global_opts$seed
  }
  global_opts$seed <- NULL

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
    finite_bounds <- .derive_finite_bounds(
      init = parameter_inits,
      data = data,
      upper_bounds = upper_bounds,
      delta = delta_for_init
    )
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
  .vmsg("Rescore control: enabled=", rescore_enabled,
        ", top=", rescore_top_n,
        ", objective_tol=", signif(rescore_obj_tol, 4),
        ", param_tol=", signif(rescore_par_tol, 4), .indent = 1L)
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
      if (isTRUE(rescore_enabled) && isTRUE(do_multires) && lvl > 1L) .vmsg("Rescoring candidates on this grid before optimizing.", .indent = 1L)

      t_lvl <- .tic()

      # In multires mode, rescore previous-level candidates on the current grid.
      if (isTRUE(rescore_enabled) && isTRUE(do_multires) && lvl > 1L) {
        prev <- stage_store[[lvl - 1L]]

        cand <- list()
        # best local
        cand[[length(cand) + 1L]] <- prev$best$solution
        # all local fits
        if (!is.null(prev$local_fits) && length(prev$local_fits)) {
          for (ff in prev$local_fits) cand[[length(cand) + 1L]] <- ff$solution
        }
        # global best + restarts
        if (!is.null(prev$global_fit)) cand[[length(cand) + 1L]] <- prev$global_fit$solution
        if (!is.null(prev$global_fits) && length(prev$global_fits)) {
          for (gg in prev$global_fits) cand[[length(cand) + 1L]] <- gg$solution
        }

        cand <- .safe_unique_params(cand)

        # rescore candidates on THIS grid, then keep a diverse subset near the best
        resc <- .rescore_on_grid(cand, grids_lvl, data_mat, upper_bounds)

        keep_n <- min(as.integer(rescore_top_n), length(resc$params))
        best_obj <- resc$obj[1]

        kept <- list(resc$params[[1]])
        for (i in 2:length(resc$params)) {
          if (length(kept) >= keep_n) break
          if (!is.finite(resc$obj[i])) next
          if (resc$obj[i] > best_obj + rescore_obj_tol) next
          # keep only if meaningfully different from all already kept
          rels <- vapply(kept, function(p0) .rel_l2(resc$params[[i]], p0), numeric(1))
          if (all(rels >= rescore_par_tol)) kept[[length(kept) + 1L]] <- resc$params[[i]]
        }

        # Optional local jitter around the best rescored start.
        # This preserves the carry-forward best while still exploring nearby starts.
        starts_list <- kept
        if (isTRUE(starts$local > 1L)) {
          if (.needs_finite_bounds(local_algorithm)) {
            lb_loc <- finite_bounds$lb
            ub_loc <- finite_bounds$ub
          } else {
            lb_loc <- ub_loc <- NULL
          }
          jit <- .jitter_inits(
            base_init = kept[[1]],
            n = starts$local,
            sd = starts$jitter_sd,
            seed = as.integer(starts$seed) + as.integer(lvl) - 1L,
            lb = lb_loc,
            ub = ub_loc,
            clamp = !is.null(lb_loc)
          )
          starts_list <- .safe_unique_params(c(starts_list, jit))
        }

        # run local polishing from rescored starts (+ optional jitter), take best
        pol <- .run_local_from_starts(
          starts_list = starts_list,
          grids_lvl = grids_lvl,
          data_mat = data_mat,
          upper_bounds = upper_bounds,
          local_algorithm = local_algorithm,
          local_options = local_opts,
          finite_bounds = finite_bounds
        )

        lvl_fit <- list(
          best = pol$best,
          local_fits = pol$fits,
          global_fit = NULL,
          global_fits = NULL
        )

        current_params <- lvl_fit$best$solution

      } else {

        # default path: original fitter (global optional + local multistart optional)
        lvl_fit <- .fit_one_level_sc(
          level_grids      = grids_lvl,
          data_mat         = data_mat,
          start_params     = current_params,
          upper_bounds     = upper_bounds,
          global_algorithm = global_algorithm,
          local_algorithm  = local_algorithm,
          global_options   = global_opts,
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

        # Optionally pass distinct global solutions into local polishing.
        if (isTRUE(do_global) && !is.null(lvl_fit$global_fits) && length(lvl_fit$global_fits) > 1L) {

          # Keep up to rescore_top global solutions (by objective), distinct-ish
          gobj <- vapply(lvl_fit$global_fits, function(f) f$objective, numeric(1))
          ord  <- order(gobj)
          gsol <- lapply(lvl_fit$global_fits[ord], function(f) f$solution)

          gsol <- .safe_unique_params(gsol)
          gsol <- gsol[seq_len(min(length(gsol), as.integer(rescore_top_n)))]

          # polish them locally (in addition to whatever we already did)
          pol <- .run_local_from_starts(
            starts_list = gsol,
            grids_lvl = grids_lvl,
            data_mat = data_mat,
            upper_bounds = upper_bounds,
            local_algorithm = local_algorithm,
            local_options = local_opts,
            finite_bounds = finite_bounds
          )

          # merge: if polished best beats existing best, replace
          if (is.finite(pol$best$objective) && pol$best$objective < lvl_fit$best$objective) {
            lvl_fit$best <- pol$best
          }
          # keep a trail of fits
          lvl_fit$local_fits <- c(lvl_fit$local_fits, pol$fits)
        }

        # Avoid bound-hugging solutions by nudging inward then re-polishing.
        if (isTRUE(avoid_bound_solutions_flag) && .needs_finite_bounds(local_algorithm)) {
          lb <- finite_bounds$lb
          ub <- finite_bounds$ub
          if (.is_on_bound(lvl_fit$best$solution, lb, ub, eps = bound_eps_val)) {
            .vmsg("Best solution is on/near a bound; nudging interior and re-polishing locally.", .indent = 1L)
            nudged <- .nudge_interior(lvl_fit$best$solution, lb, ub, eps = bound_eps_val)
            pol2 <- .run_local_from_starts(
              starts_list = list(nudged),
              grids_lvl = grids_lvl,
              data_mat = data_mat,
              upper_bounds = upper_bounds,
              local_algorithm = local_algorithm,
              local_options = local_opts,
              finite_bounds = finite_bounds
            )
            if (is.finite(pol2$best$objective) && pol2$best$objective <= lvl_fit$best$objective) {
              lvl_fit$best <- pol2$best
              lvl_fit$local_fits <- c(lvl_fit$local_fits, pol2$fits)
            }
          }
        }

        current_params <- lvl_fit$best$solution
      }

      stage_store[[lvl]] <- lvl_fit

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
      settings = list(
        strategy = strategy,
        global_algorithm = global_algorithm,
        local_algorithm = local_algorithm,
        starts = starts,
        rescore_control = rc,
        parallel = isTRUE(parallel),
        num_cores = if (isTRUE(parallel)) as.integer(num_cores) else 1L,
        set_future_plan = isTRUE(set_future_plan),
        refine_best_delta = isTRUE(refine_best_delta)
      ),
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
      global_options   = global_opts,
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

      # Optional rescoring during delta refinement.
      if (isTRUE(rescore_enabled)) {
        # candidates: previous refined best + coarse best
        cand <- list(final_params)
        # also include a few coarse candidates (if stored)
        coarse_best <- coarse_results[[best_idx]]
        if (!is.null(coarse_best$local_fits) && length(coarse_best$local_fits)) {
          for (ff in coarse_best$local_fits) cand[[length(cand) + 1L]] <- ff$solution
        }
        if (!is.null(coarse_best$global_fits) && length(coarse_best$global_fits)) {
          for (gg in coarse_best$global_fits) cand[[length(cand) + 1L]] <- gg$solution
        }
        cand <- .safe_unique_params(cand)

        resc <- .rescore_on_grid(cand, grids_lvl, mat_best, upper_bounds)
        keep_n <- min(as.integer(rescore_top_n), length(resc$params))
        best_obj <- resc$obj[1]
        kept <- list(resc$params[[1]])
        for (i in 2:length(resc$params)) {
          if (length(kept) >= keep_n) break
          if (resc$obj[i] > best_obj + rescore_obj_tol) next
          rels <- vapply(kept, function(p0) .rel_l2(resc$params[[i]], p0), numeric(1))
          if (all(rels >= rescore_par_tol)) kept[[length(kept) + 1L]] <- resc$params[[i]]
        }

        starts_list <- kept
        if (isTRUE(starts$local > 1L)) {
          if (.needs_finite_bounds(local_algorithm)) {
            lb_loc <- finite_bounds$lb
            ub_loc <- finite_bounds$ub
          } else {
            lb_loc <- ub_loc <- NULL
          }
          jit <- .jitter_inits(
            base_init = kept[[1]],
            n = starts$local,
            sd = starts$jitter_sd,
            seed = as.integer(starts$seed) + as.integer(lvl) - 1L,
            lb = lb_loc,
            ub = ub_loc,
            clamp = !is.null(lb_loc)
          )
          starts_list <- .safe_unique_params(c(starts_list, jit))
        }

        pol <- .run_local_from_starts(
          starts_list = starts_list,
          grids_lvl = grids_lvl,
          data_mat = mat_best,
          upper_bounds = upper_bounds,
          local_algorithm = local_algorithm,
          local_options = local_opts,
          finite_bounds = finite_bounds
        )

        lvl_fit <- list(best = pol$best, local_fits = pol$fits, global_fit = NULL, global_fits = NULL)

      } else {

        lvl_fit <- .fit_one_level_sc(
          level_grids      = grids_lvl,
          data_mat         = mat_best,
          start_params     = final_params,
          upper_bounds     = upper_bounds,
          global_algorithm = global_algorithm,
          local_algorithm  = local_algorithm,
          global_options   = global_opts,
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
      }

      refined_store[[lvl]] <- lvl_fit
      final_fit <- lvl_fit$best
      final_params <- lvl_fit$best$solution

      # avoid bound solutions in refinement too
      if (isTRUE(avoid_bound_solutions_flag) && .needs_finite_bounds(local_algorithm)) {
        lb <- finite_bounds$lb
        ub <- finite_bounds$ub
        if (.is_on_bound(final_params, lb, ub, eps = bound_eps_val)) {
          nudged <- .nudge_interior(final_params, lb, ub, eps = bound_eps_val)
          pol2 <- .run_local_from_starts(
            starts_list = list(nudged),
            grids_lvl = grids_lvl,
            data_mat = mat_best,
            upper_bounds = upper_bounds,
            local_algorithm = local_algorithm,
            local_options = local_opts,
            finite_bounds = finite_bounds
          )
          if (is.finite(pol2$best$objective) && pol2$best$objective <= final_fit$objective) {
            final_fit <- pol2$best
            final_params <- pol2$best$solution
            refined_store[[lvl]]$best <- pol2$best
            refined_store[[lvl]]$local_fits <- c(refined_store[[lvl]]$local_fits, pol2$fits)
          }
        }
      }

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
    settings = list(
      strategy = strategy,
      global_algorithm = global_algorithm,
      local_algorithm = local_algorithm,
      starts = starts,
      rescore_control = rc,
      parallel = isTRUE(parallel),
      num_cores = if (isTRUE(parallel)) as.integer(num_cores) else 1L,
      set_future_plan = isTRUE(set_future_plan),
      refine_best_delta = isTRUE(refine_best_delta)
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
