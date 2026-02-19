#' Check the fit of an estimated model using global envelope tests
#'
#' @description
#' Performs global envelope tests for nonparametric L, F, G, J, E, and V summary
#' functions (\code{\link[spatstat:spatstat]{spatstat}}/\code{\link[GET:GET]{GET}}).
#' These tests assess goodness-of-fit of the estimated model relative to a reference marked point pattern.
#' The reference marked point pattern can be supplied directly via \code{reference_data} (a marked \code{ppp} object),
#' or derived internally from a \code{ldmppr_fit} object.
#'
#' @param reference_data (optional) a marked \code{ppp} object for the reference dataset.
#'   If \code{NULL}, the reference pattern is derived from \code{process_fit} when
#'   \code{process_fit} is an \code{ldmppr_fit} and contains \code{data_original}
#'   (preferred) or \code{data} with columns \code{(x,y,size)}.
#' @param t_min minimum value for time.
#' @param t_max maximum value for time.
#' @param process type of process used (currently supports \code{"self_correcting"}).
#' @param process_fit either an \code{ldmppr_fit} object (from \code{estimate_process_parameters})
#'   or a numeric vector of length 8 giving the process parameters.
#' @param anchor_point (optional) vector of (x,y) coordinates of the point to condition on.
#'   If \code{NULL}, inferred from the reference data (largest mark if available) or from
#'   \code{ldmppr_fit}.
#' @param raster_list (optional) list of raster objects used for mark prediction.
#'   Required when \code{mark_mode='mark_model'} unless rasters are stored in \code{mark_model}.
#' @param scaled_rasters \code{TRUE} or \code{FALSE} indicating whether rasters are already scaled.
#'   Ignored when \code{mark_mode='time_to_size'}.
#' @param mark_model a mark model object used when \code{mark_mode='mark_model'}.
#'   May be an \code{ldmppr_mark_model}, \code{model_fit}, or \code{workflow}.
#' @param xy_bounds (optional) vector of bounds as \code{c(a_x, b_x, a_y, b_y)}. If \code{NULL}, will be
#'   inferred from \code{reference_data}'s window when \code{reference_data} is provided,
#'   otherwise from \code{ldmppr_fit} with lower bounds assumed to be 0.
#' @param include_comp_inds \code{TRUE} or \code{FALSE} indicating whether to compute competition indices.
#' @param competition_radius positive numeric distance used when \code{include_comp_inds = TRUE}.
#' @param thinning \code{TRUE} or \code{FALSE} indicating whether to use the thinned simulated values.
#' @param edge_correction type of edge correction to apply (\code{"none"} or \code{"toroidal"}).
#' @param n_sim number of simulated datasets to generate.
#' @param save_sims \code{TRUE} or \code{FALSE} indicating whether to save and return the simulated metrics.
#' @param verbose \code{TRUE} or \code{FALSE} indicating whether to show progress of model checking.
#'   When \code{TRUE}, progress is reported via \pkg{progressr} (if available) and is compatible with parallel execution.
#' @param seed integer seed for reproducibility.
#' @param parallel \code{TRUE} or \code{FALSE}. If \code{TRUE}, simulations run in parallel via \pkg{furrr}/\pkg{future}.
#'   For small simulation sizes, parallel overhead may outweigh speed gains.
#' @param num_cores number of workers to use when \code{parallel=TRUE}. Defaults to one fewer than detected cores.
#' @param set_future_plan \code{TRUE} or \code{FALSE}. If \code{TRUE} and \code{parallel=TRUE},
#'   set a temporary \pkg{future} plan internally and restore the previous plan on exit.
#' @param mark_mode mark generation mode: \code{"mark_model"} uses \code{predict()} on a mark model,
#'   while \code{"time_to_size"} maps simulated times back to sizes via \code{delta}.
#' @param fg_correction correction used for F/G/J summaries (\code{"km"} or \code{"rs"}).
#' @param max_attempts maximum number of simulation attempts when rejection occurs
#'   due to non-finite summaries.
#'
#' @return an object of class \code{"ldmppr_model_check"}.
#'
#' @details
#' This function relies on the \code{\link[spatstat:spatstat]{spatstat}} package for the calculation of the point pattern metrics
#' and the \code{\link[GET:GET]{GET}} package for the global envelope tests. The L, F, G, J, E, and V functions are a collection of
#' non-parametric summary statistics that describe the spatial distribution of points and marks in a point pattern.
#' See the documentation for \code{\link[spatstat.explore:Lest]{Lest()}}, \code{\link[spatstat.explore:Fest]{Fest()}}, \code{\link[spatstat.explore:Gest]{Gest()}},
#' \code{\link[spatstat.explore:Jest]{Jest()}}, \code{\link[spatstat.explore:Emark]{Emark()}}, and \code{\link[spatstat.explore:Vmark]{Vmark()}} for more information.
#' Also, see the \code{\link[GET:global_envelope_test]{global_envelope_test()}} function for more information on the global envelope tests.
#'
#' @references
#' Baddeley, A., Rubak, E., & Turner, R. (2015). *Spatial Point Patterns:
#' Methodology and Applications with R*. Chapman and Hall/CRC Press, London.
#' ISBN 9781482210200. Available at:
#' \url{https://www.routledge.com/Spatial-Point-Patterns-Methodology-and-Applications-with-R/Baddeley-Rubak-Turner/p/book/9781482210200}.
#'
#' Myllymäki, M., & Mrkvička, T. (2023). GET: Global envelopes in R.
#' \emph{arXiv:1911.06583 [stat.ME]}. \doi{10.48550/arXiv.1911.06583}.
#'
#' @examples
#' # Note: The example below is provided for illustrative purposes and may take some time to run.
#' \donttest{
#' data(small_example_data)
#'
#' file_path <- system.file("extdata", "example_mark_model.rds", package = "ldmppr")
#' mark_model <- load_mark_model(file_path)
#'
#' raster_paths <- list.files(system.file("extdata", package = "ldmppr"),
#'                            pattern = "\\.tif$", full.names = TRUE)
#' raster_paths <- raster_paths[!grepl("_med\\.tif$", raster_paths)]
#' rasters <- lapply(raster_paths, terra::rast)
#' scaled_raster_list <- scale_rasters(rasters)
#'
#' reference_data <- generate_mpp(
#'   locations = small_example_data[, c("x", "y")],
#'   marks = small_example_data$size,
#'   xy_bounds = c(0, 25, 0, 25)
#' )
#'
#' estimated_parameters <- c(
#'   0.05167978, 8.20702166, 0.02199940, 2.63236890,
#'   1.82729512, 0.65330061, 0.86666748, 0.04681878
#' )
#'
#' # Keep parallel=FALSE in examples to avoid setup overhead.
#' example_model_fit <- check_model_fit(
#'   reference_data = reference_data,
#'   t_min = 0,
#'   t_max = 1,
#'   process = "self_correcting",
#'   process_fit = estimated_parameters,
#'   raster_list = scaled_raster_list,
#'   scaled_rasters = TRUE,
#'   mark_model = mark_model,
#'   xy_bounds = c(0, 25, 0, 25),
#'   include_comp_inds = TRUE,
#'   thinning = TRUE,
#'   edge_correction = "none",
#'   competition_radius = 10,
#'   n_sim = 100,
#'   save_sims = FALSE,
#'   verbose = TRUE,
#'   seed = 90210,
#'   parallel = FALSE
#' )
#'
#' plot(example_model_fit, which = 'combined')
#' }
#' @export
check_model_fit <- function(reference_data = NULL,
                            t_min = 0,
                            t_max = 1,
                            process = c("self_correcting"),
                            process_fit = NULL,
                            anchor_point = NULL,
                            raster_list = NULL,
                            scaled_rasters = FALSE,
                            mark_model = NULL,
                            xy_bounds = NULL,
                            include_comp_inds = FALSE,
                            thinning = TRUE,
                            edge_correction = "none",
                            competition_radius = 15,
                            n_sim = 2500,
                            save_sims = TRUE,
                            verbose = TRUE,
                            seed = 0,
                            parallel = FALSE,
                            num_cores = max(1L, parallel::detectCores() - 1L),
                            set_future_plan = FALSE,
                            mark_mode = c("mark_model", "time_to_size"),
                            fg_correction = c("km", "rs"),
                            max_attempts = NULL) {

  process <- match.arg(process)
  mark_mode <- match.arg(mark_mode)
  fg_correction <- match.arg(fg_correction)

  # ---- argument checks ----
  if (is.null(t_min) || t_min < 0 || t_min >= t_max) stop("Provide a value >= 0 and < t_max for `t_min`.", call. = FALSE)
  if (is.null(t_max) || t_max <= t_min || t_max > 1) stop("Provide a value > t_min and <= 1 for `t_max`.", call. = FALSE)
  if (!edge_correction %in% c("none", "toroidal")) stop("Provide `edge_correction` as 'none' or 'toroidal'.", call. = FALSE)
  if (!is.logical(include_comp_inds)) stop("`include_comp_inds` must be TRUE/FALSE.", call. = FALSE)
  if (!is.logical(thinning)) stop("`thinning` must be TRUE/FALSE.", call. = FALSE)
  if (!is.logical(save_sims)) stop("`save_sims` must be TRUE/FALSE.", call. = FALSE)
  if (!is.logical(verbose)) stop("`verbose` must be TRUE/FALSE.", call. = FALSE)
  if (isTRUE(include_comp_inds) && (is.null(competition_radius) || !is.numeric(competition_radius) || competition_radius <= 0)) {
    stop("Provide a positive numeric `competition_radius` when `include_comp_inds = TRUE`.", call. = FALSE)
  }
  if (!is.numeric(n_sim) || n_sim < 1) stop("Provide a positive integer for `n_sim`.", call. = FALSE)
  if (is.na(seed) || seed < 0 || seed != as.integer(seed)) stop("Provide a non-negative integer for `seed`.", call. = FALSE)
  if (!is.null(raster_list) && !is.list(raster_list)) {
    stop("If provided, `raster_list` must be a list of rasters.", call. = FALSE)
  }
  if (!is.logical(scaled_rasters)) stop("`scaled_rasters` must be TRUE/FALSE.", call. = FALSE)
  if (!is.logical(parallel)) stop("`parallel` must be TRUE/FALSE.", call. = FALSE)

  .vstate <- new.env(parent = emptyenv())
  .vstate$header_printed <- FALSE
  .vmsg <- function(..., .indent = 0L) {
    if (!isTRUE(verbose)) return(invisible(NULL))
    if (!isTRUE(.vstate$header_printed)) {
      message("[ldmppr::check_model_fit]")
      .vstate$header_printed <- TRUE
    }
    indent <- if (.indent > 0L) paste(rep("  ", .indent), collapse = "") else ""
    message(indent, paste0(..., collapse = ""))
    invisible(NULL)
  }
  .tic <- function() proc.time()[["elapsed"]]
  .toc <- function(t0) proc.time()[["elapsed"]] - t0
  .fmt_sec <- function(sec) {
    sec <- as.numeric(sec)
    if (!is.finite(sec)) return("NA")
    if (sec < 60) return(sprintf("%.1fs", sec))
    if (sec < 3600) return(sprintf("%.1fmin", sec / 60))
    sprintf("%.2fh", sec / 3600)
  }
  .step_header <- function(i, n, label) .vmsg(sprintf("Step %d/%d: %s", i, n, label))

  # mark_model required only when mark_mode == "mark_model"
  if (identical(mark_mode, "mark_model")) {
    if (is.null(mark_model)) stop("Provide a `mark_model` when mark_mode='mark_model'.", call. = FALSE)
    mark_model <- as_mark_model(mark_model)
  }

  # ---- resolve sc params ----
  sc_params <- resolve_sc_params(process_fit)

  # ---- resolve xy_bounds if possible ----
  if (is.null(xy_bounds)) {
    if (!is.null(reference_data) && spatstat.geom::is.ppp(reference_data)) {
      xy_bounds <- infer_xy_bounds_from_ppp(reference_data)
    } else if (inherits(process_fit, "ldmppr_fit") &&
               !is.null(process_fit$grid$upper_bounds) &&
               length(process_fit$grid$upper_bounds) >= 3) {
      ub <- process_fit$grid$upper_bounds
      xy_bounds <- c(0, ub[2], 0, ub[3])
    }
  }
  if (is.null(xy_bounds) || length(xy_bounds) != 4) {
    stop("Provide `xy_bounds = c(a_x, b_x, a_y, b_y)`, or supply `reference_data` as a ppp so bounds can be inferred.", call. = FALSE)
  }
  if (xy_bounds[2] < xy_bounds[1] || xy_bounds[4] < xy_bounds[3]) {
    stop("Provide (x,y) bounds in the form (a_x, b_x, a_y, b_y) for `xy_bounds`.", call. = FALSE)
  }

  # ---- resolve reference ppp ----
  reference_data <- resolve_reference_ppp(reference_data, process_fit, xy_bounds)

  # ---- resolve anchor point ----
  if (is.null(anchor_point)) {
    if (inherits(process_fit, "ldmppr_fit") && !is.null(process_fit$data_original) && is.data.frame(process_fit$data_original)) {
      anchor_point <- infer_anchor_from_df(process_fit$data_original)
    }
    if (is.null(anchor_point)) anchor_point <- infer_anchor_from_ppp(reference_data)
  }
  anchor_point <- as.numeric(anchor_point)
  if (length(anchor_point) != 2) stop("Provide `anchor_point` as a length-2 numeric vector (x,y), or let it be inferred.", call. = FALSE)

  # ---- rasters only needed for mark_model mode ----
  inferred_rasters <- FALSE
  if (identical(mark_mode, "mark_model")) {
    if (is.null(raster_list)) {
      raster_list <- infer_rasters_from_mark_model(mark_model)
      inferred_rasters <- !is.null(raster_list)
      if (inferred_rasters) {
        mm_scaled <- infer_scaled_flag_from_mark_model(mark_model)
        if (!is.null(mm_scaled)) scaled_rasters <- mm_scaled
      }
    }
    if (is.null(raster_list) || !is.list(raster_list)) {
      stop("Provide `raster_list` or pass a `mark_model` that contains stored rasters.", call. = FALSE)
    }
    if (!isTRUE(scaled_rasters)) {
      raster_list <- scale_rasters(raster_list)
      scaled_rasters <- TRUE
    }
  }

  # ---- seed ----
  set.seed(seed)

  # ---- future plan handling ----
  will_parallelize <- isTRUE(parallel) && n_sim > 1L
  t_all <- .tic()
  .vmsg("Checking model fit")
  .vmsg("mark_mode=", mark_mode, ", n_sim=", n_sim, ", parallel=", will_parallelize, .indent = 1L)
  .vmsg("include_comp_inds=", include_comp_inds, ", edge_correction=", edge_correction, .indent = 1L)

  .step_header(1, 4, "Preparing simulation setup")
  t_step <- .tic()

  if (isTRUE(set_future_plan) && isTRUE(will_parallelize)) {
    if (!requireNamespace("future", quietly = TRUE)) {
      stop("Package 'future' is required when set_future_plan=TRUE and parallel=TRUE.", call. = FALSE)
    }
    original_plan <- future::plan()
    on.exit(future::plan(original_plan), add = TRUE)

    num_cores <- max(1L, as.integer(num_cores))

    use_multicore <- FALSE
    if (requireNamespace("parallelly", quietly = TRUE)) {
      use_multicore <- parallelly::supportsMulticore()
    } else {
      use_multicore <- (.Platform$OS.type == "unix")
    }
    if (isTRUE(use_multicore)) {
      future::plan(future::multicore, workers = num_cores)
    } else {
      future::plan(future::multisession, workers = num_cores)
    }

    # do NOT change user handlers; only silence stdout if desired
    old_opts <- options(future.stdout = FALSE)
    on.exit(options(old_opts), add = TRUE)

    .vmsg("Parallel future plan set with workers=", num_cores)
  }

  # ---- terra serialization safeguard (only relevant for mark_model + multisession) ----
  wrapped_rasters <- FALSE
  if (identical(mark_mode, "mark_model") &&
      isTRUE(will_parallelize) &&
      requireNamespace("terra", quietly = TRUE) &&
      length(raster_list) > 0 &&
      any(vapply(raster_list, inherits, logical(1), what = "SpatRaster"))) {

    raster_list <- lapply(raster_list, function(r) {
      if (inherits(r, "SpatRaster")) terra::wrap(r) else r
    })
    wrapped_rasters <- TRUE
  }

  # ---- choose r grid using reference K ----
  # (we keep L/K on the full K grid, like your original code)
  K_ref <- spatstat.explore::Kest(spatstat.geom::unmark(reference_data))
  d_L <- K_ref$r
  dL_len <- length(d_L)

  # ---- determine FGJ truncation from reference pattern ----
  ref_un <- spatstat.geom::unmark(reference_data)
  F_ref <- spatstat.explore::Fest(ref_un, correction = fg_correction, r = d_L)[[fg_correction]]
  G_ref <- spatstat.explore::Gest(ref_un, correction = fg_correction, r = d_L)[[fg_correction]]

  ok_FG <- is.finite(F_ref) & is.finite(G_ref) & (F_ref < 1) & (G_ref < 1)
  if (!any(ok_FG)) {
    # fallback: use the first finite r only
    first_ok <- which(is.finite(F_ref) & is.finite(G_ref))
    if (!length(first_ok)) stop("Reference F/G are non-finite for all r; cannot determine FGJ range.", call. = FALSE)
    FGJ_max_idx <- first_ok[1]
  } else {
    FGJ_max_idx <- max(which(ok_FG))
  }
  FGJ_max_idx <- max(1L, FGJ_max_idx)

  d_FGJ <- d_L[1:FGJ_max_idx]
  dFGJ_len <- length(d_FGJ)

  .vmsg("Using FGJ r-grid from reference: 1:", FGJ_max_idx,
        " (max r=", signif(max(d_FGJ), 4), "), correction=", fg_correction)
  .vmsg("Done in ", .fmt_sec(.toc(t_step)), ".", .indent = 1L)

  # ---- E/V r grid: keep full L grid for now (you said keep them in combined test) ----
  d_EV <- d_L
  dEV_len <- length(d_EV)

  # ---- storage ----
  K_PP <- matrix(NA_real_, nrow = dL_len, ncol = n_sim)
  F_PP <- matrix(NA_real_, nrow = dFGJ_len, ncol = n_sim)
  G_PP <- matrix(NA_real_, nrow = dFGJ_len, ncol = n_sim)
  J_PP <- matrix(NA_real_, nrow = dFGJ_len, ncol = n_sim)
  E_PP <- matrix(NA_real_, nrow = dEV_len, ncol = n_sim)
  V_PP <- matrix(NA_real_, nrow = dEV_len, ncol = n_sim)
  n_real <- numeric(n_sim)

  # ---- time_to_size inversion settings (from reference) ----
  inv_delta <- NA_real_
  smin <- smax <- NA_real_
  if (identical(mark_mode, "time_to_size")) {
    mref <- spatstat.geom::marks(reference_data)
    if (is.null(mref)) stop("reference_data must have marks (sizes) for mark_mode='time_to_size'.", call. = FALSE)
    mref <- as.numeric(mref)
    if (anyNA(mref) || !all(is.finite(mref))) stop("reference_data marks must be finite numeric.", call. = FALSE)
    smin <- min(mref); smax <- max(mref)
    inv_delta <- 1
    if (inherits(process_fit, "ldmppr_fit") && !is.null(process_fit$mapping$delta) && is.finite(process_fit$mapping$delta)) {
      inv_delta <- as.numeric(process_fit$mapping$delta)
      if (!is.finite(inv_delta) || inv_delta <= 0) inv_delta <- 1
    }
  }

  # ---- attempt limits ----
  if (is.null(max_attempts) || !is.finite(max_attempts) || max_attempts < n_sim) {
    # generous default: allow lots of rejects without hanging forever
    max_attempts <- 10L * as.integer(n_sim)
  }
  max_attempts <- as.integer(max_attempts)

  # ---- one accepted simulation ----
  one_accepted <- function() {

    # unwrap rasters inside worker if needed
    rl <- raster_list
    if (identical(mark_mode, "mark_model") &&
        isTRUE(wrapped_rasters) && requireNamespace("terra", quietly = TRUE)) {
      rl <- lapply(rl, function(r) if (inherits(r, "PackedSpatRaster")) terra::unwrap(r) else r)
    }

    sim_df <- if (isTRUE(thinning)) {
      simulate_sc(t_min = t_min, t_max = t_max, sc_params = sc_params,
                  anchor_point = anchor_point, xy_bounds = xy_bounds)$thinned
    } else {
      simulate_sc(t_min = t_min, t_max = t_max, sc_params = sc_params,
                  anchor_point = anchor_point, xy_bounds = xy_bounds)$unthinned
    }

    n_real_i <- nrow(sim_df)
    if (n_real_i < 1) return(NULL)

    # ---- marks ----
    if (identical(mark_mode, "mark_model")) {
      marks_i <- if (inherits(mark_model, "ldmppr_mark_model")) {
        stats::predict(
          mark_model,
          sim_realization = sim_df,
          raster_list = rl,
          scaled_rasters = TRUE,
          xy_bounds = xy_bounds,
          include_comp_inds = include_comp_inds,
          competition_radius = competition_radius,
          edge_correction = edge_correction
        )
      } else {
        predict_marks(
          sim_realization = sim_df,
          raster_list = rl,
          scaled_rasters = TRUE,
          mark_model = mark_model,
          xy_bounds = xy_bounds,
          include_comp_inds = include_comp_inds,
          competition_radius = competition_radius,
          edge_correction = edge_correction
        )
      }
    } else {
      # invert t_scaled = 1 - ((s - smin)/(smax - smin))^delta
      # -> s = smin + (smax - smin) * (1 - t)^(1/delta)
      if (!("time" %in% names(sim_df))) {
        stop("simulate_sc() must return a data.frame with a 'time' column for mark_mode='time_to_size'.", call. = FALSE)
      }
      tt <- sim_df$time
      tt <- pmin(pmax(as.numeric(tt), 0), 1)
      marks_i <- smin + (smax - smin) * (1 - tt)^(1 / inv_delta)
    }

    PP_xy <- generate_mpp(locations = sim_df, marks = marks_i, xy_bounds = xy_bounds)
    PP_un <- spatstat.geom::unmark(PP_xy)

    # ---- compute summaries ----
    K_i <- spatstat.explore::Kest(PP_un, correction = "isotropic", r = d_L)$iso
    F_i <- spatstat.explore::Fest(PP_un, correction = fg_correction, r = d_FGJ)[[fg_correction]]
    G_i <- spatstat.explore::Gest(PP_un, correction = fg_correction, r = d_FGJ)[[fg_correction]]
    J_i <- spatstat.explore::Jest(PP_un, correction = fg_correction, r = d_FGJ)[[fg_correction]] - 1
    E_i <- spatstat.explore::Emark(PP_xy, correction = "isotropic", r = d_EV)$iso
    V_i <- spatstat.explore::Vmark(PP_xy, correction = "isotropic", r = d_EV)$iso

    # ---- acceptance rule: must be fully finite on required ranges ----
    ok <- all(is.finite(K_i)) &&
      all(is.finite(F_i)) &&
      all(is.finite(G_i)) &&
      all(is.finite(J_i)) &&
      all(is.finite(E_i)) &&
      all(is.finite(V_i))

    if (!isTRUE(ok)) return(NULL)

    list(n = n_real_i, K = K_i, F = F_i, G = G_i, J = J_i, E = E_i, V = V_i)
  }

  # ---- simulate ----
  .step_header(2, 4, "Generating accepted simulations")
  t_step <- .tic()
  if (!isTRUE(will_parallelize)) {

    if (isTRUE(verbose)) {
      pb <- progress::progress_bar$new(
        format = "ldmppr::check_model_fit sims (accepted): [:bar] :percent in :elapsed, ETA: :eta",
        total = n_sim, clear = FALSE, width = 80
      )
      pb$tick(0)
    }

    accepted <- 0L
    attempts <- 0L
    while (accepted < n_sim && attempts < max_attempts) {
      attempts <- attempts + 1L
      res <- suppressPackageStartupMessages(suppressMessages(one_accepted()))
      if (!is.null(res)) {
        accepted <- accepted + 1L
        n_real[accepted] <- res$n
        K_PP[, accepted] <- res$K
        F_PP[, accepted] <- res$F
        G_PP[, accepted] <- res$G
        J_PP[, accepted] <- res$J
        E_PP[, accepted] <- res$E
        V_PP[, accepted] <- res$V
        if (isTRUE(verbose)) pb$tick()
      }
    }

    if (accepted < n_sim) {
      warning("Only accepted ", accepted, "/", n_sim,
              " simulations before hitting max_attempts=", max_attempts,
              ". Consider reducing r-range or increasing max_attempts.", call. = FALSE)
      # trim matrices to accepted
      K_PP <- K_PP[, 1:accepted, drop = FALSE]
      F_PP <- F_PP[, 1:accepted, drop = FALSE]
      G_PP <- G_PP[, 1:accepted, drop = FALSE]
      J_PP <- J_PP[, 1:accepted, drop = FALSE]
      E_PP <- E_PP[, 1:accepted, drop = FALSE]
      V_PP <- V_PP[, 1:accepted, drop = FALSE]
      n_real <- n_real[1:accepted]
      n_sim <- accepted
    }

  } else {

    # parallel: each task returns ONE accepted sim by looping locally
    if (!requireNamespace("future", quietly = TRUE)) stop("Need 'future' for parallel=TRUE.", call. = FALSE)
    if (!requireNamespace("furrr", quietly = TRUE)) stop("Need 'furrr' for parallel=TRUE.", call. = FALSE)

    use_progressr <- isTRUE(verbose) && requireNamespace("progressr", quietly = TRUE)

    get_one <- function(dummy) {
      attempts <- 0L
      repeat {
        attempts <- attempts + 1L
        res <- suppressPackageStartupMessages(suppressMessages(one_accepted()))
        if (!is.null(res)) return(res)
        if (attempts >= 5000L) return(NULL) # per-worker safety
      }
    }

    if (isTRUE(use_progressr)) {
      results <- progressr::with_progress({
        p <- progressr::progressor(steps = n_sim)
        furrr::future_map(
          seq_len(n_sim),
          function(i) { out <- get_one(i); p(); out },
          .options = furrr::furrr_options(seed = TRUE)
        )
      })
    } else {
      results <- furrr::future_map(
        seq_len(n_sim),
        get_one,
        .options = furrr::furrr_options(seed = TRUE)
      )
    }

    # fill results (drop NULLs if any)
    keep <- vapply(results, function(x) !is.null(x), logical(1))
    results <- results[keep]
    if (length(results) < n_sim) {
      warning("Parallel run returned only ", length(results), "/", n_sim,
              " accepted simulations (some workers hit safety limit).", call. = FALSE)
      n_sim <- length(results)
      K_PP <- K_PP[, 1:n_sim, drop = FALSE]
      F_PP <- F_PP[, 1:n_sim, drop = FALSE]
      G_PP <- G_PP[, 1:n_sim, drop = FALSE]
      J_PP <- J_PP[, 1:n_sim, drop = FALSE]
      E_PP <- E_PP[, 1:n_sim, drop = FALSE]
      V_PP <- V_PP[, 1:n_sim, drop = FALSE]
      n_real <- n_real[1:n_sim]
    }

    for (j in seq_len(n_sim)) {
      res <- results[[j]]
      n_real[j] <- res$n
      K_PP[, j] <- res$K
      F_PP[, j] <- res$F
      G_PP[, j] <- res$G
      J_PP[, j] <- res$J
      E_PP[, j] <- res$E
      V_PP[, j] <- res$V
    }
  }
  .vmsg("Done in ", .fmt_sec(.toc(t_step)), ".", .indent = 1L)

  sim_metrics <- if (isTRUE(save_sims)) {
    list(Ksim = K_PP, Fsim = F_PP, Gsim = G_PP, Jsim = J_PP, Esim = E_PP, Vsim = V_PP, n_per = n_real)
  } else {
    NULL
  }

  # ---- calculate envelopes ----
  .step_header(3, 4, "Computing envelope tests")
  t_step <- .tic()

  # L (from K) envelopes
  C_ref_L <- GET::create_curve_set(list(
    r = d_L,
    obs = sqrt(K_ref$iso / pi) - d_L,
    theo = sqrt(K_ref$theo / pi) - d_L,
    sim_m = sqrt(K_PP / pi) - d_L
  ))
  r_envL <- GET::global_envelope_test(C_ref_L, type = "rank")

  # F envelopes on truncated reference-driven r grid
  F_ref_use <- spatstat.explore::Fest(ref_un, correction = fg_correction, r = d_FGJ)[[fg_correction]]
  F_theo_use <- spatstat.explore::Fest(ref_un, correction = fg_correction, r = d_FGJ)$theo
  C_ref_F <- GET::create_curve_set(list(
    r = d_FGJ,
    obs = F_ref_use,
    theo = F_theo_use,
    sim_m = F_PP
  ))
  r_envF <- GET::global_envelope_test(C_ref_F, type = "rank")

  # G envelopes
  G_ref_use <- spatstat.explore::Gest(ref_un, correction = fg_correction, r = d_FGJ)[[fg_correction]]
  G_theo_use <- spatstat.explore::Gest(ref_un, correction = fg_correction, r = d_FGJ)$theo
  C_ref_G <- GET::create_curve_set(list(
    r = d_FGJ,
    obs = G_ref_use,
    theo = G_theo_use,
    sim_m = G_PP
  ))
  r_envG <- GET::global_envelope_test(C_ref_G, type = "rank")

  # J envelopes (scale like your earlier approach; use sim max)
  J_ref_use <- spatstat.explore::Jest(ref_un, correction = fg_correction, r = d_FGJ)[[fg_correction]] - 1
  J_theo_use <- spatstat.explore::Jest(ref_un, correction = fg_correction, r = d_FGJ)$theo - 1

  Jscale <- max(J_PP, na.rm = TRUE)
  if (!is.finite(Jscale) || Jscale <= 0) Jscale <- 1

  C_ref_J <- GET::create_curve_set(list(
    r = d_FGJ,
    obs  = J_ref_use / Jscale,
    theo = J_theo_use / Jscale,
    sim_m = J_PP / Jscale
  ))
  r_envJ <- GET::global_envelope_test(C_ref_J, type = "rank")

  # E envelopes (full grid for now)
  E_ref <- spatstat.explore::Emark(reference_data, correction = "isotropic", r = d_EV)$iso
  E_theo <- spatstat.explore::Emark(reference_data, correction = "isotropic", r = d_EV)$theo
  C_ref_E <- GET::create_curve_set(list(
    r = d_EV,
    obs = E_ref,
    theo = E_theo,
    sim_m = E_PP
  ))
  r_envE <- GET::global_envelope_test(C_ref_E, type = "rank")

  # V envelopes
  V_ref <- spatstat.explore::Vmark(reference_data, correction = "isotropic", r = d_EV)$iso
  V_theo <- spatstat.explore::Vmark(reference_data, correction = "isotropic", r = d_EV)$theo
  C_ref_V <- GET::create_curve_set(list(
    r = d_EV,
    obs = V_ref,
    theo = V_theo,
    sim_m = V_PP
  ))
  r_envV <- GET::global_envelope_test(C_ref_V, type = "rank")

  # Combined global envelope test
  r_envComb <- GET::global_envelope_test(
    curve_sets = list(L = C_ref_L, F = C_ref_F, G = C_ref_G, J = C_ref_J, E = C_ref_E, V = C_ref_V),
    type = "rank"
  )
  .vmsg("Done in ", .fmt_sec(.toc(t_step)), ".", .indent = 1L)

  envs <- list(L = r_envL, F = r_envF, G = r_envG, J = r_envJ, E = r_envE, V = r_envV)
  curve_sets <- list(L = C_ref_L, F = C_ref_F, G = C_ref_G, J = C_ref_J, E = C_ref_E, V = C_ref_V)

  settings <- list(
    t_min = t_min,
    t_max = t_max,
    thinning = thinning,
    edge_correction = edge_correction,
    include_comp_inds = include_comp_inds,
    competition_radius = competition_radius,
    n_sim = n_sim,
    seed = seed,
    mark_mode = mark_mode,
    fg_correction = fg_correction,
    parallel = parallel,
    num_cores = as.integer(num_cores),
    set_future_plan = set_future_plan,
    inferred = list(
      reference_data = is.null(match.call()$reference_data),
      xy_bounds = is.null(match.call()$xy_bounds),
      anchor_point = is.null(match.call()$anchor_point),
      raster_list = is.null(match.call()$raster_list)
    ),
    r_grids = list(
      L = d_L,
      FGJ = d_FGJ,
      EV = d_EV
    )
  )

  .step_header(4, 4, "Finalizing output object")
  .vmsg("Done in ", .fmt_sec(.toc(t_all)), ".", .indent = 1L)
  .vmsg("Model check complete.")

  new_ldmppr_model_check(
    combined_env = r_envComb,
    envs = envs,
    curve_sets = curve_sets,
    sim_metrics = sim_metrics,
    settings = settings,
    call = match.call()
  )
}
