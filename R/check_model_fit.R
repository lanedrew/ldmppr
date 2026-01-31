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
#' @param raster_list a list of raster objects used for predicting marks.
#' @param scaled_rasters \code{TRUE} or \code{FALSE} indicating whether the rasters have already been scaled.
#' @param mark_model a mark model object. May be a \code{ldmppr_mark_model} or a legacy model.
#' @param xy_bounds (optional) vector of bounds as \code{c(a_x, b_x, a_y, b_y)}. If \code{NULL}, will be
#'   inferred from \code{reference_data}'s window when \code{reference_data} is provided,
#'   otherwise from \code{ldmppr_fit} with lower bounds assumed to be 0.
#' @param include_comp_inds \code{TRUE} or \code{FALSE} indicating whether to compute competition indices.
#' @param competition_radius distance for competition radius if \code{include_comp_inds = TRUE}.
#' @param thinning \code{TRUE} or \code{FALSE} indicating whether to use the thinned simulated values.
#' @param edge_correction type of edge correction to apply (\code{"none"} or \code{"toroidal"}).
#' @param n_sim number of simulated datasets to generate.
#' @param save_sims \code{TRUE} or \code{FALSE} indicating whether to save and return the simulated metrics.
#' @param verbose \code{TRUE} or \code{FALSE} indicating whether to show progress of model checking.
#'   When \code{TRUE}, progress is reported via \pkg{progressr} (if available) and is compatible with parallel execution.
#' @param seed integer seed for reproducibility.
#' @param parallel \code{TRUE} or \code{FALSE}. If \code{TRUE}, simulations are run in parallel via \pkg{furrr}/\pkg{future}.
#' @param num_cores number of workers to use when \code{parallel=TRUE}. Defaults to one fewer than the number of detected cores.
#' @param set_future_plan \code{TRUE} or \code{FALSE}. If \code{TRUE} and \code{parallel=TRUE}, set a local \pkg{future}
#'   plan internally (default behavior uses \code{multisession}).
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
#' # NOTE: examples are run by CRAN with --run-donttest; keep parallel=FALSE here.
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
                            set_future_plan = FALSE) {

  process <- match.arg(process)

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
  if (is.null(raster_list) || !is.list(raster_list)) stop("Provide a list of rasters for `raster_list`.", call. = FALSE)
  if (!is.logical(scaled_rasters)) stop("`scaled_rasters` must be TRUE/FALSE.", call. = FALSE)
  if (is.null(mark_model)) stop("Provide a `mark_model` (ldmppr_mark_model or legacy) for predicting marks.", call. = FALSE)
  if (!is.logical(parallel)) stop("`parallel` must be TRUE/FALSE.", call. = FALSE)

  # ---- resolve mark model ----
  mark_model <- as_mark_model(mark_model)

  # ---- seed ----
  set.seed(seed)

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

  # ---- scale rasters if necessary ----
  if (!isTRUE(scaled_rasters)) {
    raster_list <- scale_rasters(raster_list)
  }

  # ---- parallel plan handling (multisession default) ----
  will_parallelize <- isTRUE(parallel) && n_sim > 1L
  if (isTRUE(set_future_plan) && isTRUE(will_parallelize)) {
    if (!requireNamespace("future", quietly = TRUE)) {
      stop("Package 'future' is required when set_future_plan=TRUE and parallel=TRUE.", call. = FALSE)
    }
    original_plan <- future::plan()
    on.exit(future::plan(original_plan), add = TRUE)

    num_cores <- max(1L, as.integer(num_cores))
    future::plan(future::multisession, workers = num_cores)

    if (isTRUE(verbose)) message("future plan set to multisession with workers = ", num_cores)
  }

  # ---- terra serialization safeguard for multisession ----
  # Wrap SpatRaster so it can be serialized to multisession workers (prevents "external pointer is not valid")
  wrapped_rasters <- FALSE
  if (isTRUE(will_parallelize) &&
      requireNamespace("terra", quietly = TRUE) &&
      length(raster_list) > 0 &&
      any(vapply(raster_list, inherits, logical(1), what = "SpatRaster"))) {

    raster_list <- lapply(raster_list, function(r) {
      if (inherits(r, "SpatRaster")) terra::wrap(r) else r
    })
    wrapped_rasters <- TRUE
  }

  # ---- choose r grid using reference K ----
  K_ref <- spatstat.explore::Kest(spatstat.geom::unmark(reference_data))
  d <- K_ref$r
  d_length <- length(d)

  # ---- storage ----
  K_PP <- matrix(NA_real_, nrow = d_length, ncol = n_sim)
  F_PP <- matrix(NA_real_, nrow = d_length, ncol = n_sim)
  G_PP <- matrix(NA_real_, nrow = d_length, ncol = n_sim)
  J_PP <- matrix(NA_real_, nrow = d_length, ncol = n_sim)
  E_PP <- matrix(NA_real_, nrow = d_length, ncol = n_sim)
  V_PP <- matrix(NA_real_, nrow = d_length, ncol = n_sim)
  n_real <- numeric(n_sim)

  sim_one <- function() {
    # Unwrap rasters inside worker if needed (namespace-only loads; no attach/startup chatter)
    rl <- raster_list
    if (isTRUE(wrapped_rasters) && requireNamespace("terra", quietly = TRUE)) {
      rl <- lapply(rl, function(r) {
        if (inherits(r, "PackedSpatRaster")) terra::unwrap(r) else r
      })
    }

    sim_df <- if (isTRUE(thinning)) {
      simulate_sc(t_min = t_min, t_max = t_max, sc_params = sc_params,
                  anchor_point = anchor_point, xy_bounds = xy_bounds)$thinned
    } else {
      simulate_sc(t_min = t_min, t_max = t_max, sc_params = sc_params,
                  anchor_point = anchor_point, xy_bounds = xy_bounds)$unthinned
    }

    n_real_i <- nrow(sim_df)
    if (n_real_i < 1) {
      return(list(n = 0,
                  K = rep(NA_real_, d_length), F = rep(NA_real_, d_length), G = rep(NA_real_, d_length),
                  J = rep(NA_real_, d_length), E = rep(NA_real_, d_length), V = rep(NA_real_, d_length)))
    }

    pred_marks <- predict_marks(
      sim_realization = sim_df,
      raster_list = rl,
      scaled_rasters = TRUE,
      mark_model = mark_model,
      xy_bounds = xy_bounds,
      include_comp_inds = include_comp_inds,
      competition_radius = competition_radius,
      edge_correction = edge_correction
    )

    PP_xy <- generate_mpp(locations = sim_df, marks = pred_marks, xy_bounds = xy_bounds)

    list(
      n = n_real_i,
      K = spatstat.explore::Kest(spatstat.geom::unmark(PP_xy), correction = "isotropic", r = d)$iso,
      F = spatstat.explore::Fest(spatstat.geom::unmark(PP_xy), correction = "rs", r = d)$rs,
      G = spatstat.explore::Gest(spatstat.geom::unmark(PP_xy), correction = "rs", r = d)$rs,
      J = spatstat.explore::Jest(spatstat.geom::unmark(PP_xy), correction = "rs", r = d)$rs - 1,
      E = spatstat.explore::Emark(PP_xy, correction = "isotropic", r = d)$iso,
      V = spatstat.explore::Vmark(PP_xy, correction = "isotropic", r = d)$iso
    )
  }

  # ---- simulate (parallel or sequential) ----
  if (isTRUE(will_parallelize)) {
    if (!requireNamespace("furrr", quietly = TRUE) || !requireNamespace("future", quietly = TRUE)) {
      stop("Parallel execution requires packages 'future' and 'furrr'.", call. = FALSE)
    }

    # Goal 1: avoid repeated "Registered S3 method overwritten ..." spam from worker startup
    #   -> do NOT attach packages on workers (packages = character(0))
    #   -> do NOT relay packageStartupMessage conditions back to the main R session
    #
    # Goal 2: progress that doesn't look "stuck" in RStudio
    #   -> use smaller future chunks so progress updates arrive sooner
    #   -> use the RStudio handler when available
    workers_for_chunks <- max(1L, as.integer(num_cores))
    chunk_size <- max(1L, floor(n_sim / (workers_for_chunks * 20L)))  # ~20 chunks/worker

    furrr_opts <- furrr::furrr_options(
      seed = TRUE,
      packages = character(0),
      stdout = FALSE,
      # Relay all conditions EXCEPT package startup messages (this is what prints the overwrite notices)
      conditions = structure("condition", exclude = c("packageStartupMessage")),
      # Smaller chunks => earlier/more frequent progress updates
      chunk_size = chunk_size
    )

    use_progressr <- isTRUE(verbose) && requireNamespace("progressr", quietly = TRUE)

    if (isTRUE(use_progressr)) {
      # Choose a handler that works well in RStudio on macOS
      handler_name <- "txtprogressbar"
      if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
        handler_name <- "rstudio"
      }

      old_handlers <- try(progressr::handlers(), silent = TRUE)
      on.exit({
        if (!inherits(old_handlers, "try-error")) progressr::handlers(old_handlers)
      }, add = TRUE)

      progressr::handlers(handler_name)

      results <- progressr::with_progress({
        p <- progressr::progressor(steps = n_sim)
        furrr::future_map(
          seq_len(n_sim),
          function(i) {
            out <- sim_one()
            p()
            out
          },
          .options = furrr_opts
        )
      })
    } else {
      results <- furrr::future_map(
        seq_len(n_sim),
        function(i) sim_one(),
        .options = furrr_opts
      )
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

  } else {
    # Sequential: keep your existing progress behavior
    if (isTRUE(verbose)) {
      pb <- progress::progress_bar$new(
        format = "Model check simulations: [:bar] :percent in :elapsed, ETA: :eta",
        total = n_sim, clear = FALSE, width = 80
      )
      pb$tick(0)
      for (j in seq_len(n_sim)) {
        res <- sim_one()
        n_real[j] <- res$n
        K_PP[, j] <- res$K
        F_PP[, j] <- res$F
        G_PP[, j] <- res$G
        J_PP[, j] <- res$J
        E_PP[, j] <- res$E
        V_PP[, j] <- res$V
        pb$tick()
      }
    } else {
      for (j in seq_len(n_sim)) {
        res <- sim_one()
        n_real[j] <- res$n
        K_PP[, j] <- res$K
        F_PP[, j] <- res$F
        G_PP[, j] <- res$G
        J_PP[, j] <- res$J
        E_PP[, j] <- res$E
        V_PP[, j] <- res$V
      }
    }
  }

  sim_metrics <- if (isTRUE(save_sims)) {
    list(Ksim = K_PP, Fsim = F_PP, Gsim = G_PP, Jsim = J_PP, Esim = E_PP, Vsim = V_PP, n_per = n_real)
  } else {
    NULL
  }

  # ---- envelopes ----
  C_ref_L <- GET::create_curve_set(list(
    r = d,
    obs = sqrt(K_ref$iso / pi) - d,
    theo = sqrt(K_ref$theo / pi) - d,
    sim_m = sqrt(K_PP / pi) - d
  ))
  r_envL <- GET::global_envelope_test(C_ref_L, type = "rank")

  safe_crossing_idx <- function(mat) {
    idx <- apply(mat, 2, function(x) {
      w <- which(x >= 1)
      if (length(w) == 0) return(NA_integer_)
      min(w)
    })
    idx <- idx[!is.na(idx)]
    if (!length(idx)) return(1L)
    max(idx)
  }

  F_val <- safe_crossing_idx(F_PP)
  C_ref_F <- GET::create_curve_set(list(
    r = d[1:F_val],
    obs = spatstat.explore::Fest(spatstat.geom::unmark(reference_data), r = d[1:F_val])$rs,
    theo = spatstat.explore::Fest(spatstat.geom::unmark(reference_data), r = d[1:F_val])$theo,
    sim_m = F_PP[1:F_val, , drop = FALSE]
  ))
  r_envF <- GET::global_envelope_test(C_ref_F, type = "rank")

  G_val <- safe_crossing_idx(G_PP)
  C_ref_G <- GET::create_curve_set(list(
    r = d[1:G_val],
    obs = spatstat.explore::Gest(spatstat.geom::unmark(reference_data), r = d[1:G_val])$rs,
    theo = spatstat.explore::Gest(spatstat.geom::unmark(reference_data), r = d[1:G_val])$theo,
    sim_m = G_PP[1:G_val, , drop = FALSE]
  ))
  r_envG <- GET::global_envelope_test(C_ref_G, type = "rank")

  J_ref <- spatstat.explore::Jest(spatstat.geom::unmark(reference_data), r = d)$rs
  J_val <- min(c(
    min(apply(F_PP, 2, function(x) sum(x < 1, na.rm = TRUE))),
    min(apply(G_PP, 2, function(x) sum(x < 1, na.rm = TRUE))),
    sum(!is.na(J_ref))
  ))
  J_val <- max(1L, J_val)
  J_scale <- apply(J_PP[1:J_val, , drop = FALSE], 1, max, na.rm = TRUE)
  J_scale[J_scale == 0 | is.na(J_scale)] <- 1

  if (any(is.infinite(J_PP[1:J_val, ]) | is.na(J_PP[1:J_val, ]))) {
    warning("J_PP contains Inf or NA values in the range used for envelopes.")
  }

  C_ref_J <- GET::create_curve_set(list(
    r = d[1:J_val],
    obs = (J_ref[1:J_val] - 1) / J_scale,
    theo = spatstat.explore::Jest(spatstat.geom::unmark(reference_data), r = d[1:J_val])$theo - 1,
    sim_m = J_PP[1:J_val, , drop = FALSE] / J_scale
  ))
  r_envJ <- GET::global_envelope_test(C_ref_J, type = "rank")

  C_ref_E <- GET::create_curve_set(list(
    r = d,
    obs = spatstat.explore::Emark(reference_data, correction = "isotropic", r = d)$iso,
    theo = spatstat.explore::Emark(reference_data, correction = "isotropic", r = d)$theo,
    sim_m = E_PP
  ))
  r_envE <- GET::global_envelope_test(C_ref_E, type = "rank")

  C_ref_V <- GET::create_curve_set(list(
    r = d,
    obs = spatstat.explore::Vmark(reference_data, correction = "isotropic", r = d)$iso,
    theo = spatstat.explore::Vmark(reference_data, correction = "isotropic", r = d)$theo,
    sim_m = V_PP
  ))
  r_envV <- GET::global_envelope_test(C_ref_V, type = "rank")

  r_envComb <- GET::global_envelope_test(
    curve_sets = list(L = C_ref_L, F = C_ref_F, G = C_ref_G, J = C_ref_J, E = C_ref_E, V = C_ref_V),
    type = "rank"
  )

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
    parallel = parallel,
    num_cores = as.integer(num_cores),
    set_future_plan = set_future_plan,
    inferred = list(
      reference_data = is.null(match.call()$reference_data),
      xy_bounds = is.null(match.call()$xy_bounds),
      anchor_point = is.null(match.call()$anchor_point)
    )
  )

  new_ldmppr_model_check(
    combined_env = r_envComb,
    envs = envs,
    curve_sets = curve_sets,
    sim_metrics = sim_metrics,
    settings = settings,
    call = match.call()
  )
}
