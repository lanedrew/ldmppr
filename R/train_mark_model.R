#' Train a flexible model for the mark distribution
#'
#' @description
#' Trains a predictive model for the mark distribution of a spatio-temporal process.
#' \code{data} may be either (1) a data.frame containing columns \code{x}, \code{y}, \code{size} and \code{time},
#' (2) a data.frame containing \code{x}, \code{y}, \code{size} (time will be derived via \code{delta}),
#' or (3) a \code{ldmppr_fit} object returned by \code{\link{estimate_process_parameters}}.
#' Allows the user to incorporate location specific information and competition indices as covariates in the mark model.
#'
#' @param data a data.frame or a \code{ldmppr_fit} object. See Description.
#' @param raster_list list of raster objects used for mark-model training.
#' @param scaled_rasters \code{TRUE} or \code{FALSE} indicating whether rasters are already scaled.
#' @param model_type the machine learning model type (\code{"xgboost"} or \code{"random_forest"}).
#' @param xy_bounds a vector of domain bounds (2 for x, 2 for y). If \code{data} is an \code{ldmppr_fit}
#'   and \code{xy_bounds} is \code{NULL}, defaults to \code{c(0, b_x, 0, b_y)} derived from fit.
#' @param delta (optional) numeric scalar used only when \code{data} contains \code{(x,y,size)} but not \code{time}.
#'   If \code{data} is an \code{ldmppr_fit} and time is missing, the function will infer the \code{delta} value from the fit.
#' @param save_model \code{TRUE} or \code{FALSE} indicating whether to save the generated model.
#' @param save_path path for saving the generated model.
#' @param parallel \code{TRUE} or \code{FALSE}. If \code{TRUE}, tuning is parallelized over resamples.
#'   For small datasets, parallel overhead may outweigh speed gains.
#' @param num_cores number of workers to use when \code{parallel=TRUE}. Ignored when \code{parallel=FALSE}.
#' @param include_comp_inds \code{TRUE} or \code{FALSE} indicating whether to generate and use competition indices as covariates.
#' @param competition_radius positive numeric distance used when \code{include_comp_inds = TRUE}.
#' @param edge_correction type of edge correction to apply (\code{"none"}, \code{"toroidal"}, or \code{"truncation"}).
#' @param selection_metric metric to use for identifying the optimal model (\code{"rmse"}, \code{"mae"}, or \code{"rsq"}).
#' @param cv_folds number of cross-validation folds to use in model training.
#'   If \code{cv_folds <= 1}, tuning is skipped and the model is fit once with default hyperparameters.
#' @param tuning_grid_size size of the tuning grid for hyperparameter tuning.
#' @param seed integer seed for reproducible resampling/tuning/model fitting.
#' @param verbose \code{TRUE} or \code{FALSE} indicating whether to show progress of model training.
#' @return an object of class \code{"ldmppr_mark_model"} containing the trained mark model.
#'
#' @examples
#' # Load the small example data
#' data(small_example_data)
#'
#' # Load example raster data
#' raster_paths <- list.files(system.file("extdata", package = "ldmppr"),
#'   pattern = "\\.tif$", full.names = TRUE
#' )
#' raster_paths <- raster_paths[!grepl("_med\\.tif$", raster_paths)]
#' rasters <- lapply(raster_paths, terra::rast)
#'
#' # Scale the rasters
#' scaled_raster_list <- scale_rasters(rasters)
#'
#'
#' # Train the model
#' mark_model <- train_mark_model(
#'   data = small_example_data,
#'   raster_list = scaled_raster_list,
#'   scaled_rasters = TRUE,
#'   model_type = "xgboost",
#'   xy_bounds = c(0, 25, 0, 25),
#'   delta = 1,
#'   parallel = FALSE,
#'   include_comp_inds = FALSE,
#'   competition_radius = 10,
#'   edge_correction = "none",
#'   selection_metric = "rmse",
#'   cv_folds = 3,
#'   tuning_grid_size = 2,
#'   verbose = TRUE
#' )
#'
#' print(mark_model)
#'
#' @export
train_mark_model <- function(data,
                             raster_list = NULL,
                             scaled_rasters = FALSE,
                             model_type = "xgboost",
                             xy_bounds = NULL,
                             delta = NULL,
                             save_model = FALSE,
                             save_path = NULL,
                             parallel = FALSE,
                             num_cores = NULL,
                             include_comp_inds = FALSE,
                             competition_radius = 15,
                             edge_correction = "none",
                             selection_metric = "rmse",
                             cv_folds = 5,
                             tuning_grid_size = 200,
                             seed = 0,
                             verbose = TRUE) {

  # Residual-bootstrap settings are intentionally internal-only.
  store_residual_bootstrap <- TRUE
  residual_source <- "oos"
  residual_transform <- "sqrt"
  residual_bin_policy <- "quantile"
  resid_bins <- NULL
  resid_min_per_bin <- NULL
  residual_seed <- 0L

  # -------------------------
  # helpers (local)
  # -------------------------
  .elapsed_sec <- function(t0) as.numeric((proc.time() - t0)[3])
  .fmt_time <- function(x) {
    if (!is.finite(x)) return("NA")
    if (x < 60) return(sprintf("%.1fs", x))
    if (x < 3600) return(sprintf("%.1fm", x / 60))
    sprintf("%.2fh", x / 3600)
  }
  .vstate <- new.env(parent = emptyenv())
  .vstate$header_printed <- FALSE
  .vcat <- function(..., .indent = 0L) {
    if (!isTRUE(verbose)) return(invisible(NULL))
    if (!isTRUE(.vstate$header_printed)) {
      message("[ldmppr::train_mark_model]")
      .vstate$header_printed <- TRUE
    }
    indent <- if (.indent > 0L) paste(rep("  ", .indent), collapse = "") else ""
    message(indent, paste0(..., collapse = ""))
    invisible(NULL)
  }
  .step_header <- function(i, n, label) .vcat(sprintf("Step %d/%d: %s", i, n, label))

  .tf <- function(y, transform = c("sqrt", "log1p", "none")) {
    transform <- match.arg(transform)
    y0 <- pmax(as.numeric(y), 0)
    switch(transform,
           sqrt  = sqrt(y0),
           log1p = log1p(y0),
           none  = as.numeric(y))
  }

  .inv_tf <- function(z, transform = c("sqrt", "log1p", "none")) {
    transform <- match.arg(transform)
    z <- as.numeric(z)
    switch(transform,
           sqrt  = pmax(z, 0)^2,
           log1p = pmax(expm1(z), 0),
           none  = pmax(z, 0))
  }

  .pred_vec <- function(fit_or_wf, new_data) {
    p <- stats::predict(fit_or_wf, new_data = new_data)
    if (is.data.frame(p)) {
      if (".pred" %in% names(p)) return(as.numeric(p[[".pred"]]))
      if (ncol(p) == 1) return(as.numeric(p[[1]]))
      stop("Prediction returned a data.frame with unexpected columns.", call. = FALSE)
    }
    as.numeric(p)
  }

  # Build residual pools on transformed scale; bins are on transformed mu-scale
  .build_resid_pools <- function(y_raw, mu_raw,
                                 transform = c("sqrt", "log1p", "none"),
                                 bin_policy = c("quantile", "equal_width"),
                                 n_bins = 8L,
                                 min_per_bin = 8L) {
    transform <- match.arg(transform)
    bin_policy <- match.arg(bin_policy)

    y_t  <- .tf(y_raw, transform)
    mu_t <- .tf(mu_raw, transform)
    r_t  <- y_t - mu_t

    ok <- is.finite(r_t) & is.finite(mu_t)
    r_t  <- r_t[ok]
    mu_t <- mu_t[ok]

    if (!length(r_t)) {
      return(list(
        enabled = TRUE,
        source = NA_character_,
        transform = transform,
        bin_policy = bin_policy,
        breaks = c(-Inf, Inf),
        pools = list("1" = numeric(0)),
        n_bins = 1L,
        min_per_bin = min_per_bin
      ))
    }

    # breaks on mu_t
    if (bin_policy == "quantile") {
      probs <- seq(0, 1, length.out = n_bins + 1L)
      brks <- as.numeric(stats::quantile(mu_t, probs = probs, na.rm = TRUE, names = FALSE))
    } else {
      rng <- range(mu_t, na.rm = TRUE)
      if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] == rng[2]) {
        brks <- c(rng[1], rng[2])
      } else {
        brks <- seq(rng[1], rng[2], length.out = n_bins + 1L)
      }
    }

    brks <- unique(brks)
    if (length(brks) < 3L) {
      return(list(
        enabled = TRUE,
        source = NA_character_,
        transform = transform,
        bin_policy = bin_policy,
        breaks = c(-Inf, Inf),
        pools = list("1" = r_t),
        n_bins = 1L,
        min_per_bin = min_per_bin
      ))
    }

    bin_id <- cut(mu_t, breaks = brks, include.lowest = TRUE, labels = FALSE)
    pools <- split(r_t, bin_id)
    sizes <- vapply(pools, length, integer(1))

    # Merge sparse bins into nearest “good” bins (by bin index distance)
    if (any(sizes < min_per_bin)) {
      good <- which(sizes >= min_per_bin)
      if (!length(good)) {
        return(list(
          enabled = TRUE,
          source = NA_character_,
          transform = transform,
          bin_policy = bin_policy,
          breaks = c(-Inf, Inf),
          pools = list("1" = r_t),
          n_bins = 1L,
          min_per_bin = min_per_bin
        ))
      }

      map_to_good <- vapply(seq_along(pools), function(i) good[which.min(abs(good - i))], integer(1))
      new_pools <- lapply(good, function(g) unlist(pools[map_to_good == g], use.names = FALSE))

      # keep original interior breakpoints, but store as a valid monotone set; findInterval will still work
      breaks_out <- c(-Inf, brks[-c(1, length(brks))], Inf)

      return(list(
        enabled = TRUE,
        source = NA_character_,
        transform = transform,
        bin_policy = bin_policy,
        breaks = breaks_out,
        pools = stats::setNames(new_pools, as.character(seq_along(new_pools))),
        n_bins = length(new_pools),
        min_per_bin = min_per_bin
      ))
    }

    list(
      enabled = TRUE,
      source = NA_character_,
      transform = transform,
      bin_policy = bin_policy,
      breaks = c(-Inf, brks[-c(1, length(brks))], Inf),
      pools = stats::setNames(pools, as.character(seq_along(pools))),
      n_bins = length(pools),
      min_per_bin = min_per_bin
    )
  }

  t_total <- proc.time()

  # -------------------------
  # argument checks
  # -------------------------
  if (is.null(raster_list) || !is.list(raster_list)) {
    stop("Provide a list of rasters for the raster_list argument.", call. = FALSE)
  }
  if (!is.logical(scaled_rasters) || length(scaled_rasters) != 1L) {
    stop("Provide a single logical value for scaled_rasters.", call. = FALSE)
  }

  if (!model_type %in% c("xgboost", "random_forest")) {
    stop("Provide a valid model type for model_type ('xgboost' or 'random_forest').", call. = FALSE)
  }
  if (!edge_correction %in% c("none", "toroidal", "truncation")) {
    stop("Provide a valid correction type for edge_correction ('none', 'toroidal', 'truncation').", call. = FALSE)
  }
  if (!selection_metric %in% c("rmse", "mae", "rsq")) {
    stop("Provide a valid metric for selection_metric ('rmse', 'mae', 'rsq').", call. = FALSE)
  }
  if (!is.logical(parallel) || length(parallel) != 1L) stop("Provide a single logical value for parallel.", call. = FALSE)
  if (isTRUE(parallel) && !is.null(num_cores)) {
    if (!is.numeric(num_cores) || length(num_cores) != 1L || num_cores < 1) stop("Provide num_cores >= 1.", call. = FALSE)
  }
  if (!is.logical(include_comp_inds) || length(include_comp_inds) != 1L) stop("Provide a single logical value for include_comp_inds.", call. = FALSE)
  if (!is.logical(save_model) || length(save_model) != 1L) stop("Provide a single logical value for save_model.", call. = FALSE)
  if (isTRUE(save_model) && is.null(save_path)) stop("Provide save_path when save_model=TRUE.", call. = FALSE)
  if (!is.logical(verbose) || length(verbose) != 1L) stop("Provide a single logical value for verbose.", call. = FALSE)
  # -------------------------
  # high-level banner
  # -------------------------
  .vcat("Training mark model")
  .vcat("Model type: ", model_type, .indent = 1L)
  .vcat("Selection metric: ", selection_metric, .indent = 1L)
  .vcat("CV folds: ", cv_folds, .indent = 1L)
  .vcat("Tuning grid size: ", tuning_grid_size, .indent = 1L)
  .vcat("Include competition indices: ", if (isTRUE(include_comp_inds)) "yes" else "no", .indent = 1L)
  .vcat("Edge correction: ", edge_correction, .indent = 1L)

  # -------------------------
  # Step 1: coerce training data
  # -------------------------
  n_steps <- 6L
  step_t <- proc.time()
  .step_header(1, n_steps, "Preparing training data...")
  coerced <- .coerce_training_df(data, delta = delta, xy_bounds = xy_bounds)
  df <- coerced$df
  xy_bounds <- coerced$xy_bounds

  if (is.null(xy_bounds) || length(xy_bounds) != 4) {
    stop("Provide xy_bounds = c(a_x, b_x, a_y, b_y), or supply an ldmppr_fit with grid$upper_bounds.", call. = FALSE)
  }
  if (xy_bounds[2] < xy_bounds[1] || xy_bounds[4] < xy_bounds[3]) {
    stop("Provide xy_bounds in the form (a_x, b_x, a_y, b_y) with b_x >= a_x and b_y >= a_y.", call. = FALSE)
  }

  cv_folds <- as.integer(cv_folds)
  if (is.na(cv_folds) || cv_folds < 1L) stop("cv_folds must be >= 1.", call. = FALSE)

  tuning_grid_size <- as.integer(tuning_grid_size)
  if (is.na(tuning_grid_size) || tuning_grid_size < 1L) stop("tuning_grid_size must be >= 1.", call. = FALSE)

  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) || !is.finite(seed)) {
    stop("Provide a single finite numeric value for `seed`.", call. = FALSE)
  }
  seed <- as.integer(seed)
  set.seed(seed)

  .vcat("Rows: ", nrow(df), .indent = 1L)
  .vcat("Seed: ", seed, .indent = 1L)
  .vcat("Done in ", .fmt_time(.elapsed_sec(step_t)), ".", .indent = 1L)

  # -------------------------
  # Step 2: parallel backend (foreach)
  # -------------------------
  step_t <- proc.time()
  .step_header(2, n_steps, "Configuring parallel backend...")
  cl <- NULL
  n_workers <- 1L
  if (isTRUE(parallel)) {
    n_workers <- if (!is.null(num_cores)) as.integer(num_cores) else max(1L, floor(parallel::detectCores() / 2))
    cl <- parallel::makePSOCKcluster(n_workers)
    doParallel::registerDoParallel(cl)
    on.exit({
      parallel::stopCluster(cl)
      foreach::registerDoSEQ()
    }, add = TRUE)
    .vcat("Parallel: on (PSOCK workers = ", n_workers, ")", .indent = 1L)
  } else {
    foreach::registerDoSEQ()
    .vcat("Parallel: off", .indent = 1L)
  }

  engine_threads <- if (isTRUE(parallel)) {
    # avoid nested parallelism when tuning over resamples
    1L
  } else if (!is.null(num_cores)) {
    max(1L, as.integer(num_cores))
  } else {
    # conservative default for small datasets
    1L
  }
  .vcat("Model engine threads: ", engine_threads, .indent = 1L)
  .vcat("Done in ", .fmt_time(.elapsed_sec(step_t)), ".", .indent = 1L)

  # -------------------------
  # Step 3: raster scaling + covariate extraction
  # -------------------------
  step_t <- proc.time()
  .step_header(3, n_steps, "Extracting raster covariates...")
  if (!isTRUE(scaled_rasters)) {
    raster_list <- scale_rasters(raster_list)
    .vcat("Scaled rasters internally (scaled_rasters = FALSE).", .indent = 1L)
  } else {
    .vcat("Using pre-scaled rasters (scaled_rasters = TRUE).", .indent = 1L)
  }

  s <- as.matrix(df[, c("x", "y")])
  X <- extract_covars(locations = s, raster_list = raster_list)

  if ("ID" %in% names(X)) X <- X[, names(X) != "ID", drop = FALSE]
  names(X) <- make.unique(names(X), sep = "__")

  X$x <- df$x
  X$y <- df$y
  X$time <- df$time

  .vcat("Extracted ", ncol(X) - 3L, " raster feature(s).", .indent = 1L)
  .vcat("Done in ", .fmt_time(.elapsed_sec(step_t)), ".", .indent = 1L)

  # -------------------------
  # Step 4: competition indices (optional)
  # -------------------------
  step_t <- proc.time()
  .step_header(4, n_steps, "Building model matrix (and competition indices if requested)...")

  if (isTRUE(include_comp_inds)) {
    .vcat("Computing competition indices (radius = ", competition_radius, ").", .indent = 1L)

    X$near_nbr_dist <- NA_real_
    X$near_nbr_num <- NA_real_
    X$avg_nbr_dist <- NA_real_
    X$near_nbr_time <- NA_real_
    X$near_nbr_time_all <- NA_real_
    X$near_nbr_time_dist_ratio <- NA_real_

    if (edge_correction %in% c("none", "truncation")) {
      distance_matrix <- as.matrix(stats::dist(s, method = "euclidean"))
    } else {
      distance_matrix <- toroidal_dist_matrix_optimized(
        s,
        xy_bounds[2] - xy_bounds[1],
        xy_bounds[4] - xy_bounds[3]
      )
    }

    n_pts <- nrow(X)
    for (i in seq_len(n_pts)) {
      if (isTRUE(verbose) && n_pts >= 2000L && (i %% 500L == 0L)) .vcat("...processed ", i, "/", n_pts, .indent = 1L)

      close_points <- unique(which(distance_matrix[i, ] < competition_radius & distance_matrix[i, ] != 0))
      close_times <- X$time[close_points]

      X$near_nbr_dist[i] <- min(distance_matrix[i, ][-i])
      X$near_nbr_num[i] <- length(close_points)
      X$avg_nbr_dist[i] <- if (length(close_points)) mean(distance_matrix[i, close_points]) else min(distance_matrix[i, ][-i])

      nn_idx <- unique(which(distance_matrix[i, ] == X$near_nbr_dist[i]))
      X$near_nbr_time[i] <- X$time[nn_idx][1]
      X$near_nbr_time_all[i] <- if (length(close_points)) mean(close_times) else X$time[nn_idx][1]
      X$near_nbr_time_dist_ratio[i] <- X$near_nbr_time[i] / X$near_nbr_dist[i]
    }
  }

  model_data <- data.frame(size = df$size, X)

  if (edge_correction == "truncation") {
    ax <- xy_bounds[1]; bx <- xy_bounds[2]
    ay <- xy_bounds[3]; by <- xy_bounds[4]
    before_n <- nrow(model_data)
    model_data <- model_data[
      model_data$x > (ax + 15) &
        model_data$x < (bx - 15) &
        model_data$y > (ay + 15) &
        model_data$y < (by - 15),
      , drop = FALSE
    ]
    .vcat("Truncation kept ", nrow(model_data), "/", before_n, " rows.", .indent = 1L)
  }

  if (nrow(model_data) < 2) stop("Not enough observations to train a model after filtering.", call. = FALSE)

  # Keep stable row identities for out-of-sample fold mapping.
  rownames(model_data) <- as.character(seq_len(nrow(model_data)))

  .vcat("Final training rows: ", nrow(model_data), .indent = 1L)
  .vcat("Final feature columns (incl x,y,time): ", ncol(model_data) - 1L, .indent = 1L)
  .vcat("Done in ", .fmt_time(.elapsed_sec(step_t)), ".", .indent = 1L)

  # -------------------------
  # Step 5: tune/fit model
  # -------------------------
  step_t <- proc.time()
  .step_header(5, n_steps, "Fitting model (with optional CV tuning)...")

  metric_set <- if (selection_metric == "rsq") {
    yardstick::metric_set(yardstick::rmse, yardstick::mae, yardstick::rsq)
  } else {
    yardstick::metric_set(yardstick::rmse, yardstick::mae)
  }

  ctrl <- tune::control_grid(
    verbose = FALSE,
    parallel_over = if (isTRUE(parallel)) "resamples" else NULL
  )

  recipe_spec <- recipes::recipe(size ~ ., data = model_data)
  do_tuning <- (cv_folds >= 2L) && (nrow(model_data) >= cv_folds)

  wf_fit <- NULL
  wf_template <- NULL
  fit_engine <- NULL
  engine <- NULL

  if (model_type == "xgboost") {

    spec <- parsnip::boost_tree(
      mode = "regression",
      trees = if (do_tuning) hardhat::tune() else 500,
      min_n = if (do_tuning) hardhat::tune() else 5,
      tree_depth = if (do_tuning) hardhat::tune() else 6,
      learn_rate = if (do_tuning) hardhat::tune() else 0.05,
      loss_reduction = if (do_tuning) hardhat::tune() else 0
    ) %>%
      parsnip::set_engine(
        "xgboost",
        objective = "reg:squarederror",
        seed = seed,
        nthread = engine_threads,
        verbose = 0
      )

    wf0 <- workflows::workflow() %>%
      workflows::add_model(spec) %>%
      workflows::add_recipe(recipe_spec)

    if (do_tuning) {
      set.seed(seed + 101L)
      folds <- rsample::vfold_cv(model_data, v = cv_folds)
      params <- dials::parameters(
        dials::trees(),
        dials::min_n(),
        dials::tree_depth(),
        dials::learn_rate(),
        dials::loss_reduction()
      )
      set.seed(seed + 102L)
      grid <- dials::grid_space_filling(params, size = tuning_grid_size)

      .vcat("foreach backend: ", foreach::getDoParName(),
            " | workers=", foreach::getDoParWorkers(), .indent = 1L)

      tuned <- tune::tune_grid(
        object = wf0,
        resamples = folds,
        grid = grid,
        metrics = metric_set,
        control = ctrl
      )

      best <- tune::select_best(tuned, metric = selection_metric)
      wf_template <- tune::finalize_workflow(wf0, best) # unfitted template
      wf_fit <- parsnip::fit(wf_template, data = model_data)
    } else {
      wf_template <- wf0
      wf_fit <- parsnip::fit(wf0, data = model_data)
    }

    fit_engine <- workflows::extract_fit_engine(wf_fit)
    engine <- "xgboost"

  } else {

    spec <- parsnip::rand_forest(
      mode = "regression",
      trees = if (do_tuning) hardhat::tune() else 500,
      min_n = if (do_tuning) hardhat::tune() else 5
    ) %>%
      parsnip::set_engine("ranger", num.threads = engine_threads, seed = seed)

    wf0 <- workflows::workflow() %>%
      workflows::add_model(spec) %>%
      workflows::add_recipe(recipe_spec)

    if (do_tuning) {
      set.seed(seed + 201L)
      folds <- rsample::vfold_cv(model_data, v = cv_folds)
      params <- dials::parameters(dials::trees(), dials::min_n())
      set.seed(seed + 202L)
      grid <- dials::grid_space_filling(params, size = tuning_grid_size)

      tuned <- tune::tune_grid(
        object = wf0,
        resamples = folds,
        grid = grid,
        metrics = metric_set,
        control = ctrl
      )

      best <- tune::select_best(tuned, metric = selection_metric)
      wf_template <- tune::finalize_workflow(wf0, best)
      wf_fit <- parsnip::fit(wf_template, data = model_data)
    } else {
      wf_template <- wf0
      wf_fit <- parsnip::fit(wf0, data = model_data)
    }

    fit_engine <- workflows::extract_fit_engine(wf_fit)
    engine <- "ranger"
  }

  .vcat("Done in ", .fmt_time(.elapsed_sec(step_t)), ".", .indent = 1L)

  # -------------------------
  # Step 6: finalize object + residual bootstrap + save
  # -------------------------
  step_t <- proc.time()
  .step_header(6, n_steps, "Finalizing output object...")

  preprocessing_recipe <- recipes::prep(recipe_spec, training = model_data, retain = TRUE)
  baked <- recipes::bake(preprocessing_recipe, new_data = model_data)
  baked$size <- NULL
  feature_names <- names(baked)

  # ---- residual bootstrap (optional) ----
  rb_info <- NULL
  if (isTRUE(store_residual_bootstrap)) {

    set.seed(as.integer(residual_seed))

    n_obs <- nrow(model_data)

    # defaults based on FINAL model_data, not `data` (fixes your error)
    if (is.null(resid_bins) || !length(resid_bins) || is.na(resid_bins)) {
      resid_bins <- if (n_obs <= 200) 6L else 10L
    }
    if (is.null(resid_min_per_bin) || !length(resid_min_per_bin) || is.na(resid_min_per_bin)) {
      # small-n friendly default
      resid_min_per_bin <- if (n_obs <= 200) 8L else 15L
    }
    resid_bins <- max(2L, as.integer(resid_bins))
    resid_min_per_bin <- max(3L, as.integer(resid_min_per_bin))

    y_raw <- model_data$size

    # Compute mean predictions to build residual pools
    mu_use <- NULL
    source_used <- residual_source

    can_oos <- (residual_source == "oos") && (cv_folds >= 2L) && (n_obs >= cv_folds)

    if (isTRUE(can_oos)) {
      folds_oos <- rsample::vfold_cv(model_data, v = cv_folds)
      mu_oos <- rep(NA_real_, n_obs)
      names(mu_oos) <- rownames(model_data)

      for (k in seq_along(folds_oos$splits)) {
        sp <- folds_oos$splits[[k]]
        tr <- rsample::analysis(sp)
        te <- rsample::assessment(sp)

        # fit template on training fold
        fit_k <- parsnip::fit(wf_template, data = tr)
        pr_k <- .pred_vec(fit_k, te)

        te_rn <- rownames(te)
        mu_oos[te_rn] <- pr_k
      }

      if (anyNA(mu_oos)) {
        # safety fallback
        source_used <- "in_sample"
        mu_use <- .pred_vec(wf_fit, model_data)
      } else {
        mu_use <- as.numeric(mu_oos)
      }

    } else {
      source_used <- "in_sample"
      mu_use <- .pred_vec(wf_fit, model_data)
    }

    rb <- .build_resid_pools(
      y_raw = y_raw,
      mu_raw = mu_use,
      transform = residual_transform,
      bin_policy = residual_bin_policy,
      n_bins = resid_bins,
      min_per_bin = resid_min_per_bin
    )
    rb$source <- source_used

    rb_info <- list(
      enabled = TRUE,
      source = rb$source,
      transform = rb$transform,
      bin_policy = rb$bin_policy,
      breaks = rb$breaks,
      pools = rb$pools,
      n_bins = rb$n_bins,
      min_per_bin = rb$min_per_bin,
      seed = as.integer(residual_seed)
    )

    .vcat("Residual bootstrap stored (source=", source_used,
          ", transform=", residual_transform,
          ", bins=", resid_bins,
          ", min/bin=", resid_min_per_bin, ").", .indent = 1L)
  }

  mm <- ldmppr_mark_model(
    engine = engine,
    fit_engine = fit_engine,
    recipe = preprocessing_recipe,
    outcome = "size",
    feature_names = feature_names,
    rasters = raster_list,
    info = list(
      model_type = model_type,
      edge_correction = edge_correction,
      include_comp_inds = include_comp_inds,
      competition_radius = competition_radius,
      selection_metric = selection_metric,
      cv_folds = cv_folds,
      tuning_grid_size = tuning_grid_size,
      seed = seed,
      scaled_rasters = scaled_rasters,
      residual_bootstrap = rb_info
    )
  )

  if (isTRUE(save_model)) {
    save_mark_model(mm, save_path)
    .vcat("Saved model to: ", save_path, .indent = 1L)
  }

  .vcat("Done in ", .fmt_time(.elapsed_sec(step_t)), ".", .indent = 1L)
  .vcat("Training complete. Total time: ", .fmt_time(.elapsed_sec(t_total)), ".")

  mm
}
