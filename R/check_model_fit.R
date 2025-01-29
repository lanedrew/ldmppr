#' Check the fit of estimated self-correcting model on the reference point pattern dataset
#'
#' @description
#' Allows the user to perform global envelope tests for the nonparametric L, F, G, J, E, and V summary functions from the \code{spatstat} package.
#' These tests serve as a goodness of fit measure for the estimated model relative to the reference dataset of interest.
#'
#'
#' @param reference_data a ppp object for the reference dataset.
#' @param t_min minimum value for time.
#' @param t_max maximum value for time.
#' @param sc_params vector of parameter values corresponding to (alpha_1, beta_1, gamma_1, alpha_2, beta_2, alpha_3, beta_3, gamma_3).
#' @param anchor_point vector of (x,y) coordinates of point to condition on.
#' @param raster_list a list of raster objects.
#' @param scaled_rasters `TRUE` or `FALSE` indicating whether the rasters have been scaled.
#' @param mark_model a model object (typically from \code{train_mark_model}).
#' @param xy_bounds a vector of domain bounds (2 for x, 2 for y).
#' @param include_comp_inds `TRUE` or `FALSE` indicating whether to generate and use competition indices as covariates.
#' @param competition_radius distance for competition radius if \code{include_comp_inds} is `TRUE`.
#' @param thinning `TRUE` or `FALSE` indicating whether to use the thinned or unthinned simulated values.
#' @param correction type of correction to apply ("none" or "toroidal").
#' @param n_sim number of simulated datasets to generate.
#' @param save_sims `TRUE` or `FALSE` indicating whether to save and return the simulated datasets.
#' @param verbose `TRUE` or `FALSE` indicating whether to show progress of model checking.
#' @param seed an integer value to set the seed for reproducibility.
#'
#' @return a list containing a combined global envelope test, individual global envelope tests for the L, F, G, J, E, and V functions, and simulated metric values (if specified).
#' @export
#'
#' @details
#' This function relies on the \code{spatstat} package for the calculation of the point pattern metrics
#' and the \code{GET} package for the global envelope tests. The L, F, G, J, E, and V functions are a collection of
#' non-parametric summary statistics that describe the spatial distribution of points and marks in a point pattern.
#' See the documentation for [spatstat.explore::Lest()], [spatstat.explore::Fest()], [spatstat.explore::Gest()],
#' [spatstat.explore::Jest()], [spatstat.explore::Emark()], and [spatstat.explore::Vmark()] for more information.
#' Also, see the [GET::global_envelope_test()] function for more information on the global envelope tests.
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
#' # Load the small example data
#' data(small_example_data)
#'
#' # Load the example mark model that previously was trained on the small example data
#' file_path <- system.file("extdata", "example_mark_model.rds", package = "ldmppr")
#' mark_model <- bundle::unbundle(readRDS(file_path))
#'
#' # Load the raster files
#' raster_paths <- list.files(system.file("extdata", package = "ldmppr"),
#'                            pattern = "\\.tif$", full.names = TRUE)
#' rasters <- lapply(raster_paths, terra::rast)
#'
#' # Scale the rasters
#' scaled_raster_list <- scale_rasters(rasters)
#'
#' # Generate the reference pattern
#' reference_data <- generate_mpp(
#'   locations = small_example_data[, c("x", "y")],
#'   marks = small_example_data$size,
#'   xy_bounds = c(0, 25, 0, 25)
#' )
#'
#' # Define an anchor point
#' M_n <- as.matrix(small_example_data[1, c("x", "y")])
#'
#' # Specify the estimated parameters of the self-correcting process
#' # Note: These would generally be estimated using estimate_parameters_sc
#' # or estimate_parameters_sc_parallel. These values are estimates from
#' # the small_example_data generating script.
#' estimated_parameters <- c(
#'   1.42936311, 8.59251417, 0.02153924, 1.89763856,
#'   2.33256061, 1.09522235, 2.66250000, 0.16499789
#' )
#'
#' # Check the model fit
#' example_model_fit <- check_model_fit(
#'   reference_data = reference_data,
#'   t_min = 0,
#'   t_max = 1,
#'   sc_params = estimated_parameters,
#'   anchor_point = M_n,
#'   raster_list = scaled_raster_list,
#'   scaled_rasters = TRUE,
#'   mark_model = mark_model,
#'   xy_bounds = c(0, 25, 0, 25),
#'   include_comp_inds = TRUE,
#'   thinning = TRUE,
#'   correction = "none",
#'   competition_radius = 10,
#'   n_sim = 100,
#'   save_sims = FALSE,
#'   verbose = TRUE,
#'   seed = 90210
#' )
#'
#' plot(example_model_fit$combined_env)
#'}
check_model_fit <- function(reference_data,
                            t_min = 0,
                            t_max = 1,
                            sc_params = NULL,
                            anchor_point = NULL,
                            raster_list = NULL,
                            scaled_rasters = FALSE,
                            mark_model = NULL,
                            xy_bounds = NULL,
                            include_comp_inds = FALSE,
                            thinning = TRUE,
                            correction = "none",
                            competition_radius = 15,
                            n_sim = 2500,
                            save_sims = TRUE,
                            verbose = TRUE,
                            seed = 0) {
  # Check the arguments
  if (!spatstat.geom::is.ppp(reference_data)) stop("Provide a ppp object containing the reference data pattern for the reference_data argument.", .call = FALSE)
  if (t_min < 0 | t_min >= t_max | is.null(t_min)) stop("Provide a value greater than 0 and less than t_max for the t_min argument.", .call = FALSE)
  if (t_max > 1 | t_min >= t_max | is.null(t_max)) stop("Provide a value greater than t_min and less than 1 for the t_max argument.", .call = FALSE)
  if (length(sc_params) != 8 | anyNA(sc_params) | any(sc_params[2:8] < 0)) stop("Provide a valid set of parameter values for the sc_params argument.", .call = FALSE)
  if (length(anchor_point) != 2) stop("Provide a vector of (x,y) coordinates for the anchor_point argument.", .call = FALSE)
  if (is.null(mark_model)) stop("Provide an unbundled mark model for the mark_model argument.", .call = FALSE)
  if (is.null(xy_bounds) | !(length(xy_bounds) == 4)) stop("Provide (x,y) bounds in the form (a_x, b_x, a_y, b_y) for the xy_bounds argument.", .call = FALSE)
  if (xy_bounds[2] < xy_bounds[1] | xy_bounds[4] < xy_bounds[3]) stop("Provide (x,y) bounds in the form (a_x, b_x, a_y, b_y) for the xy_bounds argument.", .call = FALSE)
  if (!correction %in% c("none", "toroidal")) stop("Provide a valid correction type for the correction argument.", .call = FALSE)
  if (include_comp_inds == TRUE & (is.null(competition_radius) | competition_radius < 0)) stop("Provide the desired radius for competition_indices argument.", .call = FALSE)
  if (n_sim < 1) stop("Provide a positive integer value for the n_sim argument.", .call = FALSE)
  if (!is.logical(save_sims)) stop("Provide a logical value for the save_sims argument.", .call = FALSE)
  if (!is.logical(verbose)) stop("Provide a logical value for the verbose argument.", .call = FALSE)
  if (!is.numeric(seed) | seed < 0) stop("Provide a positive integer value for the seed argument.", .call = FALSE)
  if (!is.logical(scaled_rasters)) stop("Provide a logical value for the scaled_rasters argument.", .call = FALSE)

  # Set the seed
  set.seed(seed)

  # Obtain the radius to use in obtaining the various estimates using the reference data
  K_ref <- spatstat.explore::Kest(spatstat.geom::unmark(reference_data))
  d <- K_ref$r
  d_length <- base::length(d)

  # Create storage for the simulation results for the various metrics
  K_PP <- base::matrix(0, nrow = d_length, ncol = n_sim)
  F_PP <- base::matrix(0, nrow = d_length, ncol = n_sim)
  G_PP <- base::matrix(0, nrow = d_length, ncol = n_sim)
  J_PP <- base::matrix(0, nrow = d_length, ncol = n_sim)
  E_PP <- base::matrix(0, nrow = d_length, ncol = n_sim)
  V_PP <- base::matrix(0, nrow = d_length, ncol = n_sim)
  n_real <- base::numeric()

  # Scale the rasters if necessary
  if (scaled_rasters == FALSE) {
    raster_list <- scale_rasters(raster_list)
  }

  if (verbose) {
    # Create the progress bar
    pb <- progress::progress_bar$new(
      format = "Data Simulation Progress: [:bar] :percent in :elapsed, ETA: :eta",
      total = n_sim, clear = FALSE, width = 100
    )

    base::message("Beginning Data Simulations...")
    pb$tick(0)
    base::Sys.sleep(3)

    for (j in 1:n_sim) {
      # Simulate a dataset and predict the marks
      if (thinning) {
        sim_j <- simulate_sc(
          t_min = t_min,
          t_max = t_max,
          sc_params = sc_params,
          anchor_point = anchor_point,
          xy_bounds = xy_bounds
        )$thinned
      } else {
        sim_j <- simulate_sc(
          t_min = t_min,
          t_max = t_max,
          sc_params = sc_params,
          anchor_point = anchor_point,
          xy_bounds = xy_bounds
        )$unthinned
      }
      pred_marks_j <- predict_marks(
        sim_realization = sim_j,
        raster_list = raster_list,
        scaled_rasters = TRUE,
        mark_model = mark_model,
        xy_bounds = xy_bounds,
        include_comp_inds = include_comp_inds,
        competition_radius = competition_radius,
        correction = correction
      )

      # Generate a ppp object from the locations and marks
      PP_xy <- generate_mpp(locations = sim_j, marks = pred_marks_j, xy_bounds = xy_bounds)

      # Calculate the point pattern metrics of interest
      K_PP[, j] <- spatstat.explore::Kest(spatstat.geom::unmark(PP_xy), correction = "isotropic", r = d)$iso
      F_PP[, j] <- spatstat.explore::Fest(spatstat.geom::unmark(PP_xy), correction = "rs", r = d)$rs
      G_PP[, j] <- spatstat.explore::Gest(spatstat.geom::unmark(PP_xy), correction = "rs", r = d)$rs
      J_PP[, j] <- spatstat.explore::Jest(spatstat.geom::unmark(PP_xy), correction = "rs", r = d)$rs - 1
      E_PP[, j] <- spatstat.explore::Emark(PP_xy, correction = "isotropic", r = d)$iso
      V_PP[, j] <- spatstat.explore::Vmark(PP_xy, correction = "isotropic", r = d)$iso

      pb$tick()
      base::Sys.sleep(1 / 100)
    }
    base::message("Data Simulations Complete...")
  } else {
    for (j in 1:n_sim) {
      # Simulate a dataset and predict the marks
      if (thinning) {
        # Simulate a dataset using the estimated parameters
        sim_j <- simulate_sc(
          t_min = t_min,
          t_max = t_max,
          sc_params = sc_params,
          anchor_point = anchor_point,
          xy_bounds = xy_bounds
        )$thinned
      } else {
        # Simulate a dataset using the estimated parameters
        sim_j <- simulate_sc(
          t_min = t_min,
          t_max = t_max,
          sc_params = sc_params,
          anchor_point = anchor_point,
          xy_bounds = xy_bounds
        )$unthinned
      }

      # Predict the marks using the provided mark model
      pred_marks_j <- predict_marks(
        sim_realization = sim_j,
        raster_list = raster_list,
        scaled_rasters = TRUE,
        mark_model = mark_model,
        xy_bounds = xy_bounds,
        include_comp_inds = include_comp_inds,
        competition_radius = competition_radius,
        correction = correction
      )

      # Generate a ppp object from the locations and marks
      PP_xy <- generate_mpp(locations = sim_j, marks = pred_marks_j, xy_bounds = xy_bounds)
      n_real[j] <- base::length(sim_j[, 1])

      # Calculate the point pattern metrics of interest
      K_PP[, j] <- spatstat.explore::Kest(spatstat.geom::unmark(PP_xy), correction = "isotropic", r = d)$iso
      F_PP[, j] <- spatstat.explore::Fest(spatstat.geom::unmark(PP_xy), correction = "rs", r = d)$rs
      G_PP[, j] <- spatstat.explore::Gest(spatstat.geom::unmark(PP_xy), correction = "rs", r = d)$rs
      J_PP[, j] <- spatstat.explore::Jest(spatstat.geom::unmark(PP_xy), correction = "rs", r = d)$rs - 1
      E_PP[, j] <- spatstat.explore::Emark(PP_xy, correction = "isotropic", r = d)$iso
      V_PP[, j] <- spatstat.explore::Vmark(PP_xy, correction = "isotropic", r = d)$iso
    }
  }

  if (save_sims) {
    sim_list <- base::list(
      Ksim = K_PP,
      Fsim = F_PP,
      Gsim = G_PP,
      Jsim = J_PP,
      Esim = E_PP,
      Vsim = V_PP,
      n_per = n_real
    )
  }


  # Obtain the global envelope test for the L function
  C_ref_L <- GET::create_curve_set(base::list(
    r = d,
    obs = sqrt(K_ref$iso / pi) - d,
    theo = sqrt(K_ref$theo / pi) - d,
    sim_m = sqrt(K_PP / pi) - d
  ))
  r_envL <- GET::global_envelope_test(C_ref_L, type = "rank")

  # Obtain the global envelope test for the F function
  F_val <- base::max(base::apply(F_PP, 2, function(x) base::min(base::which(x >= 1))))
  C_ref_F <- GET::create_curve_set(base::list(
    r = d[1:F_val],
    obs = spatstat.explore::Fest(spatstat.geom::unmark(reference_data),
      r = d[1:F_val]
    )$rs,
    theo = spatstat.explore::Fest(spatstat.geom::unmark(reference_data),
      r = d[1:F_val]
    )$theo,
    sim_m = F_PP[1:F_val, ]
  ))
  r_envF <- GET::global_envelope_test(C_ref_F, type = "rank")

  # Obtain the global envelope test for the G function
  G_val <- base::max(base::apply(G_PP, 2, function(x) base::min(base::which(x >= 1))))
  C_ref_G <- GET::create_curve_set(base::list(
    r = d[1:G_val],
    obs = spatstat.explore::Gest(spatstat.geom::unmark(reference_data),
      r = d[1:G_val]
    )$rs,
    theo = spatstat.explore::Gest(spatstat.geom::unmark(reference_data),
      r = d[1:G_val]
    )$theo,
    sim_m = G_PP[1:G_val, ]
  ))
  r_envG <- GET::global_envelope_test(C_ref_G, type = "rank")

  # Obtain the global envelope test for the J function
  J_val <- base::min(c(
    base::min(base::apply(F_PP, 2, function(x) base::sum(x < 1, na.rm = TRUE))),
    base::min(base::apply(G_PP, 2, function(x) base::sum(x < 1, na.rm = TRUE))),
    base::sum(!base::is.na(spatstat.explore::Jest(spatstat.geom::unmark(reference_data),
      r = d
    )$rs))
  ))
  J_scale <- base::max(base::apply(J_PP[1:J_val, ], 1, base::max))

  if (any(is.infinite(J_PP[1:J_val, ]) | is.na(J_PP[1:J_val, ]))) {
    warning("J_PP contains Inf or NA values.")
  }

  C_ref_J <- GET::create_curve_set(base::list(
    r = d[1:J_val],
    obs = (spatstat.explore::Jest(spatstat.geom::unmark(reference_data),
      r = d[1:J_val]
    )$rs - 1) / J_scale,
    theo = spatstat.explore::Jest(spatstat.geom::unmark(reference_data),
      r = d[1:J_val]
    )$theo - 1,
    sim_m = (J_PP[1:J_val, ]) / J_scale
  ))
  r_envJ <- GET::global_envelope_test(C_ref_J, type = "rank")

  # Obtain the global envelope test for the E function
  C_ref_E <- GET::create_curve_set(base::list(
    r = d,
    obs = spatstat.explore::Emark(reference_data, correction = "isotropic", r = d)$iso,
    theo = spatstat.explore::Emark(reference_data, correction = "isotropic", r = d)$theo,
    sim_m = E_PP
  ))
  r_envE <- GET::global_envelope_test(C_ref_E, type = "rank")

  # Obtain the global envelope test for the V function
  C_ref_V <- GET::create_curve_set(base::list(
    r = d,
    obs = spatstat.explore::Vmark(reference_data, correction = "isotropic", r = d)$iso,
    theo = spatstat.explore::Vmark(reference_data, correction = "isotropic", r = d)$theo,
    sim_m = V_PP
  ))
  r_envV <- GET::global_envelope_test(C_ref_V, type = "rank")

  # Obtain the global envelope test for the combined L, F, G, and J functions
  r_envComb <- GET::global_envelope_test(
    curve_sets = list(
      L = C_ref_L,
      F = C_ref_F,
      G = C_ref_G,
      J = C_ref_J,
      E = C_ref_E,
      V = C_ref_V
    ),
    type = "rank"
  )

  if (save_sims) {
    results_list <- base::list(
      global_envelope_tests = list(
        combined_env = r_envComb,
        curve_sets = list(
          L = C_ref_L,
          F = C_ref_F,
          G = C_ref_G,
          J = C_ref_J,
          E = C_ref_E,
          V = C_ref_V
        ),
        L_env = r_envL,
        F_env = r_envF,
        G_env = r_envG,
        J_env = r_envJ,
        E_env = r_envE,
        V_env = r_envV
      ),
      sim_metric_values = sim_list
    )
  } else {
    results_list <- base::list(
      combined_env = r_envComb,
      curve_sets = list(
        L = C_ref_L,
        F = C_ref_F,
        G = C_ref_G,
        J = C_ref_J,
        E = C_ref_E,
        V = C_ref_V
      ),
      L_env = r_envL,
      F_env = r_envF,
      G_env = r_envG,
      J_env = r_envJ,
      E_env = r_envE,
      V_env = r_envV
    )
  }

  return(results_list)
}
