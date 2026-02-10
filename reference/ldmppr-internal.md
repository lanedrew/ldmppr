# Internal helpers (not part of the public API)

These functions are used internally by ldmppr and are not intended to be
called directly by users.

## Usage

``` r
new_ldmppr_model_check(
  combined_env,
  envs,
  curve_sets,
  sim_metrics = NULL,
  settings = list(),
  call = NULL
)

new_ldmppr_sim(
  process,
  mpp,
  realization,
  params,
  bounds,
  anchor_point,
  thinning,
  edge_correction,
  include_comp_inds,
  competition_radius,
  call = NULL,
  meta = list()
)

new_ldmppr_mark_model(
  engine,
  fit_engine = NULL,
  xgb_raw = NULL,
  recipe = NULL,
  outcome = "size",
  feature_names = NULL,
  rasters = NULL,
  info = list()
)

new_ldmppr_fit(
  process,
  fit,
  fits = NULL,
  mapping = NULL,
  grid = NULL,
  data_summary = NULL,
  data = NULL,
  data_original = NULL,
  engine = "nloptr",
  call = NULL,
  timing = NULL
)

preprocess_new_data(object, new_data)

rehydrate_xgb(object)

as_mark_model(mark_model)

.build_sc_matrix(data, delta = NULL)

.default_sc_param_bounds(txy, upper_bounds)

a %||% b

.require_pkgs(pkgs)

.coerce_training_df(x, delta = NULL, xy_bounds = NULL)

infer_xy_bounds_from_ppp(ppp)

infer_anchor_from_ppp(ppp)

infer_anchor_from_df(df)

resolve_sc_params(process_fit)

resolve_reference_ppp(reference_data, process_fit, xy_bounds)

.as_sc_params(process_fit)

.infer_xy_bounds(process_fit)

.infer_anchor_point(process_fit)

new_ldmppr_budgets(
  global_options,
  local_budget_first_level,
  local_budget_refinement_levels = NULL
)

is_ldmppr_budgets(x)

as_ldmppr_budgets(x, ...)

.validate_ldmppr_budgets(b)

new_ldmppr_grids(levels, upper_bounds, labels = NULL, include_endpoints = TRUE)

is_ldmppr_grids(x)

as_ldmppr_grids(x, ...)

.validate_ldmppr_grids(g)

.ldmppr_make_grid_schedule(
  upper_bounds,
  levels,
  labels = NULL,
  include_endpoints = TRUE
)

infer_rasters_from_mark_model(mm)

infer_scaled_flag_from_mark_model(mm)
```
