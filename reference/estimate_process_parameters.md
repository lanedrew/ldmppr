# Estimate point process parameters using log-likelihood maximization

Estimate spatio-temporal point process parameters by maximizing the
(approximate) full log-likelihood using
[`nloptr`](https://astamm.github.io/nloptr/reference/nloptr.html). For
the self-correcting process, the arrival times must be on \\(0,1)\\ and
can either be supplied directly in `data` as `time`, or constructed from
`size` via the gentle-decay (power-law) mapping
[`power_law_mapping`](https://lanedrew.github.io/ldmppr/reference/power_law_mapping.md)
using `delta` (single fit) or `delta_values` (delta search).

## Usage

``` r
estimate_process_parameters(
  data,
  process = c("self_correcting"),
  x_grid = NULL,
  y_grid = NULL,
  t_grid = NULL,
  upper_bounds = NULL,
  parameter_inits = NULL,
  delta = NULL,
  delta_values = NULL,
  parallel = FALSE,
  num_cores = max(1L, parallel::detectCores() - 1L),
  set_future_plan = FALSE,
  strategy = c("local", "global_local", "multires_global_local"),
  grid_levels = NULL,
  refine_best_delta = TRUE,
  global_algorithm = "NLOPT_GN_CRS2_LM",
  local_algorithm = "NLOPT_LN_BOBYQA",
  global_options = list(maxeval = 150),
  local_options = list(maxeval = 300, xtol_rel = 1e-05, maxtime = NULL),
  global_n_starts = 1L,
  n_starts = 1L,
  jitter_sd = 0.35,
  seed = 1L,
  finite_bounds = NULL,
  verbose = TRUE
)
```

## Arguments

- data:

  a data.frame or matrix. Must contain either columns `(time, x, y)` or
  `(x, y, size)`. If a matrix is provided for delta search, it must have
  column names `c("x","y","size")`.

- process:

  character string specifying the process model. Currently supports
  `"self_correcting"`.

- x_grid, y_grid, t_grid:

  numeric vectors defining the integration grid for \\(x,y,t)\\.

- upper_bounds:

  numeric vector of length 3 giving `c(b_t, b_x, b_y)`.

- parameter_inits:

  (optional) numeric vector of length 8 giving initialization values for
  the model parameters. If `NULL`, sensible defaults are derived from
  `data` and `upper_bounds`.

- delta:

  (optional) numeric scalar used only when `data` contains `(x,y,size)`
  but not `time`.

- delta_values:

  (optional) numeric vector. If supplied, the function fits the model
  for each value of `delta_values` (mapping `size -> time` via
  [`power_law_mapping`](https://lanedrew.github.io/ldmppr/reference/power_law_mapping.md))
  and returns the best fit (lowest objective).

- parallel:

  logical. If `TRUE`, uses furrr/future to parallelize either (a) over
  `delta_values` (when provided) or (b) over multi-start initializations
  (when `delta_values` is `NULL` and `n_starts > 1`).

- num_cores:

  Integer number of workers to use when `set_future_plan = TRUE`.

- set_future_plan:

  `TRUE` or `FALSE`, if `TRUE`, temporarily sets
  `future::plan(multisession, workers = num_cores)` and restores the
  original plan on exit.

- strategy:

  Character string specifying the estimation strategy: - `"local"`:
  single-level local optimization from `parameter_inits`. -
  `"global_local"`: single-level global optimization (from
  `parameter_inits`) followed by local polish. -
  `"multires_global_local"`: multi-resolution fitting over `grid_levels`
  (coarsest level uses global + local; finer levels use local polish
  only).

- grid_levels:

  (optional) list defining the multi-resolution grid schedule when
  `strategy = "multires_global_local"`. Each entry can be a numeric
  vector `c(nx, ny, nt)` or a list with named entries
  `list(nx=..., ny=..., nt=...)`. If `NULL`, uses the supplied
  `(x_grid, y_grid, t_grid)` as a single level.

- refine_best_delta:

  `TRUE` or `FALSE`, if `TRUE` and `delta_values` is supplied, performs
  a final refinement fit at the best delta found using the full
  multi-resolution strategy.

- global_algorithm, local_algorithm:

  character strings specifying the NLopt algorithms to use for the
  global and local optimization stages, respectively.

- global_options, local_options:

  named lists of options to pass to
  [`nloptr::nloptr()`](https://astamm.github.io/nloptr/reference/nloptr.html)
  for the global and local optimization stages, respectively.

- global_n_starts:

  integer number of restarts to use for the global optimization stage.

- n_starts:

  integer number of multi-start initializations to use for the local
  optimization stage.

- jitter_sd:

  numeric standard deviation used to jitter the multi-start
  initializations.

- seed:

  integer random seed used for multi-start initialization jittering.

- finite_bounds:

  (optional) list with components `lb` and `ub` giving finite lower and
  upper bounds for all 8 parameters. Used only when the selected
  optimization algorithms require finite bounds.

- verbose:

  `TRUE` or `FALSE`, if `TRUE`, prints progress messages during fitting.

## Value

an object of class `"ldmppr_fit"` containing the best `nloptr` fit and
(optionally) all fits from a delta search.

## Details

For the self-correcting process, the log-likelihood integral is
approximated using the supplied grid `(x_grid, y_grid, t_grid)` over the
bounded domain `upper_bounds`. When `delta_values` is supplied, this
function performs a grid search over `delta` values, fitting the model
separately for each mapped dataset and selecting the best objective
value.

## References

Møller, J., Ghorbani, M., & Rubak, E. (2016). Mechanistic
spatio-temporal point process models for marked point processes, with a
view to forest stand data. *Biometrics*, 72(3), 687-696.
[doi:10.1111/biom.12466](https://doi.org/10.1111/biom.12466) .

## Examples

``` r
data(small_example_data)

x_grid <- seq(0, 25, length.out = 5)
y_grid <- seq(0, 25, length.out = 5)
t_grid <- seq(0, 1,  length.out = 5)

parameter_inits <- c(1.5, 8.5, .015, 1.5, 3.2, .75, 3, .08)
upper_bounds <- c(1, 25, 25)

fit <- estimate_process_parameters(
  data = small_example_data,
  process = "self_correcting",
  x_grid = x_grid,
  y_grid = y_grid,
  t_grid = t_grid,
  upper_bounds = upper_bounds,
  parameter_inits = parameter_inits,
  delta = 1,
  strategy = "global_local",
  global_algorithm = "NLOPT_GN_CRS2_LM",
  local_algorithm = "NLOPT_LN_BOBYQA",
  global_options = list(maxeval = 150),
  local_options = list(maxeval = 25, xtol_rel = 1e-2),
  verbose = TRUE
)
#> Estimating self-correcting process parameters
#>   Strategy: global_local
#>   Delta: 1
#>   Local optimizer: NLOPT_LN_BOBYQA (maxeval=25)
#>   Global optimizer: NLOPT_GN_CRS2_LM (maxeval=150)
#>   Grid levels: 1
#>   Local starts per level: 1
#>   Global starts (first level only): 1
#>   Parallel: off
#> Step 1/2: Preparing data and objective function...
#>   Prepared 121 points.
#>   Done in 0.0s.
#> Step 2/2: Optimizing parameters...
#> Single level (grid 5x5x5)
#>   Global search: 1 start(s), then local refinement.
#>   Completed in 0.0s.
#>   Best objective: 269.07245
#> Finished. Total time: 0.0s.

coef(fit)
#> [1] 1.5000000 7.6837649 0.0302285 1.5000000 3.2000000 0.7500000 3.0000000
#> [8] 0.0800000
logLik(fit)
#> 'log Lik.' -269.0724 (df=8)

# \donttest{
# Delta-search example (data has x,y,size; time is derived internally for each delta)
fit_delta <- estimate_process_parameters(
  data = small_example_data, # x,y,size
  process = "self_correcting",
  x_grid = x_grid,
  y_grid = y_grid,
  t_grid = t_grid,
  upper_bounds = upper_bounds,
  parameter_inits = parameter_inits,
  delta_values = c(0.35, 0.5, 0.65, 0.9, 1.0),
  parallel = TRUE,
  set_future_plan = TRUE,
  num_cores = 2,
  strategy = "multires_global_local",
  grid_levels = list(
  list(nx = 5, ny = 5, nt = 5),
  list(nx = 8, ny = 8, nt = 8),
  list(nx = 10, ny = 10, nt = 10)
  ),
  global_options = list(maxeval = 100),
  local_options  = list(maxeval = 100, xtol_rel = 1e-3),
  n_starts = 1,
  refine_best_delta = FALSE,
  verbose = FALSE
)
plot(fit_delta)

# }
```
