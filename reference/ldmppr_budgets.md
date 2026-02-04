# Create an optimization budget specification for estimate_process_parameters()

`ldmppr_budgets()` defines per-stage optimization options (budgets) used
by
[`estimate_process_parameters`](https://lanedrew.github.io/ldmppr/reference/estimate_process_parameters.md)
for NLopt via
[`nloptr`](https://astamm.github.io/nloptr/reference/nloptr.html).

## Usage

``` r
ldmppr_budgets(
  global_options = NULL,
  local_budget_first_level = NULL,
  local_budget_refinement_levels = NULL
)
```

## Arguments

- global_options:

  Optional list of NLopt options used for the global stage (only
  relevant when `strategy` uses a global optimizer). Examples:
  `list(maxeval = 2000, maxtime = 10)`.

- local_budget_first_level:

  Optional list of NLopt options used for the local stage at the first
  (coarsest) grid level.

- local_budget_refinement_levels:

  Optional list of NLopt options used for local refinement on subsequent
  (finer) grid levels in multi-resolution strategies. If `NULL`, the
  estimator will fall back to `local_budget_first_level`.

## Value

An object of class `"ldmppr_budgets"`.

## Details

The returned object is an S3 class. Use
[`summary()`](https://rdrr.io/r/base/summary.html) and
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) methods
(if you provide them) to inspect.

## See also

[`ldmppr_grids-class`](https://lanedrew.github.io/ldmppr/reference/ldmppr_grids-class.md)
for methods and details.

## Examples

``` r
b <- ldmppr_budgets(
  global_options = list(maxeval = 150),
  local_budget_first_level = list(maxeval = 300, xtol_rel = 1e-5),
  local_budget_refinement_levels = list(maxeval = 150, xtol_rel = 1e-5)
)
b
#> <ldmppr_budgets>
#>   global_options:
#>     - maxeval: 150
#>     - maxtime: NA
#>   local_budget_first_level:
#>     - maxeval:  300
#>     - maxtime:  NA
#>     - xtol_rel: 1e-05
#>   local_budget_refinement_levels:
#>     - maxeval:  150
#>     - maxtime:  NA
#>     - xtol_rel: 1e-05
```
