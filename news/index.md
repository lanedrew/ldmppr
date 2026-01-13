# Changelog

## ldmppr (development version)

## ldmppr 1.1.1

CRAN release: 2026-01-13

- Updated the workflow pipeline to allow passage of S3 class objects
  between functions to simplify the user experience. The four main steps
  of the workflow (fit process model, train mark model, check model fit,
  simulate realizations) now each have their own S3 class and associated
  methods and you can pass the objects forward between functions. See
  the updated documentation for details.

## ldmppr 1.1.0

CRAN release: 2026-01-08

- Introduced S3 classes and associated methods for the 4 main workflow
  steps (fit process model/train mark model/check model fit/simulate
  realizations).

- Replaced the `estimate_parameters_sc()` and
  `estimate_parameters_sc_parallel()` functions with the unified
  [`estimate_process_parameters()`](https://lanedrew.github.io/ldmppr/reference/estimate_process_parameters.md)
  function. We redesigned this function to provide multiple strategies
  for the optimization procedure and refactored the underlying C++ code
  to improve efficiency.

- Removed explicit dependence on the Bundle package and introduced the
  [`save_mark_model()`](https://lanedrew.github.io/ldmppr/reference/ldmppr_mark_model.md)
  and
  [`load_mark_model()`](https://lanedrew.github.io/ldmppr/reference/ldmppr_mark_model.md)
  functions to handle saving and loading trained mark models.

- Updated the small example dataset and example trained mark model to
  reflect changes in the package functions.

## ldmppr 1.0.4

CRAN release: 2025-02-24

- updated
  [`train_mark_model()`](https://lanedrew.github.io/ldmppr/reference/train_mark_model.md)
  and
  [`check_model_fit()`](https://lanedrew.github.io/ldmppr/reference/check_model_fit.md)
  and
  [`simulate_mpp()`](https://lanedrew.github.io/ldmppr/reference/simulate_mpp.md)
  to include `scaled_rasters` argument to determine if scaling needs to
  be performed.

- added a new example dataset entitled `medium_example_data` and
  corresponding raster files.

- updated the
  [`plot_mpp()`](https://lanedrew.github.io/ldmppr/reference/plot_mpp.md)
  function to use the operator `%>%` instead of the `|>` operator to
  ensure compatibility with older versions of R.

## ldmppr 1.0.3

CRAN release: 2024-12-02

## ldmppr 1.0.2

## ldmppr 1.0.1

## ldmppr 1.0.0

## ldmppr 0.1.0

- Initial CRAN submission.
