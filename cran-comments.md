## Resubmission
This is a small patch. In this version of the package we have:

* updated `train_mark_model()` and  `check_model_fit()` and `simulate_mpp()` to include `scaled_rasters` argument to determine if scaling needs to be performed.

* added a new example dataset entitled `medium_example_data` and corresponding raster files.

* updated the `plot_mpp()` function to use the operator `%>%` instead of the `|>` operator to ensure compatibility with older versions of R.

* Bumped version number from 1.0.3 to 1.0.4.

We have not included a reference describing the methods in this package as the reference is not yet published. We will include this reference in the next version of the package.


## R CMD check results

0 errors | 0 warnings | 1 note


## Downstream dependencies
This package has no downstream dependencies.
