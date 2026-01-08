## Resubmission
This is a minor update. In this version of the package we have:

* Introduced S3 classes and associated methods for the 4 main workflow steps (fit process model/train mark model/check model fit/simulate realizations).

* Replaced the `estimate_parameters_sc()` and `estimate_parameters_sc_parallel()` functions with the unified `estimate_process_parameters()` function. We redesigned this function to provide
multiple strategies for the optimization procedure and refactored the underlying C++ code to improve efficiency.

* Removed explicit dependence on the Bundle package and introduced the `save_mark_model()` and `load_mark_model()` functions to handle saving and loading trained mark models.

* Updated the small example dataset and example trained mark model to reflect changes in the package functions.

* Bumped version number from 1.0.4 to 1.1.0.

We have not included a reference describing the methods in this package as the reference is not yet published. The reference is under revision currently and will be included
in the next version of the package.

## R CMD check results

0 errors | 0 warnings | 0 note

## Downstream dependencies
This package has no downstream dependencies.
