## Resubmission
This is a minor update. In this version of the package we have:

* Reworked the documentation to improve clarity and fix minor issues.

* Updated the workflow pipeline to allow passage of S3 class objects between functions to simplify the user experience.
The four main steps of the workflow (fit process model, train mark model, check model fit, simulate realizations) now each have their own S3 class and associated methods and you can pass the
objects forward between functions. See the updated documentation for details.

* Bumped version number from 1.1.0 to 1.1.1.

We have not included a reference describing the methods in this package as the reference is not yet published. The reference is under revision currently and will be included
in the next version of the package.

## R CMD check results

0 errors | 0 warnings | 0 note

## Downstream dependencies
This package has no downstream dependencies.
