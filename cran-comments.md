## Resubmission
This is a targeted patch release to address an additional issue reported by CRAN
under gcc-UBSAN.

* Fixed undefined behavior in `src/self_correcting_model.cpp` (left shift of
  negative value) in the spatial hash key function used by `thin_st_fast()`
  and `interaction_st_fast()`.

* Replaced signed left-shift key construction with unsigned-safe packing
  (`uint32_t` -> `uint64_t`) to avoid UB while preserving key behavior.

* Bumped version number from 1.1.2 to 1.1.3.

We have not included a reference describing the methods in this package as the reference is not yet published. The reference is under revision currently and will be included
in the next version of the package.

## R CMD check results

0 errors | 0 warnings | 1 notes

Notes observed locally:

* "Days since last update: 1" (expected for this quick follow-up patch).

## Downstream dependencies
This package has no downstream dependencies.
