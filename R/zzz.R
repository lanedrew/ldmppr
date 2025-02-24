# This file contains the .onLoad() function to handle initialization tasks
# and ensure imports from xgboost and ranger are recognized.
.onLoad <- function(libname, pkgname) {

  # CRAN OMP THREAD LIMIT
  Sys.setenv("OMP_THREAD_LIMIT" = 1)

  # Reference the functions to ensure they are recognized
  xgboost::xgboost
  ranger::ranger
}
