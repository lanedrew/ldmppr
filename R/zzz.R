# This file contains the .onLoad() function to handle initialization tasks
# and ensure imports from xgboost and ranger are recognized.
.onLoad <- function(libname, pkgname) {
  # Reference the functions to ensure they are recognized
  xgboost::xgboost
  ranger::ranger
}
