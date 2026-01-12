#' @useDynLib ldmppr
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"
NULL

utils::globalVariables("size")
utils::globalVariables("time")
utils::globalVariables(c("x", "y", "marks"))

# Various imports
#' @importFrom stats predict
#' @importFrom dplyr select
#' @importFrom stats coef
