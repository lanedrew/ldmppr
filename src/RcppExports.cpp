// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// pdistCVec
double pdistCVec(NumericVector x, NumericVector y);
RcppExport SEXP _ldmppr_pdistCVec(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(pdistCVec(x, y));
    return rcpp_result_gen;
END_RCPP
}
// prodFullCpp
double prodFullCpp(double xgrid, double ygrid, double tgrid, NumericMatrix data, NumericVector params);
RcppExport SEXP _ldmppr_prodFullCpp(SEXP xgridSEXP, SEXP ygridSEXP, SEXP tgridSEXP, SEXP dataSEXP, SEXP paramsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type xgrid(xgridSEXP);
    Rcpp::traits::input_parameter< double >::type ygrid(ygridSEXP);
    Rcpp::traits::input_parameter< double >::type tgrid(tgridSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    rcpp_result_gen = Rcpp::wrap(prodFullCpp(xgrid, ygrid, tgrid, data, params));
    return rcpp_result_gen;
END_RCPP
}
// Ctheta2i
double Ctheta2i(NumericVector xgrid, NumericVector ygrid, double tgrid, NumericMatrix data, NumericVector params, NumericVector bounds);
RcppExport SEXP _ldmppr_Ctheta2i(SEXP xgridSEXP, SEXP ygridSEXP, SEXP tgridSEXP, SEXP dataSEXP, SEXP paramsSEXP, SEXP boundsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type xgrid(xgridSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ygrid(ygridSEXP);
    Rcpp::traits::input_parameter< double >::type tgrid(tgridSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bounds(boundsSEXP);
    rcpp_result_gen = Rcpp::wrap(Ctheta2i(xgrid, ygrid, tgrid, data, params, bounds));
    return rcpp_result_gen;
END_RCPP
}
// CondSumCpp
double CondSumCpp(NumericVector obst, double evalt, NumericVector y);
RcppExport SEXP _ldmppr_CondSumCpp(SEXP obstSEXP, SEXP evaltSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type obst(obstSEXP);
    Rcpp::traits::input_parameter< double >::type evalt(evaltSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(CondSumCpp(obst, evalt, y));
    return rcpp_result_gen;
END_RCPP
}
// CondSumCppR
double CondSumCppR(NumericVector obst, double evalt, LogicalVector y);
RcppExport SEXP _ldmppr_CondSumCppR(SEXP obstSEXP, SEXP evaltSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type obst(obstSEXP);
    Rcpp::traits::input_parameter< double >::type evalt(evaltSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(CondSumCppR(obst, evalt, y));
    return rcpp_result_gen;
END_RCPP
}
// pdistC
NumericVector pdistC(NumericVector evalu, NumericMatrix obsu);
RcppExport SEXP _ldmppr_pdistC(SEXP evaluSEXP, SEXP obsuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type evalu(evaluSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type obsu(obsuSEXP);
    rcpp_result_gen = Rcpp::wrap(pdistC(evalu, obsu));
    return rcpp_result_gen;
END_RCPP
}
// rdistC
NumericVector rdistC(double evalt, NumericVector obst);
RcppExport SEXP _ldmppr_rdistC(SEXP evaltSEXP, SEXP obstSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type evalt(evaltSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type obst(obstSEXP);
    rcpp_result_gen = Rcpp::wrap(rdistC(evalt, obst));
    return rcpp_result_gen;
END_RCPP
}
// Part2FullCpp
double Part2FullCpp(NumericVector xgrid, NumericVector ygrid, NumericVector tgrid, NumericMatrix data, NumericVector params, NumericVector bounds);
RcppExport SEXP _ldmppr_Part2FullCpp(SEXP xgridSEXP, SEXP ygridSEXP, SEXP tgridSEXP, SEXP dataSEXP, SEXP paramsSEXP, SEXP boundsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type xgrid(xgridSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ygrid(ygridSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type tgrid(tgridSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bounds(boundsSEXP);
    rcpp_result_gen = Rcpp::wrap(Part2FullCpp(xgrid, ygrid, tgrid, data, params, bounds));
    return rcpp_result_gen;
END_RCPP
}
// Part1_1FullCpp
double Part1_1FullCpp(NumericMatrix data, NumericVector paramt);
RcppExport SEXP _ldmppr_Part1_1FullCpp(SEXP dataSEXP, SEXP paramtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type paramt(paramtSEXP);
    rcpp_result_gen = Rcpp::wrap(Part1_1FullCpp(data, paramt));
    return rcpp_result_gen;
END_RCPP
}
// Part1_2FullCpp
double Part1_2FullCpp(NumericMatrix data, NumericVector params);
RcppExport SEXP _ldmppr_Part1_2FullCpp(SEXP dataSEXP, SEXP paramsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    rcpp_result_gen = Rcpp::wrap(Part1_2FullCpp(data, params));
    return rcpp_result_gen;
END_RCPP
}
// Part1_3FullCpp
double Part1_3FullCpp(NumericVector xgrid, NumericVector ygrid, NumericVector tgrid, NumericMatrix data, NumericVector params, NumericVector bounds);
RcppExport SEXP _ldmppr_Part1_3FullCpp(SEXP xgridSEXP, SEXP ygridSEXP, SEXP tgridSEXP, SEXP dataSEXP, SEXP paramsSEXP, SEXP boundsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type xgrid(xgridSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ygrid(ygridSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type tgrid(tgridSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bounds(boundsSEXP);
    rcpp_result_gen = Rcpp::wrap(Part1_3FullCpp(xgrid, ygrid, tgrid, data, params, bounds));
    return rcpp_result_gen;
END_RCPP
}
// Part1_4FullCpp
double Part1_4FullCpp(NumericMatrix data, NumericVector paramg);
RcppExport SEXP _ldmppr_Part1_4FullCpp(SEXP dataSEXP, SEXP paramgSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type paramg(paramgSEXP);
    rcpp_result_gen = Rcpp::wrap(Part1_4FullCpp(data, paramg));
    return rcpp_result_gen;
END_RCPP
}
// Part1FullCpp
double Part1FullCpp(NumericVector xgrid, NumericVector ygrid, NumericVector tgrid, NumericMatrix data, NumericVector params, NumericVector bounds);
RcppExport SEXP _ldmppr_Part1FullCpp(SEXP xgridSEXP, SEXP ygridSEXP, SEXP tgridSEXP, SEXP dataSEXP, SEXP paramsSEXP, SEXP boundsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type xgrid(xgridSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ygrid(ygridSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type tgrid(tgridSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bounds(boundsSEXP);
    rcpp_result_gen = Rcpp::wrap(Part1FullCpp(xgrid, ygrid, tgrid, data, params, bounds));
    return rcpp_result_gen;
END_RCPP
}
// full_sc_lhood
double full_sc_lhood(NumericVector xgrid, NumericVector ygrid, NumericVector tgrid, NumericMatrix data, NumericVector params, NumericVector bounds);
RcppExport SEXP _ldmppr_full_sc_lhood(SEXP xgridSEXP, SEXP ygridSEXP, SEXP tgridSEXP, SEXP dataSEXP, SEXP paramsSEXP, SEXP boundsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type xgrid(xgridSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ygrid(ygridSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type tgrid(tgridSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bounds(boundsSEXP);
    rcpp_result_gen = Rcpp::wrap(full_sc_lhood(xgrid, ygrid, tgrid, data, params, bounds));
    return rcpp_result_gen;
END_RCPP
}
// interactionCpp
double interactionCpp(NumericMatrix Hist, NumericVector newp, NumericVector par);
RcppExport SEXP _ldmppr_interactionCpp(SEXP HistSEXP, SEXP newpSEXP, SEXP parSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type Hist(HistSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type newp(newpSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type par(parSEXP);
    rcpp_result_gen = Rcpp::wrap(interactionCpp(Hist, newp, par));
    return rcpp_result_gen;
END_RCPP
}
// interactionCpp_st
NumericVector interactionCpp_st(NumericMatrix data, NumericVector paramg);
RcppExport SEXP _ldmppr_interactionCpp_st(SEXP dataSEXP, SEXP paramgSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type paramg(paramgSEXP);
    rcpp_result_gen = Rcpp::wrap(interactionCpp_st(data, paramg));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ldmppr_pdistCVec", (DL_FUNC) &_ldmppr_pdistCVec, 2},
    {"_ldmppr_prodFullCpp", (DL_FUNC) &_ldmppr_prodFullCpp, 5},
    {"_ldmppr_Ctheta2i", (DL_FUNC) &_ldmppr_Ctheta2i, 6},
    {"_ldmppr_CondSumCpp", (DL_FUNC) &_ldmppr_CondSumCpp, 3},
    {"_ldmppr_CondSumCppR", (DL_FUNC) &_ldmppr_CondSumCppR, 3},
    {"_ldmppr_pdistC", (DL_FUNC) &_ldmppr_pdistC, 2},
    {"_ldmppr_rdistC", (DL_FUNC) &_ldmppr_rdistC, 2},
    {"_ldmppr_Part2FullCpp", (DL_FUNC) &_ldmppr_Part2FullCpp, 6},
    {"_ldmppr_Part1_1FullCpp", (DL_FUNC) &_ldmppr_Part1_1FullCpp, 2},
    {"_ldmppr_Part1_2FullCpp", (DL_FUNC) &_ldmppr_Part1_2FullCpp, 2},
    {"_ldmppr_Part1_3FullCpp", (DL_FUNC) &_ldmppr_Part1_3FullCpp, 6},
    {"_ldmppr_Part1_4FullCpp", (DL_FUNC) &_ldmppr_Part1_4FullCpp, 2},
    {"_ldmppr_Part1FullCpp", (DL_FUNC) &_ldmppr_Part1FullCpp, 6},
    {"_ldmppr_full_sc_lhood", (DL_FUNC) &_ldmppr_full_sc_lhood, 6},
    {"_ldmppr_interactionCpp", (DL_FUNC) &_ldmppr_interactionCpp, 3},
    {"_ldmppr_interactionCpp_st", (DL_FUNC) &_ldmppr_interactionCpp_st, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_ldmppr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
