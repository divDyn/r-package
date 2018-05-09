// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// Counts
NumericMatrix Counts(NumericVector tax, NumericVector bin);
RcppExport SEXP _divDyn_Counts(SEXP taxSEXP, SEXP binSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type tax(taxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bin(binSEXP);
    rcpp_result_gen = Rcpp::wrap(Counts(tax, bin));
    return rcpp_result_gen;
END_RCPP
}
// CRbinwise
NumericMatrix CRbinwise(NumericVector binVar, int quota);
RcppExport SEXP _divDyn_CRbinwise(SEXP binVarSEXP, SEXP quotaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type binVar(binVarSEXP);
    Rcpp::traits::input_parameter< int >::type quota(quotaSEXP);
    rcpp_result_gen = Rcpp::wrap(CRbinwise(binVar, quota));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_divDyn_Counts", (DL_FUNC) &_divDyn_Counts, 2},
    {"_divDyn_CRbinwise", (DL_FUNC) &_divDyn_CRbinwise, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_divDyn(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
