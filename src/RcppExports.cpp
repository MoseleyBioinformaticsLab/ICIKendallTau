// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// pairwiseComparisons
DataFrame pairwiseComparisons(CharacterVector entries, int n_cores, bool include_self);
RcppExport SEXP _ICIKendallTau_pairwiseComparisons(SEXP entriesSEXP, SEXP n_coresSEXP, SEXP include_selfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type entries(entriesSEXP);
    Rcpp::traits::input_parameter< int >::type n_cores(n_coresSEXP);
    Rcpp::traits::input_parameter< bool >::type include_self(include_selfSEXP);
    rcpp_result_gen = Rcpp::wrap(pairwiseComparisons(entries, n_cores, include_self));
    return rcpp_result_gen;
END_RCPP
}
// sortedIndex
IntegerVector sortedIndex(NumericVector x);
RcppExport SEXP _ICIKendallTau_sortedIndex(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(sortedIndex(x));
    return rcpp_result_gen;
END_RCPP
}
// compare_self
IntegerVector compare_self(NumericVector x);
RcppExport SEXP _ICIKendallTau_compare_self(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(compare_self(x));
    return rcpp_result_gen;
END_RCPP
}
// compare_both
IntegerVector compare_both(IntegerVector x, IntegerVector y);
RcppExport SEXP _ICIKendallTau_compare_both(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(compare_both(x, y));
    return rcpp_result_gen;
END_RCPP
}
// which_notzero
IntegerVector which_notzero(IntegerVector x);
RcppExport SEXP _ICIKendallTau_which_notzero(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(which_notzero(x));
    return rcpp_result_gen;
END_RCPP
}
// kendall_discordant
int kendall_discordant(IntegerVector x, IntegerVector y);
RcppExport SEXP _ICIKendallTau_kendall_discordant(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(kendall_discordant(x, y));
    return rcpp_result_gen;
END_RCPP
}
// count_rank_tie
NumericVector count_rank_tie(IntegerVector ranks);
RcppExport SEXP _ICIKendallTau_count_rank_tie(SEXP ranksSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type ranks(ranksSEXP);
    rcpp_result_gen = Rcpp::wrap(count_rank_tie(ranks));
    return rcpp_result_gen;
END_RCPP
}
// ici_kt
NumericVector ici_kt(NumericVector x, NumericVector y, String perspective, String alternative, bool continuity, String output);
RcppExport SEXP _ICIKendallTau_ici_kt(SEXP xSEXP, SEXP ySEXP, SEXP perspectiveSEXP, SEXP alternativeSEXP, SEXP continuitySEXP, SEXP outputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< String >::type perspective(perspectiveSEXP);
    Rcpp::traits::input_parameter< String >::type alternative(alternativeSEXP);
    Rcpp::traits::input_parameter< bool >::type continuity(continuitySEXP);
    Rcpp::traits::input_parameter< String >::type output(outputSEXP);
    rcpp_result_gen = Rcpp::wrap(ici_kt(x, y, perspective, alternative, continuity, output));
    return rcpp_result_gen;
END_RCPP
}
// ici_kt_pairs
NumericVector ici_kt_pairs(NumericVector x, NumericVector y, String perspective, String alternative, String output);
RcppExport SEXP _ICIKendallTau_ici_kt_pairs(SEXP xSEXP, SEXP ySEXP, SEXP perspectiveSEXP, SEXP alternativeSEXP, SEXP outputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< String >::type perspective(perspectiveSEXP);
    Rcpp::traits::input_parameter< String >::type alternative(alternativeSEXP);
    Rcpp::traits::input_parameter< String >::type output(outputSEXP);
    rcpp_result_gen = Rcpp::wrap(ici_kt_pairs(x, y, perspective, alternative, output));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ICIKendallTau_pairwiseComparisons", (DL_FUNC) &_ICIKendallTau_pairwiseComparisons, 3},
    {"_ICIKendallTau_sortedIndex", (DL_FUNC) &_ICIKendallTau_sortedIndex, 1},
    {"_ICIKendallTau_compare_self", (DL_FUNC) &_ICIKendallTau_compare_self, 1},
    {"_ICIKendallTau_compare_both", (DL_FUNC) &_ICIKendallTau_compare_both, 2},
    {"_ICIKendallTau_which_notzero", (DL_FUNC) &_ICIKendallTau_which_notzero, 1},
    {"_ICIKendallTau_kendall_discordant", (DL_FUNC) &_ICIKendallTau_kendall_discordant, 2},
    {"_ICIKendallTau_count_rank_tie", (DL_FUNC) &_ICIKendallTau_count_rank_tie, 1},
    {"_ICIKendallTau_ici_kt", (DL_FUNC) &_ICIKendallTau_ici_kt, 6},
    {"_ICIKendallTau_ici_kt_pairs", (DL_FUNC) &_ICIKendallTau_ici_kt_pairs, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_ICIKendallTau(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
