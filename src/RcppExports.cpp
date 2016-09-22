// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// ggrscore
List ggrscore(NumericVector x, std::string score, int iter);
RcppExport SEXP dfphase1_ggrscore(SEXP xSEXP, SEXP scoreSEXP, SEXP iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< std::string >::type score(scoreSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    __result = Rcpp::wrap(ggrscore(x, score, iter));
    return __result;
END_RCPP
}
// ggdepthranks
List ggdepthranks(NumericVector x, int L);
RcppExport SEXP dfphase1_ggdepthranks(SEXP xSEXP, SEXP LSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    __result = Rcpp::wrap(ggdepthranks(x, L));
    return __result;
END_RCPP
}
// ggclassicmshewhart
List ggclassicmshewhart(NumericVector x, std::string stat, int L);
RcppExport SEXP dfphase1_ggclassicmshewhart(SEXP xSEXP, SEXP statSEXP, SEXP LSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< std::string >::type stat(statSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    __result = Rcpp::wrap(ggclassicmshewhart(x, stat, L));
    return __result;
END_RCPP
}
// ggscore2mshewhart
List ggscore2mshewhart(NumericVector x, std::string stat, int L);
RcppExport SEXP dfphase1_ggscore2mshewhart(SEXP xSEXP, SEXP statSEXP, SEXP LSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< std::string >::type stat(statSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    __result = Rcpp::wrap(ggscore2mshewhart(x, stat, L));
    return __result;
END_RCPP
}
// ggglrchart
List ggglrchart(NumericVector x, bool onlymean, int L);
RcppExport SEXP dfphase1_ggglrchart(SEXP xSEXP, SEXP onlymeanSEXP, SEXP LSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< bool >::type onlymean(onlymeanSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    __result = Rcpp::wrap(ggglrchart(x, onlymean, L));
    return __result;
END_RCPP
}
// ggxbars
List ggxbars(NumericMatrix x, bool aggr_with_mean, int L);
RcppExport SEXP dfphase1_ggxbars(SEXP xSEXP, SEXP aggr_with_meanSEXP, SEXP LSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< bool >::type aggr_with_mean(aggr_with_meanSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    __result = Rcpp::wrap(ggxbars(x, aggr_with_mean, L));
    return __result;
END_RCPP
}
// ggxbarsall
NumericMatrix ggxbarsall(int n, int m, bool aggr_with_mean, int rep);
RcppExport SEXP dfphase1_ggxbarsall(SEXP nSEXP, SEXP mSEXP, SEXP aggr_with_meanSEXP, SEXP repSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< bool >::type aggr_with_mean(aggr_with_meanSEXP);
    Rcpp::traits::input_parameter< int >::type rep(repSEXP);
    __result = Rcpp::wrap(ggxbarsall(n, m, aggr_with_mean, rep));
    return __result;
END_RCPP
}
// ggrank
List ggrank(NumericMatrix x);
RcppExport SEXP dfphase1_ggrank(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    __result = Rcpp::wrap(ggrank(x));
    return __result;
END_RCPP
}
// ggrankall
NumericMatrix ggrankall(int n, int m, int rep);
RcppExport SEXP dfphase1_ggrankall(SEXP nSEXP, SEXP mSEXP, SEXP repSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type rep(repSEXP);
    __result = Rcpp::wrap(ggrankall(n, m, rep));
    return __result;
END_RCPP
}