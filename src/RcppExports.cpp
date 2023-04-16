// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// Indiv_Score_Test_meta
List Indiv_Score_Test_meta(arma::vec U, arma::vec V);
RcppExport SEXP _MetaSTAAR_Indiv_Score_Test_meta(SEXP USEXP, SEXP VSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type U(USEXP);
    Rcpp::traits::input_parameter< arma::vec >::type V(VSEXP);
    rcpp_result_gen = Rcpp::wrap(Indiv_Score_Test_meta(U, V));
    return rcpp_result_gen;
END_RCPP
}
// MetaSTAAR_O_SMMAT
arma::vec MetaSTAAR_O_SMMAT(arma::vec x, arma::mat Cov, arma::mat weights_B, arma::mat weights_S, arma::mat weights_A, arma::vec mac, int mac_thres);
RcppExport SEXP _MetaSTAAR_MetaSTAAR_O_SMMAT(SEXP xSEXP, SEXP CovSEXP, SEXP weights_BSEXP, SEXP weights_SSEXP, SEXP weights_ASEXP, SEXP macSEXP, SEXP mac_thresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Cov(CovSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type weights_B(weights_BSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type weights_S(weights_SSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type weights_A(weights_ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mac(macSEXP);
    Rcpp::traits::input_parameter< int >::type mac_thres(mac_thresSEXP);
    rcpp_result_gen = Rcpp::wrap(MetaSTAAR_O_SMMAT(x, Cov, weights_B, weights_S, weights_A, mac, mac_thres));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MetaSTAAR_Indiv_Score_Test_meta", (DL_FUNC) &_MetaSTAAR_Indiv_Score_Test_meta, 2},
    {"_MetaSTAAR_MetaSTAAR_O_SMMAT", (DL_FUNC) &_MetaSTAAR_MetaSTAAR_O_SMMAT, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_MetaSTAAR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
