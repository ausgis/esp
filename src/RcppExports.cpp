// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// QuantileDisc
Rcpp::IntegerVector QuantileDisc(const arma::vec& x, int n);
RcppExport SEXP _esp_QuantileDisc(SEXP xSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(QuantileDisc(x, n));
    return rcpp_result_gen;
END_RCPP
}
// ArmaRunique
arma::ivec ArmaRunique(const arma::ivec& x);
RcppExport SEXP _esp_ArmaRunique(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::ivec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(ArmaRunique(x));
    return rcpp_result_gen;
END_RCPP
}
// ArmaDummyVec
arma::mat ArmaDummyVec(const arma::vec& vec);
RcppExport SEXP _esp_ArmaDummyVec(SEXP vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type vec(vecSEXP);
    rcpp_result_gen = Rcpp::wrap(ArmaDummyVec(vec));
    return rcpp_result_gen;
END_RCPP
}
// ArmaDummyMat
arma::mat ArmaDummyMat(const arma::imat& mat);
RcppExport SEXP _esp_ArmaDummyMat(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::imat& >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(ArmaDummyMat(mat));
    return rcpp_result_gen;
END_RCPP
}
// InteractionType
std::string InteractionType(double qv12, double qv1, double qv2);
RcppExport SEXP _esp_InteractionType(SEXP qv12SEXP, SEXP qv1SEXP, SEXP qv2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type qv12(qv12SEXP);
    Rcpp::traits::input_parameter< double >::type qv1(qv1SEXP);
    Rcpp::traits::input_parameter< double >::type qv2(qv2SEXP);
    rcpp_result_gen = Rcpp::wrap(InteractionType(qv12, qv1, qv2));
    return rcpp_result_gen;
END_RCPP
}
// GenTau
Rcpp::NumericVector GenTau(double n);
RcppExport SEXP _esp_GenTau(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(GenTau(n));
    return rcpp_result_gen;
END_RCPP
}
// ComputeR2
double ComputeR2(const arma::vec& y, const arma::vec& y_pred);
RcppExport SEXP _esp_ComputeR2(SEXP ySEXP, SEXP y_predSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y_pred(y_predSEXP);
    rcpp_result_gen = Rcpp::wrap(ComputeR2(y, y_pred));
    return rcpp_result_gen;
END_RCPP
}
// ComputeLocalR2
double ComputeLocalR2(const arma::vec& y, const arma::vec& y_pred, double Ymean);
RcppExport SEXP _esp_ComputeLocalR2(SEXP ySEXP, SEXP y_predSEXP, SEXP YmeanSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y_pred(y_predSEXP);
    Rcpp::traits::input_parameter< double >::type Ymean(YmeanSEXP);
    rcpp_result_gen = Rcpp::wrap(ComputeLocalR2(y, y_pred, Ymean));
    return rcpp_result_gen;
END_RCPP
}
// SLMQ
Rcpp::NumericVector SLMQ(const arma::vec& Y, const arma::mat& FitY);
RcppExport SEXP _esp_SLMQ(SEXP YSEXP, SEXP FitYSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type FitY(FitYSEXP);
    rcpp_result_gen = Rcpp::wrap(SLMQ(Y, FitY));
    return rcpp_result_gen;
END_RCPP
}
// SLMLocalQ
Rcpp::List SLMLocalQ(const arma::mat& FitY, const arma::vec& Y, const arma::ivec& Zones);
RcppExport SEXP _esp_SLMLocalQ(SEXP FitYSEXP, SEXP YSEXP, SEXP ZonesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type FitY(FitYSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type Zones(ZonesSEXP);
    rcpp_result_gen = Rcpp::wrap(SLMLocalQ(FitY, Y, Zones));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_esp_QuantileDisc", (DL_FUNC) &_esp_QuantileDisc, 2},
    {"_esp_ArmaRunique", (DL_FUNC) &_esp_ArmaRunique, 1},
    {"_esp_ArmaDummyVec", (DL_FUNC) &_esp_ArmaDummyVec, 1},
    {"_esp_ArmaDummyMat", (DL_FUNC) &_esp_ArmaDummyMat, 1},
    {"_esp_InteractionType", (DL_FUNC) &_esp_InteractionType, 3},
    {"_esp_GenTau", (DL_FUNC) &_esp_GenTau, 1},
    {"_esp_ComputeR2", (DL_FUNC) &_esp_ComputeR2, 2},
    {"_esp_ComputeLocalR2", (DL_FUNC) &_esp_ComputeLocalR2, 3},
    {"_esp_SLMQ", (DL_FUNC) &_esp_SLMQ, 2},
    {"_esp_SLMLocalQ", (DL_FUNC) &_esp_SLMLocalQ, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_esp(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
