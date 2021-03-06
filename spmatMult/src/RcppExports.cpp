// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;

// rcpparma_hello_world
arma::mat rcpparma_hello_world();
RcppExport SEXP _spmatMult_rcpparma_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpparma_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_outerproduct
arma::mat rcpparma_outerproduct(const arma::colvec& x);
RcppExport SEXP _spmatMult_rcpparma_outerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_outerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_innerproduct
double rcpparma_innerproduct(const arma::colvec& x);
RcppExport SEXP _spmatMult_rcpparma_innerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_innerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_bothproducts
Rcpp::List rcpparma_bothproducts(const arma::colvec& x);
RcppExport SEXP _spmatMult_rcpparma_bothproducts(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_bothproducts(x));
    return rcpp_result_gen;
END_RCPP
}
// spsp_prodSp_serial
sp_mat spsp_prodSp_serial(const sp_mat& A, const sp_mat& B, const uword& isProdSym, double p);
RcppExport SEXP _spmatMult_spsp_prodSp_serial(SEXP ASEXP, SEXP BSEXP, SEXP isProdSymSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const sp_mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const sp_mat& >::type B(BSEXP);
    Rcpp::traits::input_parameter< const uword& >::type isProdSym(isProdSymSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(spsp_prodSp_serial(A, B, isProdSym, p));
    return rcpp_result_gen;
END_RCPP
}
// spsp_prodSp_openmp
sp_mat spsp_prodSp_openmp(const sp_mat& A, const sp_mat& B, const uword& isProdSym, double p, int nthreads);
RcppExport SEXP _spmatMult_spsp_prodSp_openmp(SEXP ASEXP, SEXP BSEXP, SEXP isProdSymSEXP, SEXP pSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const sp_mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const sp_mat& >::type B(BSEXP);
    Rcpp::traits::input_parameter< const uword& >::type isProdSym(isProdSymSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(spsp_prodSp_openmp(A, B, isProdSym, p, nthreads));
    return rcpp_result_gen;
END_RCPP
}
// spsp_prodSp_fast_openmp
sp_mat spsp_prodSp_fast_openmp(const sp_mat& A, const sp_mat& B, const uword& isProdSym, int nthreads);
RcppExport SEXP _spmatMult_spsp_prodSp_fast_openmp(SEXP ASEXP, SEXP BSEXP, SEXP isProdSymSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const sp_mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const sp_mat& >::type B(BSEXP);
    Rcpp::traits::input_parameter< const uword& >::type isProdSym(isProdSymSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(spsp_prodSp_fast_openmp(A, B, isProdSym, nthreads));
    return rcpp_result_gen;
END_RCPP
}
// spsp_prodMat_openmp
mat spsp_prodMat_openmp(const sp_mat& A, const sp_mat& B, const uword& isProdSym, int nthreads);
RcppExport SEXP _spmatMult_spsp_prodMat_openmp(SEXP ASEXP, SEXP BSEXP, SEXP isProdSymSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const sp_mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const sp_mat& >::type B(BSEXP);
    Rcpp::traits::input_parameter< const uword& >::type isProdSym(isProdSymSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(spsp_prodMat_openmp(A, B, isProdSym, nthreads));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_spmatMult_rcpparma_hello_world", (DL_FUNC) &_spmatMult_rcpparma_hello_world, 0},
    {"_spmatMult_rcpparma_outerproduct", (DL_FUNC) &_spmatMult_rcpparma_outerproduct, 1},
    {"_spmatMult_rcpparma_innerproduct", (DL_FUNC) &_spmatMult_rcpparma_innerproduct, 1},
    {"_spmatMult_rcpparma_bothproducts", (DL_FUNC) &_spmatMult_rcpparma_bothproducts, 1},
    {"_spmatMult_spsp_prodSp_serial", (DL_FUNC) &_spmatMult_spsp_prodSp_serial, 4},
    {"_spmatMult_spsp_prodSp_openmp", (DL_FUNC) &_spmatMult_spsp_prodSp_openmp, 5},
    {"_spmatMult_spsp_prodSp_fast_openmp", (DL_FUNC) &_spmatMult_spsp_prodSp_fast_openmp, 4},
    {"_spmatMult_spsp_prodMat_openmp", (DL_FUNC) &_spmatMult_spsp_prodMat_openmp, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_spmatMult(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
