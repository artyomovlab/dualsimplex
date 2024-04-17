// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// alternative_derivative_stage2
Rcpp::List alternative_derivative_stage2(const arma::mat& X, const arma::mat& Omega, const arma::mat& D_w, const arma::mat& SVRt, const arma::mat& R, const arma::mat& S, const double coef_der_X, const double coef_der_Omega, const double coef_hinge_H, const double coef_hinge_W, const double coef_pos_D_h, const double coef_pos_D_w, const int cell_types, const double N, const double M, const int iterations, const double mean_radius_X, const double mean_radius_Omega, const double r_const_X, const double r_const_Omega, const double thresh);
RcppExport SEXP _DualSimplex_alternative_derivative_stage2(SEXP XSEXP, SEXP OmegaSEXP, SEXP D_wSEXP, SEXP SVRtSEXP, SEXP RSEXP, SEXP SSEXP, SEXP coef_der_XSEXP, SEXP coef_der_OmegaSEXP, SEXP coef_hinge_HSEXP, SEXP coef_hinge_WSEXP, SEXP coef_pos_D_hSEXP, SEXP coef_pos_D_wSEXP, SEXP cell_typesSEXP, SEXP NSEXP, SEXP MSEXP, SEXP iterationsSEXP, SEXP mean_radius_XSEXP, SEXP mean_radius_OmegaSEXP, SEXP r_const_XSEXP, SEXP r_const_OmegaSEXP, SEXP threshSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Omega(OmegaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type D_w(D_wSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type SVRt(SVRtSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type R(RSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S(SSEXP);
    Rcpp::traits::input_parameter< const double >::type coef_der_X(coef_der_XSEXP);
    Rcpp::traits::input_parameter< const double >::type coef_der_Omega(coef_der_OmegaSEXP);
    Rcpp::traits::input_parameter< const double >::type coef_hinge_H(coef_hinge_HSEXP);
    Rcpp::traits::input_parameter< const double >::type coef_hinge_W(coef_hinge_WSEXP);
    Rcpp::traits::input_parameter< const double >::type coef_pos_D_h(coef_pos_D_hSEXP);
    Rcpp::traits::input_parameter< const double >::type coef_pos_D_w(coef_pos_D_wSEXP);
    Rcpp::traits::input_parameter< const int >::type cell_types(cell_typesSEXP);
    Rcpp::traits::input_parameter< const double >::type N(NSEXP);
    Rcpp::traits::input_parameter< const double >::type M(MSEXP);
    Rcpp::traits::input_parameter< const int >::type iterations(iterationsSEXP);
    Rcpp::traits::input_parameter< const double >::type mean_radius_X(mean_radius_XSEXP);
    Rcpp::traits::input_parameter< const double >::type mean_radius_Omega(mean_radius_OmegaSEXP);
    Rcpp::traits::input_parameter< const double >::type r_const_X(r_const_XSEXP);
    Rcpp::traits::input_parameter< const double >::type r_const_Omega(r_const_OmegaSEXP);
    Rcpp::traits::input_parameter< const double >::type thresh(threshSEXP);
    rcpp_result_gen = Rcpp::wrap(alternative_derivative_stage2(X, Omega, D_w, SVRt, R, S, coef_der_X, coef_der_Omega, coef_hinge_H, coef_hinge_W, coef_pos_D_h, coef_pos_D_w, cell_types, N, M, iterations, mean_radius_X, mean_radius_Omega, r_const_X, r_const_Omega, thresh));
    return rcpp_result_gen;
END_RCPP
}
// find_cosine
arma::rowvec find_cosine(const arma::mat& X);
RcppExport SEXP _DualSimplex_find_cosine(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(find_cosine(X));
    return rcpp_result_gen;
END_RCPP
}
// nnls_C__
arma::mat nnls_C__(arma::mat A, arma::mat b, int max_iter, double tol);
RcppExport SEXP _DualSimplex_nnls_C__(SEXP ASEXP, SEXP bSEXP, SEXP max_iterSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type b(bSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(nnls_C__(A, b, max_iter, tol));
    return rcpp_result_gen;
END_RCPP
}
// jump_norm
arma::mat jump_norm(arma::mat& X, const double r_const_X);
RcppExport SEXP _DualSimplex_jump_norm(SEXP XSEXP, SEXP r_const_XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const double >::type r_const_X(r_const_XSEXP);
    rcpp_result_gen = Rcpp::wrap(jump_norm(X, r_const_X));
    return rcpp_result_gen;
END_RCPP
}
// update_idx
arma::uvec update_idx(const arma::mat& prev_X, const arma::mat& new_X, const double thresh);
RcppExport SEXP _DualSimplex_update_idx(SEXP prev_XSEXP, SEXP new_XSEXP, SEXP threshSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type prev_X(prev_XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type new_X(new_XSEXP);
    Rcpp::traits::input_parameter< const double >::type thresh(threshSEXP);
    rcpp_result_gen = Rcpp::wrap(update_idx(prev_X, new_X, thresh));
    return rcpp_result_gen;
END_RCPP
}
// hinge_der_proportions_C__
arma::mat hinge_der_proportions_C__(const arma::mat& H, const arma::mat& R, double precision_);
RcppExport SEXP _DualSimplex_hinge_der_proportions_C__(SEXP HSEXP, SEXP RSEXP, SEXP precision_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type H(HSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type R(RSEXP);
    Rcpp::traits::input_parameter< double >::type precision_(precision_SEXP);
    rcpp_result_gen = Rcpp::wrap(hinge_der_proportions_C__(H, R, precision_));
    return rcpp_result_gen;
END_RCPP
}
// hinge_der_basis_C__
arma::mat hinge_der_basis_C__(const arma::mat& W, const arma::mat& S, double precision_);
RcppExport SEXP _DualSimplex_hinge_der_basis_C__(SEXP WSEXP, SEXP SSEXP, SEXP precision_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type W(WSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S(SSEXP);
    Rcpp::traits::input_parameter< double >::type precision_(precision_SEXP);
    rcpp_result_gen = Rcpp::wrap(hinge_der_basis_C__(W, S, precision_));
    return rcpp_result_gen;
END_RCPP
}
// hinge_C__
double hinge_C__(const arma::mat& X);
RcppExport SEXP _DualSimplex_hinge_C__(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(hinge_C__(X));
    return rcpp_result_gen;
END_RCPP
}
// calcErrors
Rcpp::List calcErrors(const arma::mat& X, const arma::mat& Omega, const arma::mat& D_w, const arma::mat& D_h, const arma::mat& SVRt, const arma::mat& R, const arma::mat& S, const double coef_, const double coef_der_X, const double coef_der_Omega, const double coef_hinge_H, const double coef_hinge_W, const double coef_pos_D_h, const double coef_pos_D_w);
RcppExport SEXP _DualSimplex_calcErrors(SEXP XSEXP, SEXP OmegaSEXP, SEXP D_wSEXP, SEXP D_hSEXP, SEXP SVRtSEXP, SEXP RSEXP, SEXP SSEXP, SEXP coef_SEXP, SEXP coef_der_XSEXP, SEXP coef_der_OmegaSEXP, SEXP coef_hinge_HSEXP, SEXP coef_hinge_WSEXP, SEXP coef_pos_D_hSEXP, SEXP coef_pos_D_wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Omega(OmegaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type D_w(D_wSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type D_h(D_hSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type SVRt(SVRtSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type R(RSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S(SSEXP);
    Rcpp::traits::input_parameter< const double >::type coef_(coef_SEXP);
    Rcpp::traits::input_parameter< const double >::type coef_der_X(coef_der_XSEXP);
    Rcpp::traits::input_parameter< const double >::type coef_der_Omega(coef_der_OmegaSEXP);
    Rcpp::traits::input_parameter< const double >::type coef_hinge_H(coef_hinge_HSEXP);
    Rcpp::traits::input_parameter< const double >::type coef_hinge_W(coef_hinge_WSEXP);
    Rcpp::traits::input_parameter< const double >::type coef_pos_D_h(coef_pos_D_hSEXP);
    Rcpp::traits::input_parameter< const double >::type coef_pos_D_w(coef_pos_D_wSEXP);
    rcpp_result_gen = Rcpp::wrap(calcErrors(X, Omega, D_w, D_h, SVRt, R, S, coef_, coef_der_X, coef_der_Omega, coef_hinge_H, coef_hinge_W, coef_pos_D_h, coef_pos_D_w));
    return rcpp_result_gen;
END_RCPP
}
// derivative_stage2
Rcpp::List derivative_stage2(const arma::mat& X, const arma::mat& Omega, const arma::mat& D_w, const arma::mat& SVRt, const arma::mat& R, const arma::mat& S, const double coef_der_X, const double coef_der_Omega, const double coef_hinge_H, const double coef_hinge_W, const double coef_pos_D_h, const double coef_pos_D_w, const int cell_types, const double N, const double M, const int iterations, const double mean_radius_X, const double mean_radius_Omega, const double r_const_X, const double r_const_Omega, const double thresh);
RcppExport SEXP _DualSimplex_derivative_stage2(SEXP XSEXP, SEXP OmegaSEXP, SEXP D_wSEXP, SEXP SVRtSEXP, SEXP RSEXP, SEXP SSEXP, SEXP coef_der_XSEXP, SEXP coef_der_OmegaSEXP, SEXP coef_hinge_HSEXP, SEXP coef_hinge_WSEXP, SEXP coef_pos_D_hSEXP, SEXP coef_pos_D_wSEXP, SEXP cell_typesSEXP, SEXP NSEXP, SEXP MSEXP, SEXP iterationsSEXP, SEXP mean_radius_XSEXP, SEXP mean_radius_OmegaSEXP, SEXP r_const_XSEXP, SEXP r_const_OmegaSEXP, SEXP threshSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Omega(OmegaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type D_w(D_wSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type SVRt(SVRtSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type R(RSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S(SSEXP);
    Rcpp::traits::input_parameter< const double >::type coef_der_X(coef_der_XSEXP);
    Rcpp::traits::input_parameter< const double >::type coef_der_Omega(coef_der_OmegaSEXP);
    Rcpp::traits::input_parameter< const double >::type coef_hinge_H(coef_hinge_HSEXP);
    Rcpp::traits::input_parameter< const double >::type coef_hinge_W(coef_hinge_WSEXP);
    Rcpp::traits::input_parameter< const double >::type coef_pos_D_h(coef_pos_D_hSEXP);
    Rcpp::traits::input_parameter< const double >::type coef_pos_D_w(coef_pos_D_wSEXP);
    Rcpp::traits::input_parameter< const int >::type cell_types(cell_typesSEXP);
    Rcpp::traits::input_parameter< const double >::type N(NSEXP);
    Rcpp::traits::input_parameter< const double >::type M(MSEXP);
    Rcpp::traits::input_parameter< const int >::type iterations(iterationsSEXP);
    Rcpp::traits::input_parameter< const double >::type mean_radius_X(mean_radius_XSEXP);
    Rcpp::traits::input_parameter< const double >::type mean_radius_Omega(mean_radius_OmegaSEXP);
    Rcpp::traits::input_parameter< const double >::type r_const_X(r_const_XSEXP);
    Rcpp::traits::input_parameter< const double >::type r_const_Omega(r_const_OmegaSEXP);
    Rcpp::traits::input_parameter< const double >::type thresh(threshSEXP);
    rcpp_result_gen = Rcpp::wrap(derivative_stage2(X, Omega, D_w, SVRt, R, S, coef_der_X, coef_der_Omega, coef_hinge_H, coef_hinge_W, coef_pos_D_h, coef_pos_D_w, cell_types, N, M, iterations, mean_radius_X, mean_radius_Omega, r_const_X, r_const_Omega, thresh));
    return rcpp_result_gen;
END_RCPP
}
// reverse_sinkhorn_c
Rcpp::List reverse_sinkhorn_c(const arma::mat& result_H_row, const arma::mat& result_W_col, const arma::mat& D_vs_row, const arma::mat& D_vs_col, int iterations);
RcppExport SEXP _DualSimplex_reverse_sinkhorn_c(SEXP result_H_rowSEXP, SEXP result_W_colSEXP, SEXP D_vs_rowSEXP, SEXP D_vs_colSEXP, SEXP iterationsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type result_H_row(result_H_rowSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type result_W_col(result_W_colSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type D_vs_row(D_vs_rowSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type D_vs_col(D_vs_colSEXP);
    Rcpp::traits::input_parameter< int >::type iterations(iterationsSEXP);
    rcpp_result_gen = Rcpp::wrap(reverse_sinkhorn_c(result_H_row, result_W_col, D_vs_row, D_vs_col, iterations));
    return rcpp_result_gen;
END_RCPP
}
// sinkhorn_scale_c
Rcpp::List sinkhorn_scale_c(const arma::mat& V, int iterations);
RcppExport SEXP _DualSimplex_sinkhorn_scale_c(SEXP VSEXP, SEXP iterationsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type V(VSEXP);
    Rcpp::traits::input_parameter< int >::type iterations(iterationsSEXP);
    rcpp_result_gen = Rcpp::wrap(sinkhorn_scale_c(V, iterations));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_DualSimplex_alternative_derivative_stage2", (DL_FUNC) &_DualSimplex_alternative_derivative_stage2, 21},
    {"_DualSimplex_find_cosine", (DL_FUNC) &_DualSimplex_find_cosine, 1},
    {"_DualSimplex_nnls_C__", (DL_FUNC) &_DualSimplex_nnls_C__, 4},
    {"_DualSimplex_jump_norm", (DL_FUNC) &_DualSimplex_jump_norm, 2},
    {"_DualSimplex_update_idx", (DL_FUNC) &_DualSimplex_update_idx, 3},
    {"_DualSimplex_hinge_der_proportions_C__", (DL_FUNC) &_DualSimplex_hinge_der_proportions_C__, 3},
    {"_DualSimplex_hinge_der_basis_C__", (DL_FUNC) &_DualSimplex_hinge_der_basis_C__, 3},
    {"_DualSimplex_hinge_C__", (DL_FUNC) &_DualSimplex_hinge_C__, 1},
    {"_DualSimplex_calcErrors", (DL_FUNC) &_DualSimplex_calcErrors, 14},
    {"_DualSimplex_derivative_stage2", (DL_FUNC) &_DualSimplex_derivative_stage2, 21},
    {"_DualSimplex_reverse_sinkhorn_c", (DL_FUNC) &_DualSimplex_reverse_sinkhorn_c, 5},
    {"_DualSimplex_sinkhorn_scale_c", (DL_FUNC) &_DualSimplex_sinkhorn_scale_c, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_DualSimplex(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
