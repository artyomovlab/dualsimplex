// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

//' Reverse Sinkhorn scaling method
//'
//' @param result_H_row H_ss returned by DualSimplex. (row normalized X*R).
//' @param result_W_col W_gs returned by DualSimplex. (row normalized t(S)*Omega).
//' @param D_vs_row row normalizing matrices used for V in forward procedure.
//' @param D_vs_col column normalizing matrices used for V in forward procedure.
//' @param iterations how many iterations back
//' @return named list of W, H, Dv_inv_W_row, H_row, D_ws_col, D_hs_row.
//' @export
// [[Rcpp::export]]
Rcpp::List reverse_sinkhorn_c(const arma::mat& result_H_row,
                              const arma::mat& result_W_col,
                              const arma::mat& D_vs_row,
                              const arma::mat& D_vs_col,
                              int iterations);

//' Forward Sinkhorn scaling method
//'
//' @param V matrix to scale
//' @param iterations W_gs number of iterations.
//' @return named list of V_row, V_col, Dvs_row, Dvs_col
//' @export
// [[Rcpp::export]]
Rcpp::List sinkhorn_scale_c(const arma::mat& V, int iterations);
