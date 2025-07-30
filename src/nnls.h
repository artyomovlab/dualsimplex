// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

//' Get nnls solution
//'
//' @param A matrix A
//' @param b coefficients
//' @param max_iter max number of iterations
//' @param tol precision
//' @return derivative for X
//' @export
// [[Rcpp::export]]
arma::mat nnls_C__(arma::mat A, arma::mat b, int max_iter = 500, double tol = 1e-6);


//' Get nnls solution
//'
//' @param A matrix A
//' @param b coefficients
//' @param max_iter max number of iterations
//' @param tol precision
//' @return derivative for X
//' @export
// [[Rcpp::export]]
arma::mat nnls_nonzero_C__(arma::mat A, arma::mat b, int max_iter = 500, double tol = 1e-6);
