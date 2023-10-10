// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::export]]
arma::mat nnls_C__(arma::mat A, arma::mat b, int max_iter = 500, double tol = 1e-6);
