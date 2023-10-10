// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::rowvec find_cosine(const arma::mat& X);

arma::mat correctByNorm(arma::mat& X);

arma::mat getDoubleProjection(const arma::mat& V, const arma::mat& R, const arma::mat& S);

arma::uword getNegative(arma::mat X);

double getSum(arma::mat X, arma::mat M);
