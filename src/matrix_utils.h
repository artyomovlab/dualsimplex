// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//' Get cosine distance between columns
//'
//' @param input vector a.
//' @param input vector b.
//' @return double value of cosine distance
// [[Rcpp::export]]
double cosine_distance(const arma::rowvec& a,const arma::rowvec& b);

//' Get cosine distance between columns
//'
//' @param X input matrix
//' @return derivative for X
// [[Rcpp::export]]
arma::rowvec cosine_between_rows(const arma::mat& X);

arma::mat correctByNorm(arma::mat& X);

arma::mat getDoubleProjection(const arma::mat& V, const arma::mat& R, const arma::mat& S);

arma::uword getNegative(arma::mat X);

double getSum(arma::mat X, arma::mat M);


//' Get low rank approximation with SVD method.
//'
//' @param X inpit matrix.
//' @param rank desired approximation rank.
//' @param iterations number of iterations to perform.
//' @param left elements cropped (should be zero).
//' @param right elements cropped (infinity by default)
//' @return named list containing new matrix, frobenious history for negative elements and number of negative elements.
//' @export
// [[Rcpp::export]]
Rcpp::List getNonnegativeLowRankApproximationWithSVD(const arma::mat& X,  
                                                     const int rank,
                                                     const int iterations,
                                                     const double left=0,
                                                     const double right=-1);

//
//Rcpp::List getNonnegativeLowRankApproximationWithTangent(const arma::mat& X,
//                                                     const int rank,
//                                                     const int iterations,
//                                                     const double left=0);
//
//

//' Get low rank approximation with SVD method.
//'
//' @param X inpit matrix.
//' @param rank desired approximation rank.
//' @param p number of randomizations.
//' @param k number of randomized columns.
//' @param iterations number of iterations to perform.
//' @param left elements cropped (should be zero).
//' @param right elements cropped (infinity by default)
//' @return named list containing new matrix, frobenious history for negative elements and number of negative elements.
//' @export
// [[Rcpp::export]]
Rcpp::List getNonnegativeLowRankApproximationWithHMT(const arma::mat& X,
                                                         const int rank,
                                                         const int p,
                                                         const int k,
                                                         const int iterations,
                                                         const double left=0,
                                                         const double right=-1);
//
//Rcpp::List getNonnegativeLowRankApproximationWithTropp(const arma::mat& X,
//                                                     const int rank,
//                                                     const int k,
//                                                     const int l,
//                                                     const double rho,
//                                                     const int iterations,
//                                                     const double left=0);
//

//' Get low rank approximation with GN method.
//'
//' @param X inpit matrix
//' @param rank desired approximation rank
//' @param l parameter for Psi
//' @param iterations number of iterations to perform
//' @param left elements cropped (should be zero)
//' @param right elements cropped (infinity by default)
//' @return named list containing new matrix, frobenious history for negative elements and number of negative elements.
//' @export
// [[Rcpp::export]]
Rcpp::List getNonnegativeLowRankApproximationWithGN(const arma::mat& X,
                                                      const int rank,
                                                      const int l,
                                                      const int iterations,
                                                      const double left=0,
                                                      const double right=-1);


//' Get low rank approximation with Tangent method.
//'
//' @param X inpit matrix.
//' @param rank desired approximation rank.
//' @param iterations number of iterations to perform.
//' @param left elements cropped (should be zero).
//' @param right elements cropped (infinity by default)
//' @return named list containing new matrix, frobenious history for negative elements and number of negative elements.
//' @export
// [[Rcpp::export]]
Rcpp::List getNonnegativeLowRankApproximationWithTangentMethod(const arma::mat& X,  
                                                     const int rank,
                                                     const int iterations,
                                                     const double left=0, 
                                                     const double right=-1);