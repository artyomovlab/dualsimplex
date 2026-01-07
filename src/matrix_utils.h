// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//' Get cosine similarity between columns
//'
//' @param a input vector a.
//' @param b input vector b.
//' @return double value of cosine distance
// [[Rcpp::export]]
double cosine_similarity(const arma::rowvec& a,const arma::rowvec& b);

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

//' Get relative coordinates with respect to vertices.
//' This is calculated as a determinant rations of the respective simplexes.
//'
//' @param projected_points coordinates of rows/columns in svd space (number_of_poitns x K).
//' @param solution_points solution points to calculate respective relative coordinates (K x K).
//' @return arma::mat coordinates of the points with respect to simplex vertices (number_of_points x K).
//' @export
// [[Rcpp::export]]
arma::mat get_relative_coordinates(const arma::mat& projected_points,const arma::mat& solution_points);


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