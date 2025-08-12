// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//' Main function to calculate error terms
//'
//' @param X current X
//' @param Omega current Omega
//' @param D_w current D_w
//' @param SVRt current SVRt (sigma)
//' @param R current R
//' @param S current S
//' @param coef_der_X learning rate X
//' @param coef_der_Omega learning rate Omega
//' @param coef_hinge_H lambda
//' @param coef_hinge_W beta
//' @param coef_pos_D_h experimental coefficient for D. legacy not tested.
//' @param coef_pos_D_w experimental coefficient for D. legacy not tested.
//' @param cell_types number of components (K)
//' @param N current N
//' @param M current M
//' @param iterations number of iterations
//' @param mean_radius_X data dependent restriction for updates
//' @param mean_radius_Omega dependent restriction for updates
//' @param r_const_X experimental. not tested
//' @param r_const_Omega experimental. not tested
//' @param thresh experimental. not tested
//' @param X_center optimization restriction directions for X.
//' @param center_constraint value to restrict optimization from going to far
//' @return new parameters
// [[Rcpp::export]]
Rcpp::List markers_derivative_stage2(const arma::mat& X,
                             const arma::mat& Omega,
                             const arma::mat& D_w,
                             const arma::mat& SVRt,
                             const arma::mat& R,
                             const arma::mat& S,
                             const arma::mat& X_center,
                             const double coef_der_X,
                             const double coef_der_Omega,
                             const double coef_hinge_H,
                             const double coef_hinge_W,
                             const double coef_pos_D_h,
                             const double coef_pos_D_w,
                             const int cell_types,
                             const double N,
                             const double M,
                             const int iterations,
                             const double mean_radius_X,
                             const double mean_radius_Omega,
                             const double r_const_X = 0,
                             const double r_const_Omega = 0,
                             const double thresh = 0.8,
                             const double center_constraint = 0
);
