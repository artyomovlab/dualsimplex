#pragma once

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//' Main training loop with all gradient steps. This is naive algorithm from version 1.0.
//' It optimizes 3 terms (5 in experimental version). 
//' Gradient steps performed in X and Omega spaces separatelly with NNLS used to merge them together into deconvolution term.
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
//' @param convergence_tol tolerance for convergence.
//' @param stop_criteria_window how long error should be on plateu to decrease the learning rate
//' @param debug_stats wether to save grad norm values.
//' @return new parameters
// [[Rcpp::export]]
Rcpp::List optimize_coordinate_descent(const arma::mat& X,
                             const arma::mat& Omega,
                             const arma::mat& D_w,
                             const arma::mat& SVRt,
                             const arma::mat& R,
                             const arma::mat& S,
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
                             const double convergence_tol=1e-12,
                             const int stop_criteria_window = 1e+5, 
                             const bool debug_stats=false);
