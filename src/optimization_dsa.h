#pragma once

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(spdl)]]
#include <spdl.h>

#include <string>
#include <vector>
#include <memory>
#include <map>
#include <fstream>

#include "optimization_error_metrics.h"
#include "optimization_logger.h"





//' Main training loop with all gradient steps. This is the main algorithm for now. 
//' It optimizes positivity of W and regularize the size of X by dual simplex alignment. 
//' Gradient steps performed in X space.
//'
//' @param X current X
//' @param Omega current Omega
//' @param D_w current D_w
//' @param SVRt current SVRt (sigma)
//' @param R current R
//' @param S current S
//' @param coef_der_X learning rate X
//' @param coef_hinge_W beta
//' @param coef_hinge_H lambda
//' @param coef_alignment gamma?
//' @param cell_types number of components (K)
//' @param N current N
//' @param M current M
//' @param iterations number of iterations
//' @param total_regularization_weight total weight for the regularization terms
//' @param reg_X proportion / regularization coefficient for X
//' @param reg_Omega proportion / regularization coefficient for Omega.
//' @param convergence_tol tolerance for convergence.
//' @param stop_criteria_window how long error should be on plateu to decrease the learning rate
//' @param debug_stats wether to save grad norm values.
//' @return new parameters
// [[Rcpp::export]]
Rcpp::List optimize_alignment(
    const arma::mat& X,
    const arma::mat& Omega,
    const arma::mat& D_w,
    const arma::mat& SVRt,
    const arma::mat& R,
    const arma::mat& S,
    const double coef_der_X,
    double coef_hinge_W,
    double coef_hinge_H,
    double coef_alignment,
    const int cell_types,
    const double N,
    const double M,
    const int iterations,
    double total_regularization_weight,
    const double reg_X,
    const double reg_Omega,
    const double convergence_tol,
    const int stop_criteria_window,
    const bool debug_stats
);



//' Main training loop with all gradient steps.
//' This formulation eliminate the Frobenius norm in the objective function and hardcode the simplex alignment by 
//  how omega_dtilda is found. It optimizes positivity of W and possibly H as the positivity optimization. 
//' Gradient steps performed in X space.
//'
//' @param X current X
//' @param Omega current Omega
//' @param D_w current D_w
//' @param SVRt current SVRt (sigma)
//' @param R current R
//' @param S current S
//' @param coef_der_X learning rate X
//' @param coef_hinge_W beta
//' @param coef_hinge_H lambda
//' @param cell_types number of components (K)
//' @param N current N
//' @param M current M
//' @param iterations number of iterations
//' @param total_regularization_weight total weight for the regularization terms
//' @param reg_X proportion / regularization coefficient for X
//' @param reg_Omega proportion / regularization coefficient for Omega.
//' @param convergence_tol tolerance for convergence.
//' @param stop_criteria_window how long error should be on plateu to decrease the learning rate
//' @param debug_stats wether to save grad norm values.
//' @return new parameters
Rcpp::List optimize_alignment_hard(
    const arma::mat& X,
    const arma::mat& Omega,
    const arma::mat& D_w,
    const arma::mat& SVRt,
    const arma::mat& R,
    const arma::mat& S,
    const double coef_der_X,
    double coef_hinge_W,
    double coef_hinge_H,
    const int cell_types,
    const double N,
    const double M,
    const int iterations,
    double total_regularization_weight,
    const double reg_X,
    const double reg_Omega,
    const double convergence_tol,
    const int stop_criteria_window,
    const bool debug_stats
);