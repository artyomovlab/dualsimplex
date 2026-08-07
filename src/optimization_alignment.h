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





//' Optimizing simplex alignment and positivitis.
//  Gradient steps are performed in the transformed space not Sinkhorned space.
//' 
//' @param initial_X current X
//' @param initial_Omega current Omega
//' @param initial_D_w current D_w
//' @param SVRt current SVRt (sigma_ss)
//' @param R current R
//' @param S current S
//' @param coef_der_X learning rate X
//' @param coef_hinge_W beta
//' @param coef_hinge_H lambda
//' @param coef_alignment gamma
//' @param k number of components
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
    const arma::mat& initial_X,
    const arma::mat& initial_Omega,
    const arma::mat& initial_D_w,
    const arma::mat& SVRt,
    const arma::mat& R,
    const arma::mat& S,
    const double coef_der_X,
    double coef_hinge_W,
    double coef_hinge_H,
    double coef_alignment,
    const int k,
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



//' Optimizing simplex alignment and positivitis. This algorithm enforce simplex alignment property.
//  Gradient steps are performed in the transformed space not Sinkhorned space.
//'
//' @param initial_X current X
//' @param initial_Omega current Omega
//' @param initial_D_w current D_w
//' @param SVRt current SVRt (sigma_ss)
//' @param R current R
//' @param S current S
//' @param coef_der_X learning rate X
//' @param coef_hinge_W beta
//' @param coef_hinge_H lambda
//' @param k number of components
//' @param N current N
//' @param M current M
//' @param iterations number of iterations
//' @param convergence_tol tolerance for convergence.
//' @param stop_criteria_window how long error should be on plateu to decrease the learning rate
//' @param debug_stats wether to save grad norm values.
//' @return new parameters
// [[Rcpp::export]]
Rcpp::List optimize_alignment_pgd(
    const arma::mat& initial_X,
    const arma::mat& initial_Omega,
    const arma::mat& initial_D_w,
    const arma::mat& SVRt,
    const arma::mat& R,
    const arma::mat& S,
    const double coef_der_X,
    double coef_hinge_W,
    double coef_hinge_H,
    const int k,
    const double N,
    const double M,
    const int iterations,
    const double convergence_tol,
    const int stop_criteria_window,
    const bool debug_stats
);
