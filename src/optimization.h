// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::mat jump_norm(arma::mat& X, const double r_const_X = 0);

// [[Rcpp::export]]
arma::uvec update_idx(const arma::mat& prev_X, const arma::mat& new_X, const double thresh = 0.8);

// [[Rcpp::export]]
arma::mat hinge_der_proportions_C__(const arma::mat& H,
                                    const arma::mat& R,
                                    double precision_ = 1e-10);

arma::mat hinge_der_basis_C__(const arma::mat& W, const arma::mat& S, double precision_ = 1e-10);

// [[Rcpp::export]]
double hinge_C__(const arma::mat& X);

// [[Rcpp::export]]
Rcpp::List calcErrors(const arma::mat& X,
                      const arma::mat& Omega,
                      const arma::mat& D_w,
                      const arma::mat& D_h,
                      const arma::mat& SVRt,
                      const arma::mat& R,
                      const arma::mat& S,
                      const double coef_,
                      const double coef_der_X,
                      const double coef_der_Omega,
                      const double coef_hinge_H,
                      const double coef_hinge_W,
                      const double coef_pos_D_h,
                      const double coef_pos_D_w);

// [[Rcpp::export]]
Rcpp::List derivative_stage2(const arma::mat& X,
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
                             const double thresh = 0.8);
