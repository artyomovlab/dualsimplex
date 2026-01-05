// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>


//' Experimental jump norm calculation. legacy
//'
//' @param X matrix to check
//' @param r_const_X constraint
//' @return matrix
// [[Rcpp::export]]
arma::mat jump_norm(arma::mat& X, const double r_const_X = 0);


//' Experimental: Update X in the direction of new X with predefined threshold
//'
//' @param prev_X previous matrix
//' @param new_X new target matrix
//' @param thresh threshold value
//' @return vector of indicies which passed threshold
// [[Rcpp::export]]
arma::uvec update_idx(const arma::mat& prev_X, const arma::mat& new_X, const double thresh = 0.8);

//' Hinge loss derivative for proportion matrix (H)
//'
//' @param H result H matrix obtained from X
//' @param R projection vectors R
//' @param precision_ precision to calculate value
//' @return derivative for X
// [[Rcpp::export]]
arma::mat hinge_der_proportions_C__(const arma::mat& H,
                                    const arma::mat& R,
                                    double precision_ = 1e-10);

//' Hinge loss derivative for basis matrix (W)
//'
//' @param W result W matrix obtained from Omega
//' @param S projection vectors S
//' @param precision_ precision to calculate value
//' @return derivative for Omega
// [[Rcpp::export]]
arma::mat hinge_der_basis_C__(const arma::mat& W, const arma::mat& S, double precision_ = 1e-10);

//' Hinge function value for input matrix X
//'
//' @param X input matrix
//' @return hinge function value
// [[Rcpp::export]]
double hinge_C__(const arma::mat& X);


//' Squared hinge function value for input matrix X
//'
//' @param X input matrix
//' @return hinge function value
// [[Rcpp::export]]
double squared_hinge_C__(const arma::mat& X);


//' Main function to calculate error terms
//'
//' @param X current X
//' @param Omega current Omega
//' @param D_w current D_w
//' @param D_h current D_h
//' @param SVRt current SVRt (Sigma matrix)
//' @param R current R
//' @param S current S
//' @param coef_  # this argument is not used
//' @param coef_der_X learning rate for X
//' @param coef_der_Omega learning rate for Omega
//' @param coef_hinge_H lambda
//' @param coef_hinge_W beta
//' @param coef_pos_D_h experimental coefficient for D_h. legacy. not tested.
//' @param coef_pos_D_w experimental coefficient for D_w.legacy. not tested.
//' @return list with error values
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
//' @return new parameters
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
