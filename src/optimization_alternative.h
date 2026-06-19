// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>


std::tuple<arma::mat, arma::mat, arma::mat> ensure_D_integrity_c(const arma::mat& X_dtilde,
                              const arma::mat& Omega_dtilde,
                              const arma::vec sqrt_Sigma,
                              const double N,
                              const double M);

//' Transform X and Omega points enforcing the desired equality for first coordinates
//' This is done by moving magnitude from  i-th point of X to respective i-th point of the Omega and vice versa.
//'
//' @param X_dtilde current X_tilde_tilde matrix
//' @param Omega_dtilde current Omega_tilde_tilde matrix
//' @param sqrt_Sigma current sqrt of Omega
//' @param N current sqrt of Omega
//' @param M current sqrt of Omega
//' @return corrected params
//' @export
// [[Rcpp::export]]
Rcpp::List ensure_D_integrity(const arma::mat& X_dtilde,
                              const arma::mat& Omega_dtilde,
                              const arma::vec sqrt_Sigma,
                              const double N,
                              const double M);





//' Main function to calculate error terms
//'
//' @param X current X
//' @param Omega current Omega
//' @param D_w current D_w
//' @param SVRt current SVRt (sigma)
//' @param R current R
//' @param S current S
//' @param coef_der_X learning rate X
//' @param coef_hinge_H lambda
//' @param coef_hinge_W beta
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
Rcpp::List alternative_derivative_stage2(const arma::mat& X,
                             const arma::mat& Omega,
                             const arma::mat& D_w,
                             const arma::mat& SVRt,
                             const arma::mat& R,
                             const arma::mat& S,
                             const double coef_der_X,
                             double coef_hinge_H,
                             double coef_hinge_W,
                             const int cell_types,
                             const double N,
                             const double M,
                             const int iterations,
                             double total_regularization_weight = 0,
                             const double reg_X = 1,
                             const double reg_Omega=1,
                             const double convergence_tol=1e-12,
                             const int stop_criteria_window = 1e+5, 
                             const bool debug_stats=false);