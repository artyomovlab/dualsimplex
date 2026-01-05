// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

arma::mat alternative_hinge_der_basis_C__(const arma::mat& W, const arma::mat& S, double precision_ = 1e-10);

//' Hinge loss derivative for proportion matrix (H)
//'
//' @param H result H matrix obtained from X
//' @param R projection vectors R
//' @return derivative for X
// [[Rcpp::export]]
arma::mat squared_hinge_der_proportions_C__(const arma::mat& H,
                                    const arma::mat& R);

//' Hinge loss derivative for proportion matrix (H)
//'
//' @param H result H matrix obtained from X
//' @param R projection vectors R
//' @return derivative for X
// [[Rcpp::export]]
arma::mat l1_hinge_der_proportions_C__(const arma::mat& H,
                                    const arma::mat& R);

//' Hinge loss derivative for basis matrix (W)
//'
//' @param W result W matrix obtained from Omega
//' @param S projection vectors S
//' @return derivative for Omega
// [[Rcpp::export]]
arma::mat squared_hinge_der_basis_C__(const arma::mat& W, const arma::mat& S);


//' Hinge loss derivative for basis matrix (W)
//'
//' @param W result W matrix obtained from Omega
//' @param S projection vectors S
//' @return derivative for Omega
// [[Rcpp::export]]
arma::mat l1_hinge_der_basis_C__(const arma::mat& W, const arma::mat& S);


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
//' @param solution_balancing_threshold experimental. If solution is to far away we re-balance norms of the solution vectors between X and Omega
//' @return new parameters
// [[Rcpp::export]]
Rcpp::List alternative_derivative_stage2(const arma::mat& X,
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
                             const double solution_balancing_threshold = 10000);
