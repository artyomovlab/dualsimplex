// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

//' Reverse Sinkhorn scaling method
//'
//' @param result_H_row H_ss returned by DualSimplex. (row normalized X*R).
//' @param result_W_col W_gs returned by DualSimplex. (row normalized t(S)*Omega).
//' @param D_vs_row row normalizing matrices used for V in forward procedure.
//' @param D_vs_col column normalizing matrices used for V in forward procedure.
//' @param iterations how many iterations back
//' @return named list of W, H, Dv_inv_W_row, H_row, D_ws_col, D_hs_row.
//' @export
// [[Rcpp::export]]
Rcpp::List reverse_sinkhorn_c(const arma::mat& result_H_row,
                              const arma::mat& result_W_col,
                              const arma::mat& D_vs_row,
                              const arma::mat& D_vs_col,
                              int iterations);

//' Forward Sinkhorn scaling method
//'
//' @param V matrix to scale
//' @param iterations W_gs number of iterations.
//' @return named list of V_row, V_col, Dvs_row, Dvs_col
//' @export
// [[Rcpp::export]]
Rcpp::List sinkhorn_scale_c(const arma::mat& V, int iterations);



//' More efficient forward Sinkhorn scaling algorithm which check the convergence
//'
//' @param V matrix to scale.
//' @param max_iter Maximum iterations of the Sinkhorn scaling. Default is 20 iterations.
//' @param iter_start_check From which iteration should the function checks the convergence. By default, check started from iteration 5.
//' @param check_every_iter How offeten should we check the convergence. The default is check every 3 iterations.
//' @param epsilon The tolerance for convergece. Default value is 1.490116e-08, which is similar to R's built in `all.equal` function.
//' @return named list of V_row, V_col, D_row, D_col, iterations.
//' @export
// [[Rcpp::export]]
Rcpp::List efficient_sinkhorn(
    const arma::mat& V,
    const int max_iter = 20,
    const int iter_start_check = 5,
    const int check_every_iter = 3,
    const double epsilon = 1.490116e-08 // similar to R's all.equal
);

//' Extended version of sinkhorn transformation, returning all Dv, Dh, Dw matrices
//'
//' @param V matrix to scale. (V=WH)
//' @param W matrix to scale.
//' @param H matrix to scale.
//' @param n_iter exact number of iterations
//' @export
// [[Rcpp::export]]
Rcpp::List extended_sinkhorn(const arma::mat& V,
                              const arma::mat& W,
                              const arma::mat& H,
                              const int n_iter);


//' Helper function to perfrom iteration scaling
//'
//' @param V matrix to scale.
//' @param D_vs_row matrix stroing row scaling at each iteration  
//' @param D_vs_col matrix stroing column scaling at each iteration
//' @param iter iteration at which Sinkhorn is converged
//' @param return_col_norm whether to perform very last normalization
//' @return scaled matrix
//' @export
// [[Rcpp::export]]
arma::mat sinkhorn_sweep_c(
    const arma::mat& V,
    const arma::mat& D_vs_row,
    const arma::mat& D_vs_col,
    unsigned int iter,
    unsigned int return_col_norm
);


