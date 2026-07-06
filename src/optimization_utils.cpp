#include "optimization_utils.h"
#include "matrix_utils.h"
#include "nnls.h"

arma::mat jump_norm(arma::mat& X, const double r_const_X) {
    arma::mat norm_(X.n_rows, X.n_cols);
    norm_.fill(1.0);
    arma::mat X_trunc(X.n_rows, X.n_cols - 1);
    arma::uvec ids = arma::regspace<arma::uvec>(1, X.n_cols - 1);
    X_trunc = X.cols(ids);
    for (unsigned int k = 0; k < X.n_rows; k++) {
        double row_norm = norm(X_trunc.row(k), 2);
        for (unsigned int j = 1; j < X.n_cols; j++) {
            if (r_const_X > row_norm) {
                norm_.at(k, j) = r_const_X / row_norm;
            } else {
                norm_.at(k, j) = 1;
            }
        }
    }
    return norm_;
}

arma::uvec update_idx(const arma::mat& prev_X, const arma::mat& new_X, const double thresh) {
    arma::rowvec prev_values = cosine_between_rows(prev_X);
    arma::rowvec new_values = cosine_between_rows(new_X);
    arma::uvec idx2 = find(new_values >= thresh);
    arma::uvec new_idx = {};
    for (unsigned int i = 0; i < idx2.n_elem; i++) {
        if (new_values.at(idx2(i)) >= prev_values.at(idx2(i))) {
            int sz = new_idx.size();
            new_idx.resize(sz + 1);
            new_idx(sz) = idx2(i);
        }
    }
    return new_idx;
}


arma::mat squared_hinge_der_proportions_C__(const arma::mat& H,
                                    const arma::mat& R) {
    int k = H.n_rows;
    arma::mat H_neg = -2 * H;
    H_neg.elem(arma::find(H_neg < 0)).fill(0);
    arma::mat res(k, k, arma::fill::zeros);
    res = H_neg * R.t();
    return res;
}

arma::mat l1_hinge_der_proportions_C__(const arma::mat& H, const arma::mat& R) {
    int k = H.n_rows;
    arma::mat H_neg = - H;
    // H_neg.elem(arma::find(H_neg == 0)).fill(-0.1); we could add this to always have derivative
    H_neg.elem(arma::find(H_neg < 0)).fill(0);
    H_neg.elem(arma::find(H_neg > 0)).fill(-1);

    arma::mat res(k, k, arma::fill::zeros);

    res = H_neg * R.t();
    return res;
}

arma::mat l1_hinge_der_basis_C__(const arma::mat& W, const arma::mat& S) {
    // derivative should be the same as for X but W and Omega are transposed
    arma::mat res = l1_hinge_der_proportions_C__(W.t(), S);
    return res.t();
}

arma::mat squared_hinge_der_basis_C__(const arma::mat& W, const arma::mat& S) {
    // derivative should be the same as for X but W is transposed
    arma::mat res = squared_hinge_der_proportions_C__(W.t(), S);
    return res.t();
}

double hinge_C__(const arma::mat& X) {
    arma::mat X_ = -X;
    X_.elem(arma::find(X_ < 0)).fill(0);
    return accu(X_);
}

double squared_hinge_C__(const arma::mat& X) {
    arma::mat X_ = -X;
    X_.elem(arma::find(X_ < 0)).fill(0);
    X_ %= X_;
    return accu(X_);
}

Rcpp::List calcErrors(const arma::mat& X,
                      const arma::mat& Omega,
                      const arma::mat& D_w,
                      const arma::mat& D_h,
                      const arma::mat& SVRt,
                      const arma::mat& R,
                      const arma::mat& S,
                      const double coef_hinge_H,
                      const double coef_hinge_W,
                      const double coef_pos_D_h,
                      const double coef_pos_D_w) {
    arma::mat D_w_diag = diagmat(D_w);
    double deconv_error = pow(norm(SVRt - Omega * D_w_diag * X, "fro"), 2.0);
    // don't calculate since it is time consuming, should deliver the same minimum as th new one
    // double orig_deconv_error = pow(norm(V_row - S.t() * Omega * D_w_diag * X * R, "fro"), 2);
    double lambda_error = coef_hinge_H * hinge_C__(X * R); // average negative element value per matrix
    double scaled_lambda_error = coef_hinge_H * hinge_C__(X * R) / R.n_cols; // average negative element value per matrix

    double beta_error = coef_hinge_W * hinge_C__(S.t() * Omega); // average negative element value per matrix
    double scaled_beta_error = coef_hinge_W * hinge_C__(S.t() * Omega) / S.n_cols; // average negative element value per matrix

    arma::mat A = arma::sum(R, 1);
    arma::mat B = arma::sum(S, 1);
    double D_h_error = coef_pos_D_h * pow(norm(X.t() * D_h - A, "fro"), 2);
    double D_w_error = coef_pos_D_w * pow(norm(Omega * D_w - B, "fro"), 2);
    double total_error = deconv_error + lambda_error  + beta_error  + D_h_error + D_w_error;
    double scaled_total_error = deconv_error + scaled_lambda_error + scaled_beta_error + D_h_error + D_w_error;
    double average_norm_X = arma::mean(arma::vecnorm(X, 2, 1));
    double average_norm_Omega = arma::sum(arma::vecnorm(Omega, 2, 0));
    double norm_term = (average_norm_X + average_norm_Omega);

    return Rcpp::List::create(Rcpp::Named("deconv_error") = deconv_error,
                              Rcpp::Named("lambda_error") = lambda_error,
                              Rcpp::Named("beta_error") = beta_error,
                              Rcpp::Named("scaled_lambda_error") = scaled_lambda_error,
                              Rcpp::Named("scaled_beta_error") = scaled_beta_error,
                              Rcpp::Named("D_h_error") = D_h_error,
                              Rcpp::Named("D_w_error") = D_w_error,
                              Rcpp::Named("total_error") = total_error,
                              Rcpp::Named("scaled_total_error") = scaled_total_error,
                              Rcpp::Named("average_norm") = norm_term
                            );
}