#include "optimization.h"

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
        if (new_values.at(idx2[i]) >= prev_values.at(idx2[i])) {
            int sz = new_idx.size();
            new_idx.resize(sz + 1);
            new_idx(sz) = idx2[i];
        }
    }
    return new_idx;
}

arma::mat hinge_der_proportions_C__(const arma::mat& H, const arma::mat& R, double precision_) {
    int m = H.n_rows;
    int n = H.n_cols;

    arma::mat TMP(n, m * m, arma::fill::zeros);

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            if (H(i, j) < 0) {
                for (int k = 1; k < m; k++) {
                    TMP(j, k + i * m) = -R(k, j);
                }
            }
        }
    }

    return reshape(arma::sum(TMP, 0), m, m).t();
}

arma::mat hinge_der_basis_C__(const arma::mat& W, const arma::mat& S, double precision_) {
    int n = W.n_cols;

    arma::mat res(n, n, arma::fill::zeros);

    for (int j = 0; j < n; j++) {
        arma::vec t = W.col(j);
        res.col(j) = arma::sum(-S.cols(find(t < -precision_)), 1);
    }
    res.row(0).zeros();

    return res;
}

double hinge_C__(const arma::mat& X) {
    arma::mat X_(X.n_rows, X.n_cols, arma::fill::zeros);
    double elem_ = 0;

    for (unsigned int i = 0; i < X.n_rows; i++) {
        for (unsigned int j = 0; j < X.n_cols; j++) {
            elem_ = X(i, j);
            if (elem_ < 0) {
                X_(i, j) = -elem_;
            }
        }
    }
    return accu(X_);
}

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
                      const double coef_pos_D_w) {
    arma::mat D_w_diag = diagmat(D_w);

    double deconv_error = pow(norm(SVRt - Omega * D_w_diag * X, "fro"), 2.0);
    // don't calculate since it is time consuming, should deliver the same minimum as th new one
    // double orig_deconv_error = pow(norm(V_row - S.t() * Omega * D_w_diag * X * R, "fro"), 2);
    double lambda_error = coef_ * coef_hinge_H * hinge_C__(X * R);
    double beta_error = coef_ * coef_hinge_W * hinge_C__(S.t() * Omega);
    arma::mat A = arma::sum(R, 1);
    arma::mat B = arma::sum(S, 1);
    double D_h_error = coef_pos_D_h * pow(norm(X.t() * D_h - A, "fro"), 2);
    double D_w_error = coef_pos_D_w * pow(norm(Omega * D_w - B, "fro"), 2);
    double new_error = deconv_error + lambda_error + beta_error + D_h_error + D_w_error;

    return Rcpp::List::create(Rcpp::Named("deconv_error") = deconv_error,
                              Rcpp::Named("lambda_error") = lambda_error,
                              Rcpp::Named("beta_error") = beta_error,
                              Rcpp::Named("D_h_error") = D_h_error,
                              Rcpp::Named("D_w_error") = D_w_error,
                              Rcpp::Named("total_error") = new_error);
}

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
                             const double r_const_X,
                             const double r_const_Omega,
                             const double thresh) {
    arma::mat errors_statistics(iterations, 9, arma::fill::zeros);
    arma::mat points_statistics_X(iterations, cell_types * cell_types, arma::fill::zeros);
    arma::mat points_statistics_Omega(iterations, cell_types * cell_types, arma::fill::zeros);

    arma::mat new_X = X;
    arma::mat new_Omega = Omega;
    arma::mat new_D_w = D_w;
    arma::mat new_D_h = new_D_w * (N / M);
    arma::mat jump_X, jump_Omega;

    arma::vec vectorised_SVRt = arma::vectorise(SVRt);
    arma::colvec sum_rows_R = arma::sum(R, 1);
    arma::colvec sum_rows_S = arma::sum(S, 1);

    arma::mat B = arma::join_cols(vectorised_SVRt, coef_pos_D_w * sum_rows_S);
    arma::mat C = arma::join_cols(vectorised_SVRt, coef_pos_D_h * sum_rows_R);
    arma::mat der_X, der_Omega;

    for (int itr_ = 0; itr_ < iterations; itr_++) {
        // derivative X
        der_X =
            -2 * (diagmat(new_D_w) * new_Omega.t() * (SVRt - new_Omega * diagmat(new_D_w) * new_X));
        der_X += coef_hinge_H * hinge_der_proportions_C__(new_X * R, R);
        der_X += coef_pos_D_h * 2 * new_D_h * (new_X.t() * new_D_h - sum_rows_R).t();
        der_X.col(0).zeros();
        der_X = correctByNorm(der_X) * mean_radius_X;

        // cosine threshold correction if needed
        if (thresh > 0) {
            arma::mat tmp_X = (new_X - coef_der_X * der_X).t();
            arma::mat tmp_X_2 = (new_X).t();
            arma::uvec idx = update_idx(tmp_X, tmp_X_2, thresh);

            if (idx.n_elem > 0) {
                der_X.rows(idx).zeros();
            }
        }
        // Update X
        new_X = new_X - coef_der_X * der_X;
        // threshold for length of the new X
        if (r_const_X > 0) {
            jump_X = jump_norm(new_X, r_const_X);
            new_X = new_X % jump_X;
        }

        arma::mat vec_mtx(cell_types * cell_types, cell_types, arma::fill::zeros);
        for (int c = 0; c < cell_types; c++) {
            vec_mtx.col(c) = arma::vectorise(new_Omega.col(c) * new_X.row(c));
        }
        arma::mat A = arma::join_cols((M / N) * vec_mtx, coef_pos_D_h * new_X.t());

        new_D_h = nnls_C__(A, C);
        new_D_w = new_D_h * (M / N);

        // derivative Omega
        der_Omega =
            -2 * (SVRt - new_Omega * diagmat(new_D_w) * new_X) * new_X.t() * diagmat(new_D_w);
        der_Omega += coef_hinge_W * hinge_der_basis_C__(S.t() * new_Omega, S);
        der_Omega += coef_pos_D_w * 2 * (new_Omega * new_D_w - sum_rows_S) * new_D_w.t();
        der_Omega.row(0).zeros();
        der_Omega = correctByNorm(der_Omega) * mean_radius_Omega;

        if (thresh > 0) {
            arma::mat tmp_Omega = new_Omega - coef_der_Omega * der_Omega;
            arma::uvec idx2 = update_idx(tmp_Omega, new_Omega, thresh);

            if (idx2.n_elem > 0) {
                der_Omega.cols(idx2).zeros();
            }
        }

        new_Omega = new_Omega - coef_der_Omega * der_Omega;

        if (r_const_Omega > 0) {
            arma::mat t_Omega = new_Omega.t();
            jump_Omega = jump_norm(t_Omega, r_const_Omega);
            jump_Omega = jump_Omega.t();
            // has_jump_Omega = any(jump_Omega != 1);
            new_Omega = new_Omega % jump_Omega;
        }

        vec_mtx.fill(arma::fill::zeros);
        A.fill(arma::fill::zeros);

        for (int c = 0; c < cell_types; c++) {
            vec_mtx.col(c) = arma::vectorise(new_Omega.col(c) * new_X.row(c));
        }
        A = arma::join_cols(vec_mtx, coef_pos_D_w * new_Omega);

        new_D_w = nnls_C__(A, B);

        new_D_h = new_D_w * (N / M);

        arma::uword neg_props = getNegative(new_X * R);
        arma::uword neg_basis = getNegative(S.t() * new_Omega);
        double sum_ = accu(new_D_w) / M;

        Rcpp::List current_errors = calcErrors(new_X,
                                               new_Omega,
                                               new_D_w,
                                               new_D_h,
                                               SVRt,
                                               R,
                                               S,
                                               1,
                                               coef_der_X,
                                               coef_der_Omega,
                                               coef_hinge_H,
                                               coef_hinge_W,
                                               coef_pos_D_h,
                                               coef_pos_D_w);

        errors_statistics.row(itr_) = arma::rowvec{current_errors["deconv_error"],
                                                   current_errors["lambda_error"],
                                                   current_errors["beta_error"],
                                                   current_errors["D_h_error"],
                                                   current_errors["D_w_error"],
                                                   current_errors["total_error"],
                                                   static_cast<double>(neg_props),
                                                   static_cast<double>(neg_basis),
                                                   sum_};
        points_statistics_X.row(itr_) = new_X.as_row();
        points_statistics_Omega.row(itr_) = new_Omega.as_row();
    }

    return Rcpp::List::create(Rcpp::Named("new_X") = new_X,
                              Rcpp::Named("new_Omega") = new_Omega,
                              Rcpp::Named("new_D_w") = new_D_w,
                              Rcpp::Named("new_D_h") = new_D_h,
                              Rcpp::Named("errors_statistics") = errors_statistics,
                              Rcpp::Named("points_statistics_X") = points_statistics_X,
                              Rcpp::Named("points_statistics_Omega") = points_statistics_Omega);
}
