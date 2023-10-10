#include "sinkhorn.h"
#include "nnls.h"


Rcpp::List reverse_sinkhorn_c(const arma::mat& result_H_row,
                              const arma::mat& result_W_col,
                              const arma::mat& D_vs_row,
                              const arma::mat& D_vs_col,
                              int iterations) {
    arma::mat H_row = result_H_row;
    arma::mat W_col = result_W_col;
    arma::mat H_col;
    arma::mat W_row;

    arma::mat D_ws_col(W_col.n_cols, iterations);  // should be K*iterations
    arma::mat D_hs_row(H_row.n_rows, iterations);  // should be K*iterations

    // arma::mat D_ws_col(V.n_cols, iterations);

    arma::mat ones_like_W(W_col.n_rows, 1, arma::fill::ones);  // M*1
    arma::mat ones_like_H(H_row.n_cols, 1, arma::fill::ones);  // N*1

    arma::vec D_w_inv_pred_vec(W_col.n_cols, arma::fill::zeros);
    arma::vec D_h_inv_pred_vec(H_row.n_rows, arma::fill::zeros);

    D_w_inv_pred_vec = nnls_C__(W_col, ones_like_W);
    W_row = W_col * diagmat(D_w_inv_pred_vec);

    for (int i = iterations - 1; i >= 0; i--) {
        D_h_inv_pred_vec = nnls_C__(H_row.t(), ones_like_H);

        D_ws_col.col(i) = 1 / D_w_inv_pred_vec;
        D_hs_row.col(i) = 1 / D_h_inv_pred_vec;
        H_col = arma::diagmat(D_h_inv_pred_vec) * H_row;
        W_col = arma::diagmat(1 / D_vs_row.col(i)) * W_row * arma::diagmat(1 / D_h_inv_pred_vec);

        if (i != 0) {
            D_w_inv_pred_vec = nnls_C__(W_col, ones_like_W);
            W_row = W_col * diagmat(D_w_inv_pred_vec);
            H_row = diagmat(1 / D_w_inv_pred_vec) * H_col * arma::diagmat(1 / D_vs_col.col(i - 1));
        }
    }

    return Rcpp::List::create(Rcpp::Named("W") = W_col,
                              Rcpp::Named("H") = H_col,
                              Rcpp::Named("D_ws_col") = D_ws_col,
                              Rcpp::Named("D_hs_row") = D_hs_row);
}

Rcpp::List sinkhorn_scale_c(const arma::mat& V, int iterations) {
    arma::mat D_vs_row(V.n_rows, iterations);
    arma::mat D_vs_col(V.n_cols, iterations);

    arma::mat V_row = V;
    arma::mat V_column = V;
    for (int i = 0; i < iterations; i++) {
        arma::mat D_v1 = 1 / arma::sum(V_column, 1);
        D_vs_row.col(i) = D_v1;
        V_row = arma::diagmat(D_v1) * V_column;
        arma::mat D_v2 = 1 / arma::sum(V_row, 0).t();
        D_vs_col.col(i) = D_v2;
        V_column = V_row * arma::diagmat(D_v2);
    }
    return Rcpp::List::create(Rcpp::Named("V_row") = V_row,
                              Rcpp::Named("V_column") = V_column,
                              Rcpp::Named("D_vs_row") = D_vs_row,
                              Rcpp::Named("D_vs_col") = D_vs_col);
}

