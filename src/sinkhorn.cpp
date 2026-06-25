#include "sinkhorn.h"
#include "nnls.h"
#include "matrix_utils.h"

//' Reverse Sinkhorn scaling method without any nnls usage.
//'
//' @param result_H_col H_gs calculated geometrically from the solution
//' @param result_W_row W_ss calculated geometrically from the solution.
//' @param D_vs_row row normalizing matrices used for V in forward procedure.
//' @param D_vs_col column normalizing matrices used for V in forward procedure.
//' @param iterations how many iterations back
//' @return named list of W, H, Dv_inv_W_row, H_row, D_ws_col, D_hs_row.
Rcpp::List geometrical_reverse_sinkhorn_c(
                              const arma::mat& result_H_col,
                              const arma::mat& result_W_row,
                              const arma::mat& D_vs_row,
                              const arma::mat& D_vs_col,
                              arma::mat V_inf_row,    
                              int iterations) {
    Rcpp::Rcout << "Method enter " << ".\n";
    arma::mat H_col = result_H_col;
    arma::mat W_row = result_W_row;
    arma::mat H_row;
    arma::mat W_col;
    arma::mat D_w , D_h; //intermediate matrices to move some magnitude between matrices
    D_h = 1 / arma::sum(H_col, 1);

    // V_fs = W_ss * D_w * D_w_inv * H_ss * D_v_last
    H_row = arma::diagmat(D_h) * H_col;

    arma::mat current_V = V_inf_row;
    arma::mat S_t, R_t;
    arma::vec Sigma;
    arma::mat projected_points;
    arma::mat current_solution;
    int K = H_row.n_rows;
    Rcpp::Rcout << "Loop starting " << ".\n";
    // start iteratively update them
    for (int i = iterations - 2; i >= 0; i--) {
        Rcpp::Rcout << "Iter " << i << ".\n";

        // Transformation for W
        W_row = arma::diagmat(1 / D_vs_row.col(i+1)) * W_row;
        // column normalize  W_row -- is it valid ?
        D_w = 1 / arma::sum(W_row, 0);
        W_col = W_row *  arma::diagmat(D_w);

        // Transformation for V
        current_V = arma::diagmat(1 / D_vs_row.col(i+1)) * current_V;
        // V should be column normalized -- transform it accordingly
        current_V = current_V * arma::diagmat(1 / arma::sum(current_V, 0));

        // Now time consuming part. Calculate svd of new V to get H col 
        Rcpp::Rcout << "SVD for row -> col" << ".\n";

        arma::svd(S_t,Sigma, R_t, current_V);
        S_t = S_t.head_cols(K); // m*r -> r*m
        R_t = R_t.head_cols(K); // n*r  -> r*n
        Sigma = Sigma.head(K); // r
        // project current_V into S_t
        projected_points = current_V.t() * S_t;
        // project current W_col col into S_t
        current_solution = W_col.t() *  S_t;
        Rcpp::Rcout << " --- Solution Coordinates are " << ".\n";
        Rcpp::Rcout << current_solution << ".\n";
        Rcpp::Rcout << "Get relative coordinates for col " << ".\n";
        H_col = get_relative_coordinates_closest(projected_points, current_solution);
        H_col = H_col.t();

        // Now we have W_col H_col
        // tansform H
        H_col = H_col * arma::diagmat(1 / D_vs_col.col(i));
        // row normalize H_col -- is it valid ?
        D_h = 1 / arma::sum(H_col, 1);
        H_row =  arma::diagmat(D_h) * H_col;

        // transform V
        current_V = current_V *  arma::diagmat(1 / D_vs_col.col(i));
        current_V = arma::diagmat( 1 / arma::sum(current_V, 1)) * current_V;
        Rcpp::Rcout << "SVD for col -> row" << ".\n";

        // Time consuming part. again svd ...
        arma::svd(S_t,Sigma, R_t, current_V);
        S_t = S_t.head_cols(K); // m*r 
        R_t = R_t.head_cols(K); // n*r
        Sigma = Sigma.head(K); // r
         // project current_V into R_t
        projected_points = current_V * R_t;
        // project current H_row col into S_t
        current_solution = H_row *  R_t;
        Rcpp::Rcout << " --- Solution Coordinates are " << ".\n";
        Rcpp::Rcout << current_solution << ".\n";
        Rcpp::Rcout << "Get relative coordinates for row " << ".\n";

        W_row = get_relative_coordinates_closest(projected_points, current_solution);
    }
    Rcpp::Rcout << "Loop is done " << ".\n";

    // now we have W_row and H_row and 1 more normalization left
    // matrix needed to column normalize W row will be used as intermediate multiplier to avoid overflow

    // In case V is column normalized we expect points to form a simplex with H being column normalized
    // For this we just repeat col norm H matrix search once more time
    // Transformation for W
    W_row = arma::diagmat(1 / D_vs_row.col(0)) * W_row; // this is not row norm anymore.
    // We omit column normalization
    W_col = W_row ;
    // Transformation for V
    current_V = arma::diagmat(1 / D_vs_row.col(0)) * current_V;
    // V should be column normalized -- transform it accordingly
    current_V = current_V * arma::diagmat(1 / arma::sum(current_V, 0));
    Rcpp::Rcout << "Final SVD for row -> col" << ".\n";
    arma::svd(S_t,Sigma, R_t, current_V);
    S_t = S_t.head_cols(K); // m*r -> r*m
    R_t = R_t.head_cols(K); // n*r  -> r*n
    Sigma = Sigma.head(K); // r
    // project current_V into S_t
    projected_points = current_V.t() * S_t;
        // project current W_col col into S_t
    current_solution = W_col.t() *  S_t;
    Rcpp::Rcout << " --- Solution Coordinates are " << ".\n";
    Rcpp::Rcout << current_solution << ".\n";
    Rcpp::Rcout << "Get relative coordinates for col " << ".\n";
    H_col = get_relative_coordinates_closest(projected_points, current_solution);
    H_col = H_col.t();
    // Get row normalizing matrix for H_col
    D_h = 1 / arma::sum(H_col, 1);
    // This matrix goes to W_row
    W_col = W_row  * arma::diagmat(D_h) ;
    // Now H_col is column normalized coordinates of points between W_col
    return Rcpp::List::create( 
                            // W_row and H_row represent true general matrix factorization
                            Rcpp::Named("W") = W_row,
                            Rcpp::Named("H") = H_row,
                            // H_col represent closest column normalized coefficients for points
                            Rcpp::Named("H_col") = H_col,
                            Rcpp::Named("W_corrected") = W_col
    );
}


Rcpp::List clean_reverse_sinkhorn_c(const arma::mat& result_H_col,
                              const arma::mat& result_W_row,
                              const arma::mat& D_vs_row,
                              const arma::mat& D_vs_col,
                              int iterations) {

    arma::mat H_col = result_H_col;
    arma::mat W_row = result_W_row;
    arma::mat H_row;
    arma::mat W_col;
    arma::mat D_w , D_h; //intermediate matrices to move some magnitude between matrices

    D_h = 1 / arma::sum(H_col, 1);

    // V_fs = W_ss * D_w * D_w_inv * H_ss * D_v_last
    H_row = arma::diagmat(D_h) * H_col;
    // now we have W_ss * H_ss

    // start iteratively update them
    for (int i = iterations - 2; i >= 0; i--) {
        W_row = arma::diagmat(1 / D_vs_row.col(i+1)) * W_row;
        // matrix needed to column normalize W row will be used as intermediate multiplier to avoid overflow
        D_w = 1 / arma::sum(W_row, 0);

        W_col = W_row *  arma::diagmat(D_w);
        H_col =  arma::diagmat(1 / D_w) *  H_row;


        // now deal with col norm matrices
        H_col = H_col * arma::diagmat(1 / D_vs_col.col(i));
        // matrix needed to row normalize H_col will be used as intermediate multiplier to avoid overflow
        D_h = 1 / arma::sum(H_col, 1);
        H_row =  arma::diagmat(D_h) * H_col;
        W_row = W_col * diagmat(1 / D_h);

    }
    // now we have W_row and H_row and 1 more normalization left
    W_row = arma::diagmat(1 / D_vs_row.col(0)) * W_row; // this is not row norm anymore.

    // matrix needed to column normalize W row will be used as intermediate multiplier to avoid overflow
    D_w = 1 / arma::sum(W_row, 0);
    W_col = W_row *  arma::diagmat(D_w);
    H_col =  diagmat(1 / D_w) *  H_row;
    return Rcpp::List::create(Rcpp::Named("W") = W_col,
                              Rcpp::Named("H") = H_col,
                              Rcpp::Named("Dv_inv_W_row") = W_row,
                              Rcpp::Named("H_row") = H_row);
}


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

    D_w_inv_pred_vec = nnls_nonzero_C__(W_col, ones_like_W);
    W_row = W_col * diagmat(D_w_inv_pred_vec);

    for (int i = iterations - 1; i >= 0; i--) {
        D_h_inv_pred_vec = nnls_nonzero_C__(H_row.t(), ones_like_H);

        D_ws_col.col(i) = 1 / D_w_inv_pred_vec;
        D_hs_row.col(i) = 1 / D_h_inv_pred_vec;
        H_col = arma::diagmat(D_h_inv_pred_vec) * H_row;
        W_col = arma::diagmat(1 / D_vs_row.col(i)) * W_row * arma::diagmat(1 / D_h_inv_pred_vec);

        if (i != 0) {
            D_w_inv_pred_vec = nnls_nonzero_C__(W_col, ones_like_W);
            W_row = W_col * diagmat(D_w_inv_pred_vec);
            H_row = diagmat(1 / D_w_inv_pred_vec) * H_col * arma::diagmat(1 / D_vs_col.col(i - 1));
        }
        else {
            // fair estimate of factorization of V using W_row and H_row. For NMF we don't want H to be colnorm
            W_row = arma::diagmat(1 / D_vs_row.col(0)) * W_row; // this is not row norm anymore.
        }
    }

    return Rcpp::List::create(Rcpp::Named("W") = W_col,
                              Rcpp::Named("H") = H_col,
                              Rcpp::Named("Dv_inv_W_row") = W_row,
                              Rcpp::Named("H_row") = H_row,
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


Rcpp::List efficient_sinkhorn(const arma::mat& V,
                              const int max_iter,
                              const int iter_start_check,
                              const int check_every_iter,
                              const double epsilon) {
    // Setup
    double M = V.n_rows;
    double N = V.n_cols;
    double delta = M / N;

    // Initial normalization matrix with 1s.
    arma::vec D_row_sum_current(M);
    arma::mat D_row(M, max_iter + 1, arma::fill::ones);

    arma::rowvec D_col_sum_current(N);
    arma::mat D_col(max_iter + 1, N, arma::fill::ones);

    arma::mat V_ = V;

    // for convergence check
    arma::rowvec converged_col_sum(N, arma::fill::value(1/delta));
    bool converged = false;

    // Main algorithm
    int i;
    for (i = 0; i < max_iter; i++) {
        // Row normalize
        D_row_sum_current = 1 / arma::sum(V_, 1);
        D_row.col(i) = D_row_sum_current;
        V_.each_col() %= D_row_sum_current;

        // Column normalize
        D_col_sum_current = 1 / arma::sum(V_, 0);
        D_col.row(i) = D_col_sum_current;
        V_.each_row() %= D_col_sum_current;

        // Check convergence
        if ((i+1) >= iter_start_check && ((i + 1 -iter_start_check) % check_every_iter) == 0) {
            converged = arma::approx_equal(D_col_sum_current, converged_col_sum, "absdiff", epsilon);
            if (converged) break;  // Only check convergence as the state is updated 
        }
    }

    if (converged) {
        Rcpp::Rcout << "Sinkhorn transformation converge at iteration: " << i << ".\n";
    } else {
        Rcpp::warning("Sinkhorn transformation does not converge at iteration %i", i);
    }

    // will return all 1 columns for D_vs_row and D_vs_col if no normalizations performed
    return Rcpp::List::create(Rcpp::Named("D_vs_row") = (i > 0) ? D_row.cols(0, i - 1) :  D_row.cols(0,0),
                              Rcpp::Named("D_vs_col") = (i > 0) ? D_col.rows(0, i - 1).t() : D_col.rows(0,0).t(),
                              Rcpp::Named("iterations") = i);

}

arma::mat sinkhorn_sweep_c(const arma::mat& V,
                           const arma::mat& D_vs_row,
                           const arma::mat& D_vs_col,
                           unsigned int iter,
                           unsigned int do_last_step) {
    // This function simply does D_v_2n-2 * D_v_2n-4 * ... * D_v_0 * V * D_v_1 * ... * D_v_2n-1
    // iteratively. Iterative scaling is to prevent potential floating-point underflow (remember
    // that during scaling, the value getting smaller and smaller).

    arma::mat V_ = V;
    if (iter > 0) {
        iter -= 1;  // adjust 1-base to 0-base
        // scaling till iter-1 round of normalization
        for (unsigned int i = 0; i < iter; i++) {
            V_.each_col() %= D_vs_row.col(i);
            V_.each_row() %= D_vs_col.col(i).t();
        }
        // Last normalization, returning V_row or V_column is controled by the size of D_vs_row and D_vs_col
        V_.each_col() %= D_vs_row.col(iter);
        if (do_last_step == 1) {
            V_.each_row() %= D_vs_col.col(iter).t();
        }
    }

    return(V_);
}

Rcpp::List extended_sinkhorn(const arma::mat& V,
                             const arma::mat& W,
                             const arma::mat& H,
                             const int n_iter) {
    // Setup
    double M = V.n_rows;
    double N = V.n_cols;
    double K = W.n_cols;

    arma::mat D_v_row(M, n_iter + 1, arma::fill::ones);
    arma::mat D_v_col(N, n_iter + 1, arma::fill::ones);
    arma::mat D_w_col(K, n_iter + 1, arma::fill::ones);
    arma::mat D_h_row(K, n_iter + 1, arma::fill::ones);


    arma::vec D_v_row_sum_current(M);
    arma::vec D_v_col_sum_current(N);
    arma::vec D_w_col_sum_current(K);
    arma::vec D_h_row_sum_current(K);

    arma::mat V_col = V;
    arma::mat V_row = V;
    arma::mat W_row = W;
    arma::mat W_col = W;
    arma::mat H_row = H;
    arma::mat H_col = H;

    int i;
    for (i = 0; i < n_iter; i++) {
        // Row normalize
        D_v_row_sum_current = 1 / arma::sum(V_col, 1);
        D_h_row_sum_current = 1 / arma::sum(H_col, 1);

        D_v_row.col(i) = D_v_row_sum_current;
        D_h_row.col(i) = D_h_row_sum_current;

        V_row = V_col;
        V_row.each_col() %= D_v_row_sum_current;
        H_row = H_col;
        H_row.each_col() %= D_h_row_sum_current;
        W_row = W_col;
        W_row.each_col() %= D_v_row_sum_current;
        W_row.each_row() /= D_h_row_sum_current.t();

        // Column normalize

        D_v_col_sum_current = 1 / arma::sum(V_row, 0).t();
        D_w_col_sum_current = 1 / arma::sum(W_row, 0).t();

        D_w_col.col(i) = D_w_col_sum_current;
        D_v_col.col(i) = D_v_col_sum_current;

        V_col = V_row;
        V_col.each_row() %= D_v_col_sum_current.t();

        W_col = W_row;
        W_col.each_row() %= D_w_col_sum_current.t();
        H_col = H_row;
        H_col.each_row() %= D_v_col_sum_current.t();
        H_col.each_col() /= D_w_col_sum_current;
    }

    // will return all 1 columns for D_vs_row and D_vs_col if no normalizations performed
    return Rcpp::List::create(Rcpp::Named("V_row") = V_row,
                              Rcpp::Named("V_col") = V_col,
                              Rcpp::Named("W_row") = W_row,
                              Rcpp::Named("W_col") = W_col,
                              Rcpp::Named("H_row") = H_row,
                              Rcpp::Named("H_col") = H_col,
                              Rcpp::Named("D_vs_row") = (i > 0) ? D_v_row.cols(0, i - 1) : D_v_row.cols(0,0),
                              Rcpp::Named("D_vs_col") = (i > 0) ? D_v_col.cols(0, i - 1) : D_v_col.cols(0,0),
                              Rcpp::Named("D_hs_row") = (i > 0) ? D_h_row.cols(0, i - 1) : D_h_row.cols(0,0),
                              Rcpp::Named("D_ws_col") = (i > 0) ? D_w_col.cols(0, i - 1) : D_w_col.cols(0,0)
                              );
}

