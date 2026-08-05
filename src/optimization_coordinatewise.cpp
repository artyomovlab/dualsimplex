#include "optimization_utils.h"
#include "matrix_utils.h"
#include "nnls.h"
#include "optimization_coordinatewise.h"




Rcpp::List optimize_coordinate_descent(const arma::mat& X,
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
                             const double thresh,
                             const double convergence_tol,
                             const int stop_criteria_window,
                             const bool debug_stats) {
    arma::mat errors_statistics(iterations + 1, 23, arma::fill::zeros);
    arma::mat points_statistics_X(iterations + 1, cell_types * cell_types, arma::fill::zeros);
    arma::mat points_statistics_Omega(iterations + 1, cell_types * cell_types, arma::fill::zeros);
    arma::mat points_statistics_Dw(iterations + 1, cell_types, arma::fill::zeros);


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
    double current_learning_rate_X = coef_der_X;
    double current_learning_rate_Omega = coef_der_Omega;
    double best_error_value = 10000;
    int best_error_iteration = 0;
    double current_error_value;
    double average_gradient_norm = 0;
    double average_hinge_H_gradient_norm = 0;
    double average_hinge_W_gradient_norm = 0;
    double average_reg_X_gradient_norm = 0;
    double average_reg_Omega_gradient_norm = 0;

    arma::mat hinge_term_H, hinge_term_W, der_term_deconv_X, der_term_deconv_Omega;

                        
    int itr_ = 0;
    if ((current_learning_rate_X <= convergence_tol) & (current_learning_rate_Omega <= convergence_tol)) {
        spdl::error("Initial learning rate can not be smaller than convergence_tol");
    }                            
    while ((itr_ < iterations + 1) & (current_learning_rate_X > convergence_tol) & (current_learning_rate_Omega > convergence_tol)) {
        if (itr_ > 0) {
            // derivative X
            der_term_deconv_X =  -2 * (diagmat(new_D_w) * new_Omega.t() * (SVRt - new_Omega * diagmat(new_D_w) * new_X));
            hinge_term_H = l1_hinge_der_proportions_C__(new_X * R, R);
            der_X = der_term_deconv_X + coef_hinge_H * hinge_term_H;
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
            new_X = new_X - current_learning_rate_X * der_X;
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
            der_term_deconv_Omega =  -2 * (SVRt - new_Omega * diagmat(new_D_w) * new_X) * new_X.t() * diagmat(new_D_w);
            hinge_term_W = l1_hinge_der_basis_C__(S.t() * new_Omega, S);
            der_Omega = der_term_deconv_Omega + coef_hinge_W * hinge_term_W;
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

            // Update Omega
            new_Omega = new_Omega - current_learning_rate_Omega * der_Omega;
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
        }

        arma::uword neg_props = getNegative(new_X * R);
        arma::uword neg_basis = getNegative(S.t() * new_Omega);
        double sum_ = accu(new_D_w) / M;

        if (debug_stats) {
            average_gradient_norm = arma::mean(arma::vecnorm((der_X + der_Omega) / 2, 2, 1));
            average_hinge_H_gradient_norm = arma::mean(arma::vecnorm(hinge_term_H, 2, 1));
            average_hinge_W_gradient_norm = arma::mean(arma::vecnorm(hinge_term_W, 2, 1));
            average_reg_X_gradient_norm = arma::mean(arma::vecnorm(der_term_deconv_X, 2, 1));
            average_reg_Omega_gradient_norm = arma::mean(arma::vecnorm(der_term_deconv_Omega, 2, 1));
        }

        Rcpp::List current_errors = calcErrors(new_X,
                                               new_Omega,
                                               new_D_w,
                                               new_D_h,
                                               SVRt,
                                               R,
                                               S,
                                               coef_hinge_H,
                                               coef_hinge_W,
                                               coef_pos_D_h,
                                               coef_pos_D_w);

        current_error_value = current_errors["total_error"];
        if (current_error_value < best_error_value) {
            best_error_iteration = itr_;
            best_error_value = current_error_value;
        }
        if (itr_ - best_error_iteration > stop_criteria_window) {
            // looks like best solution was not updated for stop_criteria_window iterations. reducing step size.
            if (current_learning_rate_X > convergence_tol) current_learning_rate_X = current_learning_rate_X / 2;
            if (current_learning_rate_Omega > convergence_tol) current_learning_rate_Omega = current_learning_rate_Omega / 2;
            best_error_iteration = itr_; // reset iteration counter
        }


        errors_statistics.row(itr_) = arma::rowvec{current_errors["deconv_error"], //1
                                                   current_errors["lambda_error"], //2
                                                   current_errors["beta_error"],   //3
                                                   current_errors["D_h_error"],    //4
                                                   current_errors["D_w_error"],    //5
                                                   current_errors["total_error"],  //6
                                                   current_errors["scaled_total_error"], //7

                                                   static_cast<double>(neg_props), //8
                                                   static_cast<double>(neg_basis), //9
                                                   sum_,  //10
                                                   current_errors["average_norm"], //11
                                                   (current_learning_rate_X + current_learning_rate_Omega) / 2 , //12
                                                   average_gradient_norm, //13
                                                   average_hinge_H_gradient_norm, //14
                                                   average_hinge_W_gradient_norm, //15
                                                   average_reg_X_gradient_norm, //16
                                                   average_reg_Omega_gradient_norm, //17
                                                   best_error_value, //18
                                                   static_cast<double>(best_error_iteration), //19,
                                                   current_errors["scaled_lambda_error"],            //20
                                                   current_errors["scaled_beta_error"], //21
                                                   average_gradient_norm, //22
                                                   0
                                                   };
        points_statistics_X.row(itr_) = new_X.as_row();
        points_statistics_Omega.row(itr_) = new_Omega.as_row();
        points_statistics_Dw.row(itr_) = new_D_w.as_row();
        itr_++;

    }
    spdl::info("Optimization completed with number of iterations perfomed: {}", itr_ - 1);
    if (itr_ < iterations + 1) {
        points_statistics_X.resize(itr_, points_statistics_X.n_cols);
        points_statistics_Omega.resize(itr_, points_statistics_Omega.n_cols);
        points_statistics_Dw.resize(itr_, points_statistics_Dw.n_cols);
        errors_statistics.resize(itr_, errors_statistics.n_cols);
    }

    return Rcpp::List::create(Rcpp::Named("new_X") = new_X,
                              Rcpp::Named("new_Omega") = new_Omega,
                              Rcpp::Named("new_D_w") = new_D_w,
                              Rcpp::Named("new_D_h") = new_D_h,
                              Rcpp::Named("errors_statistics") = errors_statistics,
                              Rcpp::Named("points_statistics_X") = points_statistics_X,
                              Rcpp::Named("points_statistics_Omega") = points_statistics_Omega,
                              Rcpp::Named("points_statistics_Dw") = points_statistics_Dw);
}
