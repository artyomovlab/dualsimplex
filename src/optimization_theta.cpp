#include "optimization_theta.h"

#include "matrix_utils.h"
#include "nnls.h"
#include "optimization.h"


Rcpp::List theta_derivative_stage2(const arma::mat& X,
                             const arma::mat& Omega,
                             const arma::mat& D_w,
                             const arma::mat& SVRt,
                             const arma::mat& R,
                             const arma::mat& S,
                             const arma::mat& X_center,
                             const arma::mat& Omega_center,
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
                             const double theta_threshold
                            ) {
    Rcpp::Rcout << "sucessfull init training" << "\n";
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

    arma::mat B = join_cols(vectorised_SVRt, coef_pos_D_w * sum_rows_S);
    arma::mat C = join_cols(vectorised_SVRt, coef_pos_D_h * sum_rows_R);
    arma::mat der_X, der_Omega;
    arma::mat tmp_X, tmp_Omega;
    double cos_theta = std::cos(theta_threshold);


    for (int itr_ = 0; itr_ < iterations; itr_++) {
        // derivative X
        //der_X = -2 * (new_Omega.t() * (SVRt - new_Omega * new_X));
        //der_X +=  coef_hinge_H * hinge_der_proportions_C__(new_X  * R, R);

        der_X =
            -2 * (diagmat(new_D_w) * new_Omega.t() * (SVRt - new_Omega * diagmat(new_D_w) * new_X));
        der_X += coef_hinge_H * hinge_der_proportions_C__(new_X * R, R);
        der_X += coef_pos_D_h * 2 * new_D_h * (new_X.t() * new_D_h - sum_rows_R).t();
        der_X.col(0).zeros(); // no movement needed for first coordinate
        der_X = correctByNorm(der_X) * mean_radius_X;

        if (theta_threshold > 0) {
            tmp_X = (new_X - coef_der_X * der_X);
            for (int c=0; c < cell_types; c++) {
                if (!(X_center.row(c).subvec(1, cell_types - 1).is_zero())) {
                    double cos_distance_result = cosine_distance(tmp_X.row(c).subvec(1, cell_types - 1), X_center.row(c).subvec(1, cell_types - 1));
                    if (cos_distance_result < cos_theta) {
                        // start shrinking derivative to be inside
                        int shrink_iteration = 0;
                        while(cos_distance_result < cos_theta) {
                            der_X.row(c) /=  2;
                            tmp_X = (new_X - coef_der_X * der_X);
                            cos_distance_result = cosine_distance(tmp_X.row(c).subvec(1, cell_types - 1), X_center.row(c).subvec(1, cell_types - 1));
                            // Rcout << "Now cos is  : " << cos_distance_result << "\n";
                            shrink_iteration++;
                          }
                    }
                }
            }
        }
        // continue conventional optimization
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
        der_Omega =  -2 * (SVRt - new_Omega * diagmat(new_D_w) * new_X) * new_X.t() * diagmat(new_D_w);
        der_Omega += coef_hinge_W * hinge_der_basis_C__(S.t() * new_Omega, S);
        der_Omega += coef_pos_D_w * 2 * (new_Omega * new_D_w - sum_rows_S) * new_D_w.t();
        der_Omega.row(0).zeros();
        der_Omega = correctByNorm(der_Omega) * mean_radius_Omega;
        Rcpp::Rcout << new_Omega << "\n";
        if (theta_threshold > 0) {
            Rcpp::Rcout << "Start theta business for Omega" << "\n";
            tmp_Omega = new_Omega - coef_der_Omega * der_Omega;
            Rcpp::Rcout << "Current tmp Omega " << tmp_Omega <<"\n";
            for (int c=0; c < cell_types; c++) {
                if (!(Omega_center.col(c).subvec(1, cell_types - 1).is_zero())) {
                    Rcpp::Rcout << "Center point is defineed for " << c << "\n";
                    double cos_distance_result = cosine_distance(tmp_Omega.col(c).subvec(1, cell_types - 1), Omega_center.col(c).subvec(1, cell_types - 1));
                    Rcpp::Rcout << "Cosine  is  " << cos_distance_result  << "\n";
                    if (cos_distance_result < cos_theta) {
                        // start shrinking derivative to be inside
                        Rcpp::Rcout << "Distance is lower than threshold  of " << cos_theta  << ". Start shrink"<< "\n";
                        int shrink_iteration = 0;
                        while(cos_distance_result < cos_theta) {
                            der_Omega.col(c) /=  2;
                            tmp_Omega = new_Omega - coef_der_Omega * der_Omega;
                            cos_distance_result = cosine_distance(tmp_Omega.col(c).subvec(1, cell_types - 1), Omega_center.col(c).subvec(1, cell_types - 1));
                            // Rcout << "Now cos is  : " << cos_distance_result << "\n";
                            shrink_iteration++;
                          }
                          Rcpp::Rcout << "Omega shrink iterations performed : " << shrink_iteration << "\n";
                    }
                }
            }
        }
        // Continue conventional optimization
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
