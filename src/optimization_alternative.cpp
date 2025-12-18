#include "optimization_alternative.h"

#include "matrix_utils.h"
#include "nnls.h"
#include "optimization.h"
#include <tuple>


arma::mat alternative_hinge_der_basis_C__(const arma::mat& W, const arma::mat& S, double precision_) {
    int n = W.n_cols;

    arma::mat res(n, n, arma::fill::zeros);

    for (int j = 0; j < n; j++) {
        arma::vec t = W.col(j);
        res.col(j) = arma::sum(-S.cols(find(t < -precision_)), 1);
    }
    return res;
}


std::tuple<arma::mat, arma::mat, arma::mat> ensure_D_integrity_c(
                              const arma::mat& X_dtilde,
                              const arma::mat& Omega_dtilde,
                              const arma::vec sqrt_Sigma,
                              const double N,
                              const double M) {
    arma::mat new_D_w, new_D_h;
    arma::mat new_D_w_sqrt;
    arma::mat new_D_w_x_sqrt, new_D_w_omega_sqrt;
    arma::mat corrected_X_dtilde, corrected_Omega_dtilde;
    // Get matrix D estimate from X
    new_D_w_x_sqrt =  X_dtilde.col(0) * sqrt_Sigma.at(0) * sqrt(N);
    // Get matrix D estimate from Omega
    new_D_w_omega_sqrt =  Omega_dtilde.row(0).as_col() * sqrt_Sigma.at(0) * sqrt(M);
    // Combine estimates into single matrix
    new_D_w = new_D_w_x_sqrt % new_D_w_omega_sqrt;
    new_D_w_sqrt = arma::sqrt(new_D_w);
    corrected_X_dtilde = arma::diagmat(1/new_D_w_x_sqrt)* arma::diagmat(new_D_w_sqrt) * X_dtilde;
    corrected_Omega_dtilde = Omega_dtilde * arma::diagmat(1/new_D_w_omega_sqrt) * arma::diagmat(new_D_w_sqrt);

    return {corrected_X_dtilde, corrected_Omega_dtilde, new_D_w_sqrt};
}

Rcpp::List ensure_D_integrity(const arma::mat& X_dtilde,
                              const arma::mat& Omega_dtilde,
                              const arma::vec sqrt_Sigma,
                              const double N,
                              const double M) {
    arma::mat new_X = X_dtilde;
    arma::mat new_Omega = Omega_dtilde;
    arma::mat new_D_w_sqrt;
    std::tie(new_X, new_Omega, new_D_w_sqrt) = ensure_D_integrity_c(new_X, new_Omega, sqrt_Sigma, N, M);

    return Rcpp::List::create(Rcpp::Named("corrected_X_dtilde") = new_X,
                              Rcpp::Named("corrected_Omega_dtilde") = new_Omega,
                              Rcpp::Named("new_D_w_sqrt") = new_D_w_sqrt);
}




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
                             const double r_const_X,
                             const double r_const_Omega,
                             const double thresh,
                             const double solution_balancing_threshold) {
    arma::mat errors_statistics(iterations, 9, arma::fill::zeros);
    arma::mat points_statistics_X(iterations, cell_types * cell_types, arma::fill::zeros);
    arma::mat points_statistics_Omega(iterations, cell_types * cell_types, arma::fill::zeros);

    arma::mat new_X = X;
    arma::mat new_Omega = Omega;
    arma::mat temporary_new_X = X;
    arma::mat temporary_new_Omega = Omega;
    arma::mat final_X = X;
    arma::mat final_Omega = Omega;
    arma::mat new_D_w = D_w;
    arma::mat new_D_w_sqrt = arma::sqrt(new_D_w);
    arma::mat temporary_new_D_w_sqrt =arma::sqrt(new_D_w);

    arma::mat new_D_h = new_D_w * (N / M);

    arma::vec Sigma = arma::diagvec(SVRt);
    arma::vec sqrt_Sigma = arma::sqrt(Sigma);


    new_X =  arma::diagmat(new_D_w_sqrt) * new_X * arma::diagmat(1 / sqrt_Sigma);
    new_Omega =  arma::diagmat(1 / sqrt_Sigma) *  new_Omega  * arma::diagmat(new_D_w_sqrt);
    arma::mat jump_X, jump_Omega;

    arma::vec vectorised_SVRt = arma::vectorise(SVRt);
    arma::colvec sum_rows_R = arma::sum(R, 1);
    arma::colvec sum_rows_S = arma::sum(S, 1);

    arma::mat B = join_cols(vectorised_SVRt, coef_pos_D_w * sum_rows_S);
    arma::mat C = join_cols(vectorised_SVRt, coef_pos_D_h * sum_rows_R);
    arma::mat der_X, der_Omega;
    arma::mat tmp_X, tmp_Omega;
    double shrink_limit = 500;
    double mean_norm_solution_X;
    double mean_norm_solution_Omega;

    for (int itr_ = 0; itr_ < iterations; itr_++) {
        //  derivative X
        //  der_X = -2 * (new_Omega.t() * (SVRt - new_Omega * new_X));
        //  der_X +=  coef_hinge_H * hinge_der_proportions_C__(new_X  * R, R);
        //  der_X = -2 * new_X;
        der_X =  coef_hinge_H * hinge_der_proportions_C__(new_X  * arma::diagmat(sqrt_Sigma)  * R, R) * arma::diagmat(1 / sqrt_Sigma);
        der_X += 2 * new_X;
        //  Rcpp::Rcout << "original der X"  << std::endl;
        //  Rcpp::Rcout << der_X << std::endl;

        mean_norm_solution_X = arma::mean(arma::vecnorm(new_X, 2, 1));
        der_X = correctByNorm(der_X) * mean_norm_solution_X;

        //  Rcpp::Rcout << " der X"  << std::endl;
        //  Rcpp::Rcout << der_X << std::endl;
        tmp_X = (new_X - coef_der_X * der_X); // estimate new X given derivative
        //  Rcpp::Rcout << " X candidate"  << std::endl;
        //  Rcpp::Rcout << tmp_X << std::endl;

        if (any( tmp_X.col(0) <= 0)) {
        //  Rcpp::Rcout << "Derrivative caused negative for X \n"  << std::endl;
            for (int c=0; c < cell_types; c++) {
                double matrix_value =  tmp_X(c,0);
                 if (matrix_value <= 0) {
                   int shrink_iteration = 0;
                   while(matrix_value <= 0) {
                    der_X.row(c) /=  2;
                    tmp_X = (new_X - coef_der_X * der_X);
                    matrix_value =  tmp_X(c,0);
                    shrink_iteration++;
                   }
                 }
            }
            if  (any( tmp_X.col(0) <= 0)) {
                Rcpp::Rcout << "Any gradient step gives bad X, probably X was bad before\n"  << std::endl;

            }
        //  Rcpp::Rcout << "X shrink iterations performed for component " << c << " is: " << shrink_iteration << "\n";
        }
       //  Now estimate omega
       tmp_Omega = arma::pinv(tmp_X);
       //  Rcpp::Rcout << " Omega candidate"  << std::endl;
       //  Rcpp::Rcout << tmp_Omega << std::endl;
        if (any( tmp_Omega.row(0) <= 0)) {
       //  Rcpp::Rcout << "Inverse of X caused negative for Omega \n"  << std::endl;
            for (int c=0; c < cell_types; c++) {
                double matrix_value =  tmp_Omega(0,c);
                if (matrix_value <= 0) {
                    int shrink_iteration = 0;
                    while((matrix_value <= 0)& (shrink_iteration < shrink_limit)) {
                    der_X /=  2;
                    der_X.row(c) *= 2;
                    tmp_X = (new_X - coef_der_X * der_X);
                    tmp_Omega = arma::pinv(tmp_X);
                    matrix_value =  tmp_Omega(0,c);
                    shrink_iteration++;
                   }
                if (shrink_iteration != shrink_limit) {
                    // if we were able to find the solution. accept these new X and Omega
                    new_Omega = tmp_Omega;
                    new_X = tmp_X;
                    } else {
                        Rcpp::Rcout << "Couldn't find good inverse X, reject X for one of the steps\n"  << std::endl;
                    }
                }
            }
        } else {
            new_Omega = tmp_Omega;
            new_X = tmp_X;
        }
        // continue conventional optimization
        // at this stage we are sure that first first row/col of X/omega are positive
        // Correct X and Omega to have corresponding first row/column and get D matrix
//       Rcpp::Rcout << " Before correction X"  << std::endl;
//       Rcpp::Rcout << new_X << std::endl;
//       Rcpp::Rcout << " Before correction Omega"  << std::endl;
//       Rcpp::Rcout << new_Omega << std::endl;
       std::tie(new_X, new_Omega, new_D_w_sqrt) = ensure_D_integrity_c(new_X, new_Omega, sqrt_Sigma, N, M);
//       Rcpp::Rcout << " After correction X"  << std::endl;
//       Rcpp::Rcout << new_X << std::endl;
//       Rcpp::Rcout << " After correction Omega"  << std::endl;
//       Rcpp::Rcout << new_Omega << std::endl;

        // Ensure Omega and X vectors norms are  adequate
        // Get temporary final X and Omega
       final_X = arma::diagmat(1/new_D_w_sqrt) * new_X * arma::diagmat(sqrt_Sigma);
       final_Omega = arma::diagmat(sqrt_Sigma)* new_Omega * arma::diagmat(1/new_D_w_sqrt);

//       Rcpp::Rcout << " Final  X for X update"  << std::endl;
//       Rcpp::Rcout << final_X << std::endl;
//       Rcpp::Rcout << " Final  Omega for X update"  << std::endl;
//       Rcpp::Rcout << final_Omega << std::endl;

        for (int c=0; c < cell_types; c++) {
            double col_omega_norm = arma::norm(final_Omega.col(c).subvec(1, cell_types - 1), 2);
            double row_x_norm = arma::norm(final_X.row(c).subvec(1, cell_types - 1), 2);
            if (col_omega_norm > solution_balancing_threshold * mean_radius_Omega) {
                Rcpp::Rcout << "Looks like Omega points are way far away after inverse of X. \n"  << std::endl;
                Rcpp::Rcout << "We will balance solution by moving some magnitude from Omega to X. \n"  << std::endl;
                double ratio_x = row_x_norm / mean_radius_X;
                double ratio_omega = col_omega_norm / mean_radius_Omega;
                // stretch the space for all elements
                double multiplier_x = sqrt(ratio_omega)/sqrt(ratio_x);
                double multiplier_omega =  sqrt(ratio_x)/sqrt(ratio_omega);
                Rcpp::Rcout << " X cols multiplied by " << multiplier_x << std::endl;
                Rcpp::Rcout << " Omega cols multiplied by " << multiplier_omega << std::endl;

                for (int column=1; column < cell_types; column++) {
                    std::tie(temporary_new_X, temporary_new_Omega, temporary_new_D_w_sqrt) = ensure_D_integrity_c(new_X, new_Omega, sqrt_Sigma, N, M);
                    final_X = arma::diagmat(1/new_D_w_sqrt) * new_X * arma::diagmat(sqrt_Sigma);
                    final_Omega = arma::diagmat(sqrt_Sigma)* new_Omega * arma::diagmat(1/new_D_w_sqrt);
                    Rcpp::Rcout << "row " << column << " balanced " << std::endl;
                    new_X.col(column) *= multiplier_x;
                    new_Omega.row(column) *= multiplier_omega;
                }
            }
            Rcpp::Rcout << "Finish norm checking. \n"  << std::endl;

        }
// Correct X and Omega to have corresponding first row/column and derrive D
        std::tie(new_X, new_Omega, new_D_w_sqrt) = ensure_D_integrity_c(new_X, new_Omega, sqrt_Sigma, N, M);

//       Rcpp::Rcout << " After final correction X"  << std::endl;
//       Rcpp::Rcout << new_X << std::endl;
//       Rcpp::Rcout << " After final correction Omega"  << std::endl;
//       Rcpp::Rcout << new_Omega << std::endl;

        // derivative Omega
        // der_Omega = -2 * new_Omega;
        der_Omega = coef_hinge_W * arma::diagmat(1 / sqrt_Sigma) * alternative_hinge_der_basis_C__(S.t() * arma::diagmat(sqrt_Sigma) * new_Omega, S);
        mean_norm_solution_Omega = arma::mean(arma::vecnorm(new_Omega, 2, 2));
        der_Omega = correctByNorm(der_Omega) * mean_norm_solution_Omega;
        der_Omega += 2 * new_Omega;

        tmp_Omega = new_Omega - coef_der_Omega * der_Omega;
//        Rcpp::Rcout << " Omega candidate"  << std::endl;
//        Rcpp::Rcout << tmp_Omega << std::endl;

        if (any(tmp_Omega.row(0) <= 0)) {
        //Rcpp::Rcout << "Derrivative caused negative for Omega \n"  << std::endl;
        // start shrinking derivative to be inside the range
            for (int c=0; c < cell_types; c++) {
                double matrix_value =  tmp_Omega(0,c);
                if (matrix_value <= 0) {
                    int shrink_iteration = 0;
                    while(matrix_value <= 0) {
                        der_Omega.col(c) /=  2;
                        tmp_Omega = new_Omega - coef_der_Omega * der_Omega;
                        matrix_value =  tmp_Omega(0,c);
                        shrink_iteration++;
                       }
                    }
            }

            if  (any( tmp_Omega.row(0) <= 0)) {
                Rcpp::Rcout << "Any gradient step gives bad Omega, probably Omega was bad before\n"  << std::endl;
            }
        }
        // Now try estimate X as inverse
        tmp_X = arma::pinv(tmp_Omega);
        if (any(tmp_X.col(0) <= 0)) {
            for (int c=0; c < cell_types; c++) {
                double matrix_value =  tmp_X(c,0);
                if (matrix_value < 0) {
                    int shrink_iteration = 0;
                    while((matrix_value <= 0) & (shrink_iteration < shrink_limit)) {
                        der_Omega /=  2;
                        der_Omega.col(c) *=  2;
                        tmp_Omega = new_Omega - coef_der_Omega * der_Omega;
                        tmp_X = arma::pinv(tmp_Omega);
                        matrix_value =  tmp_X(c,0);
                        shrink_iteration++;
                    }
                    if (shrink_iteration != shrink_limit) {
                        // if we were able to find the solution. accept these new X and Omega
                        new_Omega = tmp_Omega;
                        new_X = tmp_X;
                    } else {
                        Rcpp::Rcout << "Couldn't find good inverse Omega, reject Omega for one of the steps \n"  << std::endl;
                    }
                }
            }
        } else {
            new_Omega = tmp_Omega;
            new_X = tmp_X;
        }
        // Correct X and Omega to have corresponding first row/column and derrive D
        std::tie(new_X, new_Omega, new_D_w_sqrt) = ensure_D_integrity_c(new_X, new_Omega, sqrt_Sigma, N, M);
        // continue conventional optimization
        // at this stage we are sure that first first row/col of X/omega are positive
        // Ensure Omega vectors norm is adequate
        // Get temporary final X and Omega
        final_X = arma::diagmat(1/new_D_w_sqrt) * new_X * arma::diagmat(sqrt_Sigma);
        final_Omega = arma::diagmat(sqrt_Sigma)* new_Omega * arma::diagmat(1/new_D_w_sqrt);
//        for (int c=0; c < cell_types; c++) {
//            double col_omega_norm = arma::norm(final_Omega.col(c).subvec(1, cell_types - 1), 2);
//            double row_x_norm = arma::norm(final_X.row(c).subvec(1, cell_types - 1), 2);
//            if (row_x_norm > solution_balancing_threshold * mean_radius_X) {
//                Rcpp::Rcout << "Looks like X point are way far away after inverse of Omega. \n"  << std::endl;
//                Rcpp::Rcout << "We will balance solution by moving some magnitude from X to Omega. \n"  << std::endl;
//                double ratio_x = row_x_norm / mean_radius_X;
//                double ratio_omega = col_omega_norm / mean_radius_Omega;
//                double multiplier_x = sqrt(ratio_omega)/sqrt(ratio_x);
//                double multiplier_omega =  sqrt(ratio_x)/sqrt(ratio_omega);
//
////                Rcpp::Rcout << " X row" << final_X.row(c) << std::endl;
////                Rcpp::Rcout << " Omega col" << new_Omega.col(c)<< std::endl;
////                Rcpp::Rcout << " ratio X" << ratio_x<< std::endl;
////                Rcpp::Rcout << " ratio Omega" << ratio_omega<< std::endl;
//                for (int column=1; column < cell_types; column++) {
//                    std::tie(temporary_new_X, temporary_new_Omega, temporary_new_D_w_sqrt) = ensure_D_integrity_c(new_X, new_Omega, sqrt_Sigma, N, M);
//                    final_X = arma::diagmat(1/new_D_w_sqrt) * new_X * arma::diagmat(sqrt_Sigma);
//                    final_Omega = arma::diagmat(sqrt_Sigma)* new_Omega * arma::diagmat(1/new_D_w_sqrt);
//                    if (final_X.col(column).max()/mean_radius_X > solution_balancing_threshold * mean_radius_X) {
//                        Rcpp::Rcout << "row " << column << " balanced " << std::endl;
//                        new_X.col(column) *= multiplier_x;
//                        new_Omega.row(column) *= multiplier_omega;
//                    }
//
//                }
//            }
//        }
        // Correct X and Omega to have corresponding first row/column and derrive D
        std::tie(new_X, new_Omega, new_D_w_sqrt) = ensure_D_integrity_c(new_X, new_Omega, sqrt_Sigma, N, M);
        new_D_w = arma::pow(new_D_w_sqrt, 2);
        new_D_h = new_D_w * (N / M);

        double sum_ = accu(new_D_w) / M;

        // result X and omega
        final_X = arma::diagmat(1/new_D_w_sqrt) * new_X * arma::diagmat(sqrt_Sigma);
        final_Omega = arma::diagmat(sqrt_Sigma)* new_Omega * arma::diagmat(1/new_D_w_sqrt);

        arma::uword neg_props = getNegative(final_X * R);
        arma::uword neg_basis = getNegative(S.t() * final_Omega);

        Rcpp::List current_errors = calcErrors(final_X,
                                               final_Omega,
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


        points_statistics_X.row(itr_) = final_X.as_row();
        points_statistics_Omega.row(itr_) = final_Omega.as_row();
    }


    return Rcpp::List::create(Rcpp::Named("new_X") = final_X,
                              Rcpp::Named("new_Omega") = final_Omega,
                              Rcpp::Named("new_D_w") = new_D_w,
                              Rcpp::Named("new_D_h") = new_D_h,
                              Rcpp::Named("errors_statistics") = errors_statistics,
                              Rcpp::Named("points_statistics_X") = points_statistics_X,
                              Rcpp::Named("points_statistics_Omega") = points_statistics_Omega);
}
