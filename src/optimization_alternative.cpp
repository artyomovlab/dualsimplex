#include "optimization_alternative.h"

#include "matrix_utils.h"
#include "nnls.h"
#include "optimization.h"

arma::mat alternative_hinge_der_basis_C__(const arma::mat& W, const arma::mat& S, double precision_) {
    int n = W.n_cols;

    arma::mat res(n, n, arma::fill::zeros);

    for (int j = 0; j < n; j++) {
        arma::vec t = W.col(j);
        res.col(j) = arma::sum(-S.cols(find(t < -precision_)), 1);
    }

    return res;
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
    arma::mat final_X = X;
    arma::mat final_Omega = Omega;
    arma::mat new_D_w = D_w;
    arma::mat new_D_w_sqrt = arma::sqrt(new_D_w);
    arma::mat new_D_w_x = D_w;
    arma::mat new_D_w_x_sqrt = arma::sqrt(new_D_w_x);
    arma::mat new_D_w_omega = D_w;
    arma::mat new_D_w_omega_sqrt = arma::sqrt(new_D_w_omega);
    arma::mat new_D_h = new_D_w * (N / M);

    arma::vec Sigma = arma::diagvec(SVRt);
    arma::vec sqrt_Sigma = arma::sqrt(Sigma);


    new_X =  arma::diagmat(new_D_w_x_sqrt) * new_X * arma::diagmat(1 / sqrt_Sigma);
    new_Omega =  arma::diagmat(1 / sqrt_Sigma) *  new_Omega * arma::diagmat(new_D_w_omega_sqrt);



    arma::mat jump_X, jump_Omega;

    arma::vec vectorised_SVRt = arma::vectorise(SVRt);
    arma::colvec sum_rows_R = arma::sum(R, 1);
    arma::colvec sum_rows_S = arma::sum(S, 1);

    arma::mat B = join_cols(vectorised_SVRt, coef_pos_D_w * sum_rows_S);
    arma::mat C = join_cols(vectorised_SVRt, coef_pos_D_h * sum_rows_R);
    arma::mat der_X, der_Omega;
    arma::mat tmp_X, tmp_Omega;
    double shrink_limit = 500;

    for (int itr_ = 0; itr_ < iterations; itr_++) {
        // derivative X
        //der_X = -2 * (new_Omega.t() * (SVRt - new_Omega * new_X));
        //der_X +=  coef_hinge_H * hinge_der_proportions_C__(new_X  * R, R);

        der_X =  coef_hinge_H * hinge_der_proportions_C__(new_X  * arma::diagmat(sqrt_Sigma)  * R, R) * arma::diagmat(1 / sqrt_Sigma);

        der_X = correctByNorm(der_X) * mean_radius_X;


       tmp_X = (new_X - coef_der_X * der_X); // estimate new X given derivative

       if (any( tmp_X.col(0) <= 0)) {
//           Rcpp::Rcout << "Derrivative caused negative for X \n"  << std::endl;
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
           //Rcpp::Rcout << "X shrink iterations performed for component " << c << " is: " << shrink_iteration << "\n";
        }
        // Now estimate omega
       tmp_Omega = arma::pinv(tmp_X);
        if (any( tmp_Omega.row(0) <= 0)) {
//            Rcpp::Rcout << "Inverse of X caused negative for Omega \n"  << std::endl;
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

        // Ensure Omega and X vectors norms are  adequate
        for (int c=0; c < cell_types; c++) {
            double col_omega_norm = arma::norm(new_Omega.col(c).subvec(1, cell_types - 1), 2);
            double row_x_norm = arma::norm(new_X.row(c).subvec(1, cell_types - 1), 2);
            if (col_omega_norm > solution_balancing_threshold * mean_radius_Omega) {
               Rcpp::Rcout << "Looks like Omega points are way far away after inverse of X. \n"  << std::endl;
                Rcpp::Rcout << "We will balance solution by moving some magnitude from Omega to X. \n"  << std::endl;
                double ratio_x = row_x_norm / mean_radius_X;
                double ratio_omega = col_omega_norm / mean_radius_Omega;
                new_X.row(c) *= ratio_omega/ratio_x;
                new_Omega.col(c) *= ratio_x/ratio_omega;
            }

        }

        // Get matrix D
        new_D_w_x_sqrt =  new_X.col(0) * sqrt_Sigma.at(0) * sqrt(N);
        new_D_w_x = arma::pow(new_D_w_x_sqrt, 2);
        new_D_w_omega_sqrt =  new_Omega.row(0).as_col() * sqrt_Sigma.at(0) * sqrt(M);
        new_D_w_omega = arma::pow(new_D_w_omega_sqrt, 2);
        new_D_w = new_D_w_x_sqrt % new_D_w_omega_sqrt;
        new_D_w_sqrt = arma::sqrt(new_D_w);

        // Fix X and omega accordingly
        new_X = arma::diagmat(1/new_D_w_x_sqrt)* arma::diagmat(new_D_w_sqrt) * new_X;
        new_Omega = new_Omega * arma::diagmat(1/new_D_w_omega_sqrt) * arma::diagmat(new_D_w_sqrt);
        new_D_w_x_sqrt = new_D_w_sqrt;
        new_D_w_omega_sqrt = new_D_w_sqrt;





        // derivative Omega

        der_Omega = coef_hinge_W * arma::diagmat(1 / sqrt_Sigma) * alternative_hinge_der_basis_C__(S.t() * arma::diagmat(sqrt_Sigma) * new_Omega, S);
        der_Omega = correctByNorm(der_Omega) * mean_radius_Omega;

        tmp_Omega = new_Omega - coef_der_Omega * der_Omega;

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
        // continue conventional optimization
        // at this stage we are sure that first first row/col of X/omega are positive
        // Ensure Omega vectors norm is adequate
        for (int c=0; c < cell_types; c++) {
            double col_omega_norm = arma::norm(new_Omega.col(c).subvec(1, cell_types - 1), 2);
            double row_x_norm = arma::norm(new_X.row(c).subvec(1, cell_types - 1), 2);
            if (row_x_norm > solution_balancing_threshold * mean_radius_X) {
//                Rcpp::Rcout << "Looks like X point are way far away after inverse of Omega. \n"  << std::endl;
//                Rcpp::Rcout << "We will balance solution by moving some magnitude from X to Omega. \n"  << std::endl;
                double ratio_x = row_x_norm / mean_radius_X;
                double ratio_omega = col_omega_norm / mean_radius_Omega;
                new_X.row(c) *= ratio_omega/ratio_x;
                new_Omega.col(c) *= ratio_x/ratio_omega;
            }
        }

       // Rcpp::Rcout << "going to get D_w from first column" << std::endl;
        new_D_w_omega_sqrt = new_Omega.row(0).as_col() * sqrt_Sigma.at(0) * sqrt(M);
        new_D_w_omega = arma::pow(new_D_w_omega_sqrt, 2);


        new_D_w_x_sqrt =  new_X.col(0) * sqrt_Sigma.at(0) * sqrt(N);
        new_D_w_x = arma::pow(new_D_w_x_sqrt, 2);


        // correct D_w to make it the same matrix
        new_D_w = new_D_w_x_sqrt % new_D_w_omega_sqrt;
        new_D_w_sqrt = arma::sqrt(new_D_w);
        new_X = arma::diagmat(1/new_D_w_x_sqrt)* arma::diagmat(new_D_w_sqrt) * new_X;
        new_Omega = new_Omega * arma::diagmat(1/new_D_w_omega_sqrt) * arma::diagmat(new_D_w_sqrt);
        new_D_w_x_sqrt = new_D_w_sqrt;
        new_D_w_omega_sqrt = new_D_w_sqrt;

        new_D_h = new_D_w * (N / M);

        double sum_ = accu(new_D_w) / M;

        // result X and omega
        final_X = arma::diagmat(1/new_D_w_x_sqrt) * new_X * arma::diagmat(sqrt_Sigma);

        final_Omega = arma::diagmat(sqrt_Sigma)* new_Omega * arma::diagmat(1/new_D_w_omega_sqrt);
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
