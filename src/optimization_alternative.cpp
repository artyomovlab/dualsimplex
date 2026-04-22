#include "optimization_alternative.h"

#include "matrix_utils.h"
#include "nnls.h"
#include "optimization.h"
#include <tuple>


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
                             const double coef_hinge_H,
                             const double coef_hinge_W,
                             const int cell_types,
                             const double N,
                             const double M,
                             const int iterations,
                             const double mean_radius_X,
                             const double mean_radius_Omega,
                             const double solution_balancing_threshold,
                             const double reg_X,
                             const double reg_Omega,
                             const double convergence_tol) {
    arma::mat errors_statistics(iterations, 10, arma::fill::zeros);
    arma::mat points_statistics_X(iterations, cell_types * cell_types, arma::fill::zeros);
    arma::mat points_statistics_Omega(iterations, cell_types * cell_types, arma::fill::zeros);
    // arma::mat points_statistics_X_dtilda_uncorrected(iterations, cell_types * cell_types, arma::fill::zeros);
    // arma::mat points_statistics_Omega_dtilda_uncorrected(iterations, cell_types * cell_types, arma::fill::zeros);
    //arma::mat points_statistics_X_dtilda_corrected(iterations, cell_types * cell_types, arma::fill::zeros);
    //arma::mat points_statistics_Omega_dtilda_corrected(iterations, cell_types * cell_types, arma::fill::zeros);
    arma::mat points_statistics_Dw(iterations, cell_types, arma::fill::zeros);

    arma::mat new_X = X;
    arma::mat new_Omega = Omega;
    arma::mat temporary_new_X = X;
    arma::mat temporary_new_Omega = Omega;
    arma::mat final_X = X;
    arma::mat final_Omega = Omega;
    arma::mat new_D_w = D_w;
    arma::mat new_D_w_sqrt = arma::sqrt(new_D_w);
    arma::mat temporary_new_D_w_sqrt = arma::sqrt(new_D_w);

    arma::mat new_D_h = new_D_w * (N / M);

    arma::vec Sigma = arma::diagvec(SVRt);
    arma::vec sqrt_Sigma = arma::sqrt(Sigma);
    double mean_rmse_X  =  std::sqrt(arma::accu(arma::pow(Sigma, 2))) / M;

    new_X =  arma::diagmat(new_D_w_sqrt) * new_X * arma::diagmat(1 / sqrt_Sigma);
    new_Omega =  arma::diagmat(1 / sqrt_Sigma) *  new_Omega  * arma::diagmat(new_D_w_sqrt);
    arma::mat der_X, der_reg;
    arma::mat hinge_term_H, hinge_term_W;
    arma::mat tmp_X, tmp_Omega;
    double shrink_limit = 500;
    double limit_for_learning_rate = 1e-15;
    double current_learning_rate = coef_der_X;
    //double mean_norm_solution_X;
    double previous_error_value;
    double error_difference;



    

//    Rcpp::Rcout << "Start X"  << std::endl;
//    Rcpp::Rcout << new_X  << std::endl;
//    Rcpp::Rcout << "Start Omega"  << std::endl;
//    Rcpp::Rcout << new_Omega  << std::endl;

    // Start initial inverse search
  //  Rcpp::Rcout << "Check initial inverse matrix properties"  << std::endl;
    tmp_Omega = arma::pinv(new_X);
    if (arma::any( tmp_Omega.row(0) <= 0)) {
       Rcpp::Rcout << "Couldn't find good initial inverse of X provided\n"  << std::endl;
       Rcpp::Rcout << "Try with Omega \n"  << std::endl;
       tmp_X = arma::pinv(new_Omega);
       if (arma::any( tmp_X.col(0) <= 0)) {
            Rcpp::Rcout << "Couldn't find good initial inverse of Omega provided\n"  << std::endl;
            Rcpp::stop("!!Start with different initialization or ensure X and Omega are inverse!! (try `random_invertible`)");
        }
        else {
            new_X = tmp_X;
        }
    }
    else {
        new_Omega = tmp_Omega;
    }

//    Rcpp::Rcout << "Initial X\n"  << std::endl;
//    Rcpp::Rcout << new_X << std::endl;
//    Rcpp::Rcout << "Initial Omega\n"  << std::endl;
//    Rcpp::Rcout << new_Omega << std::endl;
//    Rcpp::Rcout << "Initial D\n"  << std::endl;
//    Rcpp::Rcout << new_D_w << std::endl;
//    Rcpp::Rcout << "Initial sqrt D\n"  << std::endl;
//    Rcpp::Rcout << new_D_w_sqrt << std::endl;
//
//    Rcpp::Rcout << "X and Omega are acceptable. Continue with the optimization\n"  << std::endl;

    // here we assume X and Omega are inverse of each other and positive as needed
    int itr_ = 0;
    while ((itr_ < iterations) & (current_learning_rate > limit_for_learning_rate)) {
        //mean_norm_solution_X = arma::mean(arma::vecnorm(new_X, 2, 1));

        hinge_term_H = l1_hinge_der_proportions_C__(new_X  * arma::diagmat(sqrt_Sigma)  * R, R) * arma::diagmat(sqrt_Sigma);
        hinge_term_W = (-new_Omega.t())  * arma::diagmat(sqrt_Sigma) * l1_hinge_der_basis_C__(S.t() * arma::diagmat(sqrt_Sigma) * new_Omega, S) * (new_Omega.t());
      
        der_X =  coef_hinge_H * hinge_term_H;
        der_X += coef_hinge_W * hinge_term_W;
        der_reg = reg_X * 2 * new_X; //regularization for X
        // Regularization here is advised but not mandatory since X and Omega regularize each other.
        der_reg +=  reg_Omega * (-new_Omega.t()) * 2 * new_Omega * (new_Omega.t()); //regularization for Omega
        der_X = correctByNorm(der_X); 
        der_reg = correctByNorm(der_reg); 
        der_X = der_X + der_reg;
        //der_X = correctByNorm(der_X)  * mean_radius_X; // arma::diagmat(new_D_w_sqrt)  * arma::diagmat(1 / sqrt_Sigma)  * mean_radius_X;
        der_X.each_row() %= new_D_w_sqrt.t() * mean_radius_X;
        der_X.each_col() /= sqrt_Sigma;

        tmp_X = (new_X - current_learning_rate * der_X); // estimate new X given derivative

        
        // Check if first column of X is all-positive
        if (arma::any(tmp_X.col(0) <= 0)) {
            for (int c=0; c < cell_types; c++) {
                double matrix_value =  tmp_X(c,0);
                 if (matrix_value <= 0) {
                   int shrink_iteration = 0;
                   while((matrix_value <= 0) & (shrink_iteration < shrink_limit)) {
                    der_X.row(c) /=  2;
                    tmp_X = (new_X - current_learning_rate * der_X);
                    matrix_value =  tmp_X(c,0);
                    shrink_iteration++;
                   }
                 }
            }
            if  (arma::any( tmp_X.col(0) <= 0)) {
                Rcpp::Rcout << "Any gradient step gives bad X, probably X was bad before\n"  << std::endl;
            }
        }

        tmp_Omega = arma::pinv(tmp_X);
        // Check if first row of Omega is all positive
        if (arma::any( tmp_Omega.row(0) <= 0)) {
            for (int c=0; c < cell_types; c++) {
                double matrix_value =  tmp_Omega(0,c);
                if (matrix_value <= 0) {
                    int shrink_iteration = 0;
                    while((matrix_value <= 0)& (shrink_iteration < shrink_limit)) {
                    der_X /=  2;
                    der_X.row(c) *= 2;
                    tmp_X = (new_X - current_learning_rate * der_X);
                    tmp_Omega = arma::pinv(tmp_X);
                    matrix_value =  tmp_Omega(0,c);
                    shrink_iteration++;
                   }
                if (shrink_iteration != shrink_limit) {
                    // if we were able to find the solution. accept these new X and Omega
                    // Do nothing its ok
                    } else {
                        Rcpp::Rcout << "Iteration \n"<<  itr_ << "\n" << std::endl;
                        Rcpp::Rcout << "Couldn't find good inverse X for the row " <<  c << ", reject X update for this row \n"  << std::endl;
                        arma::rowvec only_good_row  = der_X.row(c);
                        der_X.fill(0);
                        der_X.row(c) = only_good_row;
                    }
                }
            }
            new_Omega = tmp_Omega;
            new_X = tmp_X;
        } else {
            new_Omega = tmp_Omega;
            new_X = tmp_X;
        }
        // optional track X dtilda statistics
        // points_statistics_X_dtilda_uncorrected.row(itr_) = new_X.as_row();
        // points_statistics_Omega_dtilda_uncorrected.row(itr_) = new_Omega.as_row();

        std::tie(new_X, new_Omega, new_D_w_sqrt) = ensure_D_integrity_c(new_X, new_Omega, sqrt_Sigma, N, M);
        final_X = arma::diagmat(1/new_D_w_sqrt) * new_X * arma::diagmat(sqrt_Sigma);
        final_Omega = arma::diagmat(sqrt_Sigma)* new_Omega * arma::diagmat(1/new_D_w_sqrt);

        //// Theoretically this should not happen. But we keep this code for now to ensure everything is good.
    //    for (int c=0; c < cell_types; c++) {
    //        double col_omega_norm = arma::norm(final_Omega.col(c).subvec(1, cell_types - 1), 2);
    //        double row_x_norm = arma::norm(final_X.row(c).subvec(1, cell_types - 1), 2);
    //        if (col_omega_norm > solution_balancing_threshold * mean_radius_Omega) {
    //            Rcpp::Rcout << "Looks like Omega points are way far away after inverse of X. \n"  << std::endl;
    //            Rcpp::Rcout << "We will balance solution by moving some magnitude from Omega to X. \n"  << std::endl;
    //            double ratio_x = row_x_norm / mean_radius_X;
    //            double ratio_omega = col_omega_norm / mean_radius_Omega;
    //            // stretch the space for all elements
    //            double multiplier_x = sqrt(ratio_omega)/sqrt(ratio_x);
    //            double multiplier_omega =  sqrt(ratio_x)/sqrt(ratio_omega);
    //            Rcpp::Rcout << " X cols multiplied by " << multiplier_x << std::endl;
    //            Rcpp::Rcout << " Omega cols multiplied by " << multiplier_omega << std::endl;

    //             for (int column=1; column < cell_types; column++) {
    //                 std::tie(temporary_new_X, temporary_new_Omega, temporary_new_D_w_sqrt) = ensure_D_integrity_c(new_X, new_Omega, sqrt_Sigma, N, M);
    //                 final_X = arma::diagmat(1/new_D_w_sqrt) * new_X * arma::diagmat(sqrt_Sigma);
    //                 final_Omega = arma::diagmat(sqrt_Sigma)* new_Omega * arma::diagmat(1/new_D_w_sqrt);
    //                 new_X.col(column) *= multiplier_x;
    //                 new_Omega.row(column) *= multiplier_omega;
    //             }
    //             Rcpp::Rcout << "Finish norm correction. \n"  << std::endl;
    //         }

    //     }
    //     // Correct X and Omega to have corresponding first row/column and calculate D
    //     std::tie(new_X, new_Omega, new_D_w_sqrt) = ensure_D_integrity_c(new_X, new_Omega, sqrt_Sigma, N, M);
           // result X and omega
    //     final_X = arma::diagmat(1/new_D_w_sqrt) * new_X * arma::diagmat(sqrt_Sigma);
    //     final_Omega = arma::diagmat(sqrt_Sigma)* new_Omega * arma::diagmat(1/new_D_w_sqrt);
        new_D_w = arma::pow(new_D_w_sqrt, 2);
        new_D_h = new_D_w * (N / M);
        double sum_ = accu(new_D_w) / M;
        arma::uword neg_props = getNegative(final_X * R);
        arma::uword neg_basis = getNegative(S.t() * final_Omega);
        
        Rcpp::List current_errors = calcErrors(final_X,
                                               final_Omega,
                                               new_D_w,
                                               new_D_h,
                                               SVRt,
                                               R,
                                               S,
                                               coef_hinge_H,
                                               coef_hinge_W);
        double current_error_value = current_errors["total_error"];
        error_difference = std::abs(previous_error_value - current_error_value);
        if (error_difference < convergence_tol) {
            current_learning_rate = current_learning_rate / 2;
        }
        Rcpp::Rcout << "Error difference: "<< error_difference << std::endl;
        previous_error_value = current_error_value;
                                               
        errors_statistics.row(itr_) = arma::rowvec{current_errors["deconv_error"],
                                                   current_errors["lambda_error"],
                                                   current_errors["beta_error"],
                                                   current_errors["D_h_error"],
                                                   current_errors["D_w_error"],
                                                   current_errors["total_error"],
                                                   static_cast<double>(neg_props),
                                                   static_cast<double>(neg_basis),
                                                   sum_,
                                                   current_errors["average_norm"]};
        
        //points_statistics_X_dtilda_corrected.row(itr_) = new_X.as_row();
        //points_statistics_Omega_dtilda_corrected.row(itr_) = new_Omega.as_row();
        points_statistics_X.row(itr_) = final_X.as_row();
        points_statistics_Omega.row(itr_) = final_Omega.as_row();
        points_statistics_Dw.row(itr_) = new_D_w.as_row();

        itr_++;
    }
    if (itr_ < iterations) {
        points_statistics_X.resize(itr_, points_statistics_X.n_cols);
        points_statistics_Omega.resize(itr_, points_statistics_Omega.n_cols);
        points_statistics_Dw.resize(itr_, points_statistics_Dw.n_cols);
        errors_statistics.resize(itr_, errors_statistics.n_cols);
    }



    return Rcpp::List::create(Rcpp::Named("new_X") = final_X,
                              Rcpp::Named("new_Omega") = final_Omega,
                              Rcpp::Named("new_D_w") = new_D_w,
                              Rcpp::Named("new_D_h") = new_D_h,
                              Rcpp::Named("errors_statistics") = errors_statistics,
                              Rcpp::Named("points_statistics_X") = points_statistics_X,
                              Rcpp::Named("points_statistics_Omega") = points_statistics_Omega,
                              Rcpp::Named("points_statistics_Dw") = points_statistics_Dw
                            //   ,Rcpp::Named("points_statistics_Omega_dtilda_uncorrected") = points_statistics_Omega_dtilda_uncorrected,
                            //   Rcpp::Named("points_statistics_Omega_dtilda") = points_statistics_Omega_dtilda_corrected,
                            //   Rcpp::Named("points_statistics_X_dtilda_uncorrected") = points_statistics_X_dtilda_uncorrected,
                            //   Rcpp::Named("points_statistics_X_dtilda") = points_statistics_X_dtilda_corrected
                              );
}
