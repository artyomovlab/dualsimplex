#include "optimization_utils.h"
#include "matrix_utils.h"
#include "nnls.h"
#include <tuple>
#include "optimization_positivity.h"

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

Rcpp::List optimize_positivity(const arma::mat& X,
                             const arma::mat& Omega,
                             const arma::mat& D_w,
                             const arma::mat& SVRt,
                             const arma::mat& R,
                             const arma::mat& S,
                             const double coef_der_X,
                             double coef_hinge_H,
                             double coef_hinge_W,
                             const int cell_types,
                             const double N,
                             const double M,
                             const int iterations,
                             double total_regularization_weight,
                             const double reg_X,
                             const double reg_Omega,
                             const double convergence_tol,
                             const int stop_criteria_window,
                             const bool debug_stats
                             ) {
    arma::mat errors_statistics(iterations + 1, 21, arma::fill::zeros);

    arma::mat points_statistics_X(iterations + 1, cell_types * cell_types, arma::fill::zeros);
    arma::mat points_statistics_Omega(iterations + 1, cell_types * cell_types, arma::fill::zeros);
    // arma::mat points_statistics_X_dtilda_uncorrected(iterations, cell_types * cell_types, arma::fill::zeros);
    // arma::mat points_statistics_Omega_dtilda_uncorrected(iterations, cell_types * cell_types, arma::fill::zeros);
    //arma::mat points_statistics_X_dtilda_corrected(iterations, cell_types * cell_types, arma::fill::zeros);
    //arma::mat points_statistics_Omega_dtilda_corrected(iterations, cell_types * cell_types, arma::fill::zeros);
    arma::mat points_statistics_Dw(iterations + 1, cell_types, arma::fill::zeros);

    arma::mat new_X = X;
    arma::mat new_Omega = Omega;
    arma::mat final_X = X;
    arma::mat final_Omega = Omega;
    arma::mat new_D_w = D_w;
    arma::mat new_D_w_sqrt = arma::sqrt(new_D_w);
    arma::mat new_D_h = new_D_w * (N / M);

    arma::vec Sigma = arma::diagvec(SVRt);
    arma::vec sqrt_Sigma = arma::sqrt(Sigma);

    new_X =  arma::diagmat(new_D_w_sqrt) * new_X * arma::diagmat(1 / sqrt_Sigma);
    new_Omega =  arma::diagmat(1 / sqrt_Sigma) *  new_Omega  * arma::diagmat(new_D_w_sqrt);
    arma::mat der_X, der_reg;
    arma::mat hinge_term_H, hinge_term_W, reg_X_term, reg_Omega_term;
    arma::mat tmp_X, tmp_Omega;
    // make all coefficients sum to 1 for simplicity
    double coef_sum = coef_hinge_H + coef_hinge_W + total_regularization_weight;
    coef_hinge_H = coef_hinge_H / coef_sum;
    coef_hinge_W = coef_hinge_W / coef_sum;
    total_regularization_weight = total_regularization_weight / coef_sum; 
    double shrink_limit = 500;
    double current_learning_rate = coef_der_X;
    double best_error_value = 10000;
    int best_error_iteration = 0;
    double current_error_value;
    double average_gradient_norm = 0;
    double average_hinge_H_gradient_norm = 0;
    double average_hinge_W_gradient_norm = 0;
    double average_hinge_reg_X_gradient_norm = 0;
    double average_hinge_reg_Omega_gradient_norm = 0;
              
//    Rcpp::Rcout << "Start X"  << std::endl;
//    Rcpp::Rcout << new_X  << std::endl;
//    Rcpp::Rcout << "Start Omega"  << std::endl;
//    Rcpp::Rcout << new_Omega  << std::endl;

    // Start initial inverse search
  //  Rcpp::Rcout << "Check initial inverse matrix properties"  << std::endl;
                        
    tmp_Omega = arma::pinv(new_X);
    if (arma::any( tmp_Omega.row(0) <= 0)) {
       spdl::warn("Couldn't find good initial inverse of X provided. Will try with Omega");
       tmp_X = arma::pinv(new_Omega);
       if (arma::any( tmp_X.col(0) <= 0)) {
            spdl::warn("Couldn't find good initial inverse of Omega provided");
            Rcpp::stop("!!Start with different initialization or ensure X and Omega are inverse!! (try `random_invertible`)");
        }
        else {
            new_X = tmp_X;
        }
    }
    else {
        new_Omega = tmp_Omega;
    }

    // here we assume X and Omega are inverse of each other and positive as needed
    spdl::warn("start with X");
    Rcpp::Rcout << new_X  << std::endl;
    spdl::warn("start with Omega");
    Rcpp::Rcout << new_Omega  << std::endl;
    spdl::warn("start with D");
    Rcpp::Rcout << new_D_w  << std::endl;

    int itr_ = 0;
    while ((itr_ < iterations + 1) & (current_learning_rate > convergence_tol)) {
        if (itr_ > 0) { // in order to save initial errors, skip first step.
        hinge_term_H = l1_hinge_der_proportions_C__(new_X  * arma::diagmat(sqrt_Sigma)  * R, R) * arma::diagmat(sqrt_Sigma);
        hinge_term_W = (-new_Omega.t())  * arma::diagmat(sqrt_Sigma) * l1_hinge_der_basis_C__(S.t() * arma::diagmat(sqrt_Sigma) * new_Omega, S) * (new_Omega.t());
        der_X =  coef_hinge_H * hinge_term_H;
        der_X += coef_hinge_W * hinge_term_W;


        reg_X_term = 2 * new_X;
        reg_Omega_term = (-new_Omega.t()) * 2 * new_Omega * (new_Omega.t());
        der_reg = reg_X * reg_X_term +  reg_Omega * reg_Omega_term ; //regularization for X
        
        // Regularization here is advised but not mandatory since X and Omega regularize each other.
        der_X = der_X + total_regularization_weight * der_reg;
        spdl::warn("derrivative X");
        Rcpp::Rcout << der_X  << std::endl;
        if (debug_stats) {
            average_gradient_norm = arma::mean(arma::vecnorm(der_X, 2, 1));
            average_hinge_H_gradient_norm = arma::mean(arma::vecnorm(hinge_term_H, 2, 1));
            average_hinge_W_gradient_norm = arma::mean(arma::vecnorm(hinge_term_W, 2, 1));
            average_hinge_reg_X_gradient_norm = arma::mean(arma::vecnorm(reg_X_term, 2, 1));
            average_hinge_reg_Omega_gradient_norm = arma::mean(arma::vecnorm(reg_Omega_term, 2, 1));
        }
        tmp_X = (new_X - current_learning_rate * der_X); // estimate new X given derivative
        spdl::warn("candidate X");
        Rcpp::Rcout << tmp_X  << std::endl;
        // Ensure if first column of X is all-positive
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
                spdl::warn("Any gradient step gives bad X, probably X was bad before");
            }
        }
        spdl::warn("corrected X");
        Rcpp::Rcout << tmp_X  << std::endl;

        tmp_Omega = arma::pinv(tmp_X);
        spdl::warn("candidate Omega");
        Rcpp::Rcout << tmp_Omega  << std::endl;
        // Ensure if first row of Omega is all positive
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
                        spdl::warn("Iteration {} Couldn't find good inverse X for the row {}, reject X update for this row", itr_, c);
                        arma::rowvec only_good_row  = der_X.row(c);
                        der_X.fill(0);
                        der_X.row(c) = only_good_row;
                        tmp_X = (new_X - current_learning_rate * der_X);
                    }
                }
            }
            spdl::warn("corrected X");
            Rcpp::Rcout << tmp_X  << std::endl;
            spdl::warn("corrected Omega");
            Rcpp::Rcout << tmp_Omega  << std::endl;
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

        new_D_w = arma::pow(new_D_w_sqrt, 2);
        new_D_h = new_D_w * (N / M);
        }
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

        current_error_value = current_errors["total_error"];
        if (current_error_value < best_error_value) {
            best_error_iteration = itr_;
            best_error_value = current_error_value;
        }
        if (itr_ - best_error_iteration > stop_criteria_window) {
            // looks like best solution was not updated for stop_criteria_window iterations. reducing step size.
            current_learning_rate = current_learning_rate / 2;
            best_error_iteration = itr_; // reset iteration counter
        }
                                               
        errors_statistics.row(itr_) = arma::rowvec{current_errors["deconv_error"],            //1
                                                   current_errors["lambda_error"], //2
                                                   current_errors["beta_error"], //3
                                                   current_errors["D_h_error"], //4
                                                   current_errors["D_w_error"], //5
                                                   current_errors["total_error"], //6
                                                   current_errors["scaled_total_error"], //7
                                                   static_cast<double>(neg_props), //8
                                                   static_cast<double>(neg_basis), //9
                                                   sum_, //10
                                                   current_errors["average_norm"], //11
                                                   current_learning_rate, //12
                                                   average_gradient_norm, //13
                                                   average_hinge_H_gradient_norm, //14
                                                   average_hinge_W_gradient_norm, //15
                                                   average_hinge_reg_X_gradient_norm, //16
                                                   average_hinge_reg_Omega_gradient_norm, //17
                                                   best_error_value, //18
                                                   static_cast<double>(best_error_iteration), //19,
                                                   current_errors["scaled_lambda_error"],            //20
                                                   current_errors["scaled_beta_error"] //21
                                                };
        
        //points_statistics_X_dtilda_corrected.row(itr_) = new_X.as_row();
        //points_statistics_Omega_dtilda_corrected.row(itr_) = new_Omega.as_row();
        points_statistics_X.row(itr_) = final_X.as_row();
        points_statistics_Omega.row(itr_) = final_Omega.as_row();
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
