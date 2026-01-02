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
    // derivative should be the same as for X but W is transposed
    arma::mat res = l1_hinge_der_proportions_C__(W.t(), S);
    return res.t();
}


arma::mat squared_hinge_der_basis_C__(const arma::mat& W, const arma::mat& S) {
    // derrivative should be the same as for X but W is transposed
    arma::mat res = squared_hinge_der_proportions_C__(W.t(), S);
    return res.t();
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
                             const double solution_balancing_threshold,
                             const double coef_norm) {
    arma::mat errors_statistics(iterations, 10, arma::fill::zeros);
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
    arma::mat hinge_term_H, hinge_term_W;
    arma::mat tmp_X, tmp_Omega;
    double shrink_limit = 500;
    double mean_norm_solution_X;

    Rcpp::Rcout << "Start X"  << std::endl;
    Rcpp::Rcout << new_X  << std::endl;
    Rcpp::Rcout << "Start Omega"  << std::endl;
    Rcpp::Rcout << new_Omega  << std::endl;

    // Start initial inverse search
    Rcpp::Rcout << "Start initial inverse search"  << std::endl;
    tmp_Omega = arma::pinv(new_X);
    if (any( tmp_Omega.row(0) <= 0)) {
       Rcpp::Rcout << "Couldn't find good initial inverse of X \n"  << std::endl;
       Rcpp::Rcout << "Try with Omega \n"  << std::endl;
       tmp_X = arma::pinv(new_Omega);
       if (any( tmp_X.col(0) <= 0)) {
            Rcpp::Rcout << "Couldn't find good initial inverse of Omega\n"  << std::endl;
            Rcpp::Rcout << "!!Start with different initialization!!\n"  << std::endl;
    }
    else {
        new_X = tmp_X;
    }}
    else {
        new_Omega = tmp_Omega;
    }

    // here we assume X and Omega are inverse of each other and positive as needed
    for (int itr_ = 0; itr_ < iterations; itr_++) {
        hinge_term_H = l1_hinge_der_proportions_C__(new_X  * arma::diagmat(sqrt_Sigma)  * R, R) * arma::diagmat(sqrt_Sigma);
//        hinge_term_H = correctByNorm(hinge_term_H);
        hinge_term_W = (-new_Omega.t())  * arma::diagmat(sqrt_Sigma) * l1_hinge_der_basis_C__(S.t() * arma::diagmat(sqrt_Sigma) * new_Omega, S) * (new_Omega.t());
//        hinge_term_W = correctByNorm(hinge_term_W);

        der_X =  coef_hinge_H * hinge_term_H;
        der_X += coef_hinge_W * hinge_term_W;
        der_X += 2 * new_X; //regularization for X
        der_X += (-new_Omega.t()) * 2 * new_Omega * (new_Omega.t()); //regularization for Omega


        mean_norm_solution_X = arma::mean(arma::vecnorm(new_X, 2, 1));
        der_X = correctByNorm(der_X) * mean_norm_solution_X; // arma::diagmat(new_D_w_sqrt)  * arma::diagmat(1 / sqrt_Sigma)  * mean_radius_X;

        tmp_X = (new_X - coef_der_X * der_X); // estimate new X given derivative
        // Check if first column of X is all-positive
        if (any( tmp_X.col(0) <= 0)) {
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
        }
       tmp_Omega = arma::pinv(tmp_X);
       // Check if first row of Omega is all positive
        if (any( tmp_Omega.row(0) <= 0)) {
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

       std::tie(new_X, new_Omega, new_D_w_sqrt) = ensure_D_integrity_c(new_X, new_Omega, sqrt_Sigma, N, M);
       final_X = arma::diagmat(1/new_D_w_sqrt) * new_X * arma::diagmat(sqrt_Sigma);
       final_Omega = arma::diagmat(sqrt_Sigma)* new_Omega * arma::diagmat(1/new_D_w_sqrt);

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
                    new_X.col(column) *= multiplier_x;
                    new_Omega.row(column) *= multiplier_omega;
                }
                Rcpp::Rcout << "Finish norm correction. \n"  << std::endl;
            }

        }
// Correct X and Omega to have corresponding first row/column and derrive D
        std::tie(new_X, new_Omega, new_D_w_sqrt) = ensure_D_integrity_c(new_X, new_Omega, sqrt_Sigma, N, M);

//       Rcpp::Rcout << " After final correction X"  << std::endl;
//       Rcpp::Rcout << new_X << std::endl;
//       Rcpp::Rcout << " After final correction Omega"  << std::endl;
//       Rcpp::Rcout << new_Omega << std::endl;


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
                                               coef_pos_D_w,
                                               coef_norm);

        errors_statistics.row(itr_) = arma::rowvec{current_errors["deconv_error"],
                                                   current_errors["squared_lambda_error"],
                                                   current_errors["squared_beta_error"],
                                                   current_errors["D_h_error"],
                                                   current_errors["D_w_error"],
                                                   current_errors["total_error"],
                                                   static_cast<double>(neg_props),
                                                   static_cast<double>(neg_basis),
                                                   sum_,
                                                   current_errors["average_norm"]};


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
