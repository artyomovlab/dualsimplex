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
                             const double thresh) {
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
    arma::mat old_X, old_Omega;
    arma::mat old_X_X, old_Omega_X;
    arma::mat new_X_after_X, new_Omega_after_X;
    for (int itr_ = 0; itr_ < iterations; itr_++) {
        // derivative X
        //der_X = -2 * (new_Omega.t() * (SVRt - new_Omega * new_X));
        //der_X +=  coef_hinge_H * hinge_der_proportions_C__(new_X  * R, R);

        der_X =  coef_hinge_H * hinge_der_proportions_C__(new_X  * arma::diagmat(sqrt_Sigma)  * R, R) * arma::diagmat(1 / sqrt_Sigma);

        der_X = correctByNorm(der_X) * mean_radius_X;


        // Update X
        old_X_X = new_X;
        old_Omega_X = new_Omega;
        new_X = new_X - coef_der_X * der_X;

        // replace negative values with very small number
        if (any( new_X.col(0) < 0)) {
           Rcpp::Rcout << "Derrivative caused negative for X \n"  << std::endl;
           Rcpp::Rcout << "---New X---" << std::endl;
           Rcpp::Rcout << new_X << "\n";
           new_X.col(0).clamp(arma::datum::eps,arma::datum::inf);
           Rcpp::Rcout << "---Clamped X---" << std::endl;
           Rcpp::Rcout << new_X << "\n";
        }



        try {
            new_Omega = arma::pinv(new_X);
            // replace negative values with very small number
           if (any( new_Omega.row(0) < 0)) {
               Rcpp::Rcout << "Inverse caused negative for Omega \n"  << std::endl;
               Rcpp::Rcout << "---New Omega---" << std::endl;
               Rcpp::Rcout << new_Omega << "\n";
               new_Omega.row(0).clamp(arma::datum::eps,arma::datum::inf);
               Rcpp::Rcout << "---Clamped Omega---" << std::endl;
               Rcpp::Rcout << new_Omega << "\n";
           }


        }
        catch (const std::runtime_error& e)
        {
         Rcpp::Rcout << "Error in inverse \n" << e.what() << std::endl;
         Rcpp::Rcout << "---Old X---"  << std::endl;
         Rcpp::Rcout << old_X_X << "\n";
         Rcpp::Rcout << "---Derrivative ---" << std::endl;
         Rcpp::Rcout << der_X << "\n";
         Rcpp::Rcout << "---New X---" << std::endl;
         Rcpp::Rcout << new_X << "\n";
         Rcpp::Rcout << "---Old Omega---" << std::endl;
         Rcpp::Rcout << old_Omega_X << "\n";
         Rcpp::Rcout << "---Current Dx---" << std::endl;
         Rcpp::Rcout << new_D_w_x << "\n";
         Rcpp::Rcout << "---Current Domega---"<< std::endl;
         Rcpp::Rcout << new_D_w_omega << "\n";
        }
        new_X_after_X = new_X;
        new_Omega_after_X = new_Omega;

        new_D_w_x_sqrt =  new_X.col(0) * sqrt_Sigma.at(0) * sqrt(N);
        new_D_w_x = arma::pow(new_D_w_x_sqrt, 2);

        new_D_w_omega_sqrt =  new_Omega.row(0).as_col() * sqrt_Sigma.at(0) * sqrt(M);
        new_D_w_omega = arma::pow(new_D_w_omega_sqrt, 2);

        new_D_w = new_D_w_x_sqrt % new_D_w_omega_sqrt;

        new_D_w_sqrt = arma::sqrt(new_D_w);
        new_X = arma::diagmat(1/new_D_w_x_sqrt)* arma::diagmat(new_D_w_sqrt) * new_X;
        new_Omega = new_Omega * arma::diagmat(1/new_D_w_omega_sqrt) * arma::diagmat(new_D_w_sqrt);
        new_D_w_x_sqrt = new_D_w_sqrt;
        new_D_w_omega_sqrt = new_D_w_sqrt;



        // derivative Omega

        der_Omega = coef_hinge_W * arma::diagmat(1 / sqrt_Sigma) * alternative_hinge_der_basis_C__(S.t() * arma::diagmat(sqrt_Sigma) * new_Omega, S);


        der_Omega = correctByNorm(der_Omega) * mean_radius_Omega;

        old_X = new_X;
        old_Omega = new_Omega;


        new_Omega = new_Omega - coef_der_Omega * der_Omega;
        if (any( new_Omega.row(0) < 0)) {
               Rcpp::Rcout << "Derrivative caused negative for Omega \n"  << std::endl;
               Rcpp::Rcout << "---New Omega---" << std::endl;
               Rcpp::Rcout << new_Omega << "\n";
               new_Omega.row(0).clamp(arma::datum::eps,arma::datum::inf);
               Rcpp::Rcout << "---Clamped Omega---" << std::endl;
               Rcpp::Rcout << new_Omega << "\n";
        }


        try {
        new_X = arma::pinv(new_Omega);
                if (any( new_X.col(0) < 0)) {
           Rcpp::Rcout << "Inverse caused negative for X \n"  << std::endl;
           Rcpp::Rcout << "---New X---" << std::endl;
           Rcpp::Rcout << new_X << "\n";
           new_X.col(0).clamp(arma::datum::eps,arma::datum::inf);
           Rcpp::Rcout << "---Clamped X---" << std::endl;
           Rcpp::Rcout << new_X << "\n";
        }

          }
        catch (const std::runtime_error& e)
        {
         Rcpp::Rcout << "Error in inverse \n" << e.what() << std::endl;


          Rcpp::Rcout << "---Old X_X---"  << std::endl;
         Rcpp::Rcout << old_X_X << "\n";
         Rcpp::Rcout << "---Old Omega_X---"  << std::endl;
         Rcpp::Rcout << old_Omega_X << "\n";
         Rcpp::Rcout << "---New X_X---"  << std::endl;
         Rcpp::Rcout << new_X_after_X << "\n";
         Rcpp::Rcout << "---New Omega X---"  << std::endl;
         Rcpp::Rcout << new_Omega_after_X << "\n";
         Rcpp::Rcout << "---Old X---"  << std::endl;
         Rcpp::Rcout << old_X << "\n";
         Rcpp::Rcout << "---Old Omega---"  << std::endl;
         Rcpp::Rcout << old_Omega << "\n";
         Rcpp::Rcout << "---Derrivative X---" << std::endl;
         Rcpp::Rcout << der_X << "\n";
         Rcpp::Rcout << "---Derrivative Omega---" << std::endl;
         Rcpp::Rcout << der_Omega << "\n";
         Rcpp::Rcout << "---New X---" << std::endl;
         Rcpp::Rcout << new_X << "\n";
         Rcpp::Rcout << "---New Omega---" << std::endl;
         Rcpp::Rcout << new_Omega << "\n";
         Rcpp::Rcout << "---Current Dx---" << std::endl;
         Rcpp::Rcout << new_D_w_x << "\n";
         Rcpp::Rcout << "---Current Domega---"<< std::endl;
         Rcpp::Rcout << new_D_w_omega << "\n";
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
        arma::uword neg_props = getNegative(new_X * R);
        arma::uword neg_basis = getNegative(S.t() * new_Omega);
        double sum_ = accu(new_D_w) / M;

        // result X and omega
        final_X = arma::diagmat(1/new_D_w_x_sqrt) * new_X * arma::diagmat(sqrt_Sigma);

        final_Omega = arma::diagmat(sqrt_Sigma)* new_Omega * arma::diagmat(1/new_D_w_omega_sqrt);


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
