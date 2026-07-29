#include "optimization_utils.h"
#include "matrix_utils.h"
#include "nnls.h"
#include <tuple>
#include "optimization_positivity.h"
#include "optimization_dsa.h"

struct HingeStats {
    double hinge_sum{0.0};
    arma::uword neg_count{0};
};

inline HingeStats compute_hinge_stats(const arma::mat& M) {
    arma::uvec neg_indices = arma::find(M < 0);
    if (neg_indices.is_empty()) {
        return {0.0, 0};
    }
    double hinge_sum = -arma::accu(M.elem(neg_indices));
    return {hinge_sum, neg_indices.n_elem};
}

// Strategy Calculators Implementation
std::vector<Metric> DeconvErrorCalculator::calculate(const OptimizationState& state) const {
    arma::mat D_w_diag = arma::diagmat(state.D_w);
    double deconv_error = std::pow(arma::norm(state.SVRt - state.Omega * D_w_diag * state.X, "fro"), 2.0);
    return {{"deconv_error", deconv_error}};
}

std::vector<Metric> HingeProportionsErrorCalculator::calculate(const OptimizationState& state) const {
    arma::mat H = state.X * state.R;
    HingeStats stats = compute_hinge_stats(H);
    double lambda_error = state.coef_hinge_H * stats.hinge_sum;
    double scaled_lambda_error = (state.R.n_cols > 0) ? (lambda_error / state.R.n_cols) : lambda_error;
    return {
        {"lambda_error", lambda_error},
        {"scaled_lambda_error", scaled_lambda_error},
        {"neg_props", static_cast<double>(stats.neg_count)}
    };
}

std::vector<Metric> HingeBasisErrorCalculator::calculate(const OptimizationState& state) const {
    arma::mat W = state.S.t() * state.Omega;
    HingeStats stats = compute_hinge_stats(W);
    double beta_error = state.coef_hinge_W * stats.hinge_sum;
    double scaled_beta_error = (state.S.n_cols > 0) ? (beta_error / state.S.n_cols) : beta_error;
    return {
        {"beta_error", beta_error},
        {"scaled_beta_error", scaled_beta_error},
        {"neg_basis", static_cast<double>(stats.neg_count)}
    };
}

std::vector<Metric> ScaleNormErrorCalculator::calculate(const OptimizationState& state) const {
    double average_norm_X = arma::mean(arma::vecnorm(state.X, 2, 1));
    double average_norm_Omega = arma::sum(arma::vecnorm(state.Omega, 2, 0));
    return {{"average_norm", average_norm_X + average_norm_Omega}};
}

Rcpp::List optimize_alignment(
    const arma::mat& X,
    const arma::mat& Omega,
    const arma::mat& D_w,
    const arma::mat& SVRt,
    const arma::mat& R,
    const arma::mat& S,
    const double coef_der_X,
    double coef_hinge_W,
    double coef_hinge_H,
    double coef_alignment,
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
    arma::mat points_statistics_X(iterations + 1, cell_types * cell_types, arma::fill::zeros);
    arma::mat points_statistics_Omega(iterations + 1, cell_types * cell_types, arma::fill::zeros);
    arma::mat points_statistics_Dw(iterations + 1, cell_types, arma::fill::zeros);

    CompositeErrorCalculator composite_calc;
    composite_calc.add_calculator(std::make_shared<DeconvErrorCalculator>());
    composite_calc.add_calculator(std::make_shared<HingeProportionsErrorCalculator>());
    composite_calc.add_calculator(std::make_shared<HingeBasisErrorCalculator>());
    composite_calc.add_calculator(std::make_shared<ScaleNormErrorCalculator>());

    RcppLogger rcpp_logger;
    std::map<std::string, std::string> metadata;
    metadata["algorithm"] = "optimize_alignment";
    metadata["iterations"] = std::to_string(iterations);
    metadata["coef_hinge_H"] = std::to_string(coef_hinge_H);
    metadata["coef_hinge_W"] = std::to_string(coef_hinge_W);
    metadata["coef_alignment"] = std::to_string(coef_alignment);
    rcpp_logger.log_metadata(metadata);

    arma::mat new_X = X;
    arma::mat new_Omega = Omega;
    arma::mat final_X = X;
    arma::mat final_Omega = Omega;
    arma::mat new_D_w = D_w;
    arma::mat new_D_w_sqrt = arma::sqrt(new_D_w);
    arma::mat new_D_h = new_D_w * (N / M);
    arma::mat Y;  // alignment difference
    arma::mat tmp_X, tmp_Omega;

    arma::vec Sigma = arma::diagvec(SVRt);  // this should be sigma_ss
    arma::vec sqrt_Sigma = arma::sqrt(Sigma);
    arma::vec sigma_fs = Sigma * N / M;

    new_X =  arma::diagmat(new_D_w_sqrt) * new_X * arma::diagmat(1 / sqrt_Sigma);
    new_Omega =  arma::diagmat(1 / sqrt_Sigma) *  new_Omega  * arma::diagmat(new_D_w_sqrt);

    // derevitives
    const int k = cell_types;
    arma::mat der_X(k, k, arma::fill::zeros);
    arma::mat der_reg(k, k, arma::fill::zeros);
    arma::mat hinge_term_H(k, k, arma::fill::zeros);
    arma::mat hinge_term_W(k, k, arma::fill::zeros);
    arma::mat reg_X_term(k, k, arma::fill::zeros);
    arma::mat reg_Omega_term(k, k, arma::fill::zeros);
    arma::mat der_dsa(k, k, arma::fill::zeros);

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
              

    // Start initial inverse search
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

    // TODO: log the error for initialized X and Omega, so we can start iteration from iter_ = 1


    // here we assume X and Omega are inverse of each other and positive as needed
    int itr_ = 0;
    while ((itr_ < iterations + 1) & (current_learning_rate > convergence_tol)) {
        if (itr_ > 0) { // in order to save initial errors, skip first step.
            // For data with only pure sample, we could use dual sipmlex alignment to find verticex of sample
            // simplex instead of negativity term.
            Y = new_Omega - arma::diagmat(sigma_fs) * new_X.t();  // alignment differenct X_dtilda^-1 - sigma_fs * X_dtilda^T
            der_dsa = -2 * new_X.t() * Y * new_X.t() - 2 * Y.t() * arma::diagmat(sigma_fs);
            der_X = coef_alignment * der_dsa;

            // negativity terms for gene and cell simplex. 
            // To save computation time, we only compute it if the coefficient is non-zero.
            if (0 < coef_hinge_W) {  // vertices of cell simplex
                hinge_term_W = (-new_Omega.t()) * arma::diagmat(sqrt_Sigma) * l1_hinge_der_basis_C__(S.t() * arma::diagmat(sqrt_Sigma) * new_Omega, S) * (new_Omega.t());
                der_X += coef_hinge_W * hinge_term_W;
            }
            if (0 < coef_hinge_H) {  // vertices of gene simplex
                hinge_term_H = l1_hinge_der_proportions_C__(new_X * arma::diagmat(sqrt_Sigma) * R, R) * arma::diagmat(sqrt_Sigma);
                der_X +=  coef_hinge_H * hinge_term_H;
            }


            //    // Regularization here is advised but not mandatory since X and Omega regularize each other.
            //     reg_X_term = 2 * new_X;
            //     reg_Omega_term = (-new_Omega.t()) * 2 * new_Omega * (new_Omega.t());
            //     der_reg = reg_X * reg_X_term +  reg_Omega * reg_Omega_term ; //regularization for X 
            //     der_X = der_X + total_regularization_weight * der_reg;

            if (debug_stats) {
                average_gradient_norm = arma::mean(arma::vecnorm(der_X, 2, 1));
                average_hinge_H_gradient_norm = arma::mean(arma::vecnorm(hinge_term_H, 2, 1));
                average_hinge_W_gradient_norm = arma::mean(arma::vecnorm(hinge_term_W, 2, 1));
                // average_hinge_reg_X_gradient_norm = arma::mean(arma::vecnorm(reg_X_term, 2, 1));
                // average_hinge_reg_Omega_gradient_norm = arma::mean(arma::vecnorm(reg_Omega_term, 2, 1));
            }
            tmp_X = (new_X - current_learning_rate * der_X); // estimate new X given derivative

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

            tmp_Omega = arma::pinv(tmp_X);
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
                new_Omega = tmp_Omega;
                new_X = tmp_X;
            } else {
                new_Omega = tmp_Omega;
                new_X = tmp_X;
            }

            std::tie(new_X, new_Omega, new_D_w_sqrt) = ensure_D_integrity_c(new_X, new_Omega, sqrt_Sigma, N, M);
            final_X = arma::diagmat(1/new_D_w_sqrt) * new_X * arma::diagmat(sqrt_Sigma);
            final_Omega = arma::diagmat(sqrt_Sigma)* new_Omega * arma::diagmat(1/new_D_w_sqrt);

            new_D_w = arma::pow(new_D_w_sqrt, 2);
            new_D_h = new_D_w * (N / M);
        }
        
        double sum_ = accu(new_D_w) / M;
        
        OptimizationState current_state{
            final_X, final_Omega, new_D_w, new_D_h, SVRt, R, S, coef_hinge_H, coef_hinge_W
        };
        std::vector<Metric> metrics = composite_calc.calculate(current_state);

        // [CJLee] Extract metrics for control flow logic
        double deconv_err = 0.0, lambda_err = 0.0, beta_err = 0.0;
        double scaled_lambda_err = 0.0, scaled_beta_err = 0.0;
        for (const auto& m : metrics) {
            if (m.name == "deconv_error") deconv_err = m.value;
            else if (m.name == "lambda_error") lambda_err = m.value;
            else if (m.name == "beta_error") beta_err = m.value;
            else if (m.name == "scaled_lambda_error") scaled_lambda_err = m.value;
            else if (m.name == "scaled_beta_error") scaled_beta_err = m.value;
        }

        current_error_value = deconv_err + lambda_err + beta_err;
        double scaled_total_error = deconv_err + scaled_lambda_err + scaled_beta_err;

        if (current_error_value < best_error_value) {
            best_error_iteration = itr_;
            best_error_value = current_error_value;
        }
        if (itr_ - best_error_iteration > stop_criteria_window) {
            // looks like best solution was not updated for stop_criteria_window iterations. reducing step size.
            current_learning_rate = current_learning_rate / 2;
            best_error_iteration = itr_; // reset iteration counter
        }

        // [CJLee] Mostly for backward compatibility
        metrics.push_back({"total_error", current_error_value});
        metrics.push_back({"scaled_total_error", scaled_total_error});
        metrics.push_back({"learning_rate", current_learning_rate});
        metrics.push_back({"sum_", sum_});
        metrics.push_back({"average_gradient_norm", average_gradient_norm});
        metrics.push_back({"average_hinge_H_gradient_norm", average_hinge_H_gradient_norm});
        metrics.push_back({"average_hinge_W_gradient_norm", average_hinge_W_gradient_norm});
        metrics.push_back({"average_hinge_reg_X_gradient_norm", average_hinge_reg_X_gradient_norm});
        metrics.push_back({"average_hinge_reg_Omega_gradient_norm", average_hinge_reg_Omega_gradient_norm});
        metrics.push_back({"best_error_value", best_error_value});
        metrics.push_back({"best_error_iteration", static_cast<double>(best_error_iteration)});

        rcpp_logger.log_metrics(itr_, metrics);
        
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
    }

    arma::mat errors_statistics = rcpp_logger.to_legacy_matrix();

    return Rcpp::List::create(
        Rcpp::Named("new_X") = final_X,
        Rcpp::Named("new_Omega") = final_Omega,
        Rcpp::Named("new_D_w") = new_D_w,
        Rcpp::Named("new_D_h") = new_D_h,
        Rcpp::Named("errors_statistics") = errors_statistics,
        Rcpp::Named("points_statistics_X") = points_statistics_X,
        Rcpp::Named("points_statistics_Omega") = points_statistics_Omega,
        Rcpp::Named("points_statistics_Dw") = points_statistics_Dw,
        Rcpp::Named("history") = rcpp_logger.to_dataframe(),
        Rcpp::Named("metadata") = rcpp_logger.to_metadata_list()
    );
}




Rcpp::List optimize_alignment_hard(
    const arma::mat& X,
    const arma::mat& Omega,
    const arma::mat& D_w,
    const arma::mat& SVRt,
    const arma::mat& R,
    const arma::mat& S,
    const double coef_der_X,
    double coef_hinge_W,
    double coef_hinge_H,
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
              
    // Start initial inverse search
    // TODO (cjlee): re-write how tmp_Omega is found. With DSA, we don't need to calculate inverse
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
        if (debug_stats) {
            average_gradient_norm = arma::mean(arma::vecnorm(der_X, 2, 1));
            average_hinge_H_gradient_norm = arma::mean(arma::vecnorm(hinge_term_H, 2, 1));
            average_hinge_W_gradient_norm = arma::mean(arma::vecnorm(hinge_term_W, 2, 1));
            average_hinge_reg_X_gradient_norm = arma::mean(arma::vecnorm(reg_X_term, 2, 1));
            average_hinge_reg_Omega_gradient_norm = arma::mean(arma::vecnorm(reg_Omega_term, 2, 1));
        }
        tmp_X = (new_X - current_learning_rate * der_X); // estimate new X given derivative

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

        // Notice that we don't use inverse here. Moreover, the first row of Omega is automatically positive
        // if X's first column is positive.
        // TODO (cjlee): calculate tmp_Omega is found by DSA.
        // tmp_Omega = arma::pinv(tmp_X);

        // Ensure if first row of Omega is all positive
        // if (arma::any( tmp_Omega.row(0) <= 0)) {
        //     for (int c=0; c < cell_types; c++) {
        //         double matrix_value =  tmp_Omega(0,c);
        //         if (matrix_value <= 0) {
        //             int shrink_iteration = 0;
        //             while((matrix_value <= 0)& (shrink_iteration < shrink_limit)) {
        //             der_X /=  2;
        //             der_X.row(c) *= 2;
        //             tmp_X = (new_X - current_learning_rate * der_X);
        //             tmp_Omega = arma::pinv(tmp_X);
        //             matrix_value =  tmp_Omega(0,c);
        //             shrink_iteration++;
        //            }
        //         if (shrink_iteration != shrink_limit) {
        //             // if we were able to find the solution. accept these new X and Omega
        //             // Do nothing its ok
        //             } else {
        //                 spdl::warn("Iteration {} Couldn't find good inverse X for the row {}, reject X update for this row", itr_, c);
        //                 arma::rowvec only_good_row  = der_X.row(c);
        //                 der_X.fill(0);
        //                 der_X.row(c) = only_good_row;
        //                 tmp_X = (new_X - current_learning_rate * der_X);
        //             }
        //         }
        //     }
        //     new_Omega = tmp_Omega;
        //     new_X = tmp_X;
        // } else {
        //     new_Omega = tmp_Omega;
        //     new_X = tmp_X;
        // }
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
                                               
        errors_statistics.row(itr_) = arma::rowvec{
            current_errors["deconv_error"],            //1
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



    return Rcpp::List::create(
        Rcpp::Named("new_X") = final_X,
        Rcpp::Named("new_Omega") = final_Omega,
        Rcpp::Named("new_D_w") = new_D_w,
        Rcpp::Named("new_D_h") = new_D_h,
        Rcpp::Named("errors_statistics") = errors_statistics,
        Rcpp::Named("points_statistics_X") = points_statistics_X,
        Rcpp::Named("points_statistics_Omega") = points_statistics_Omega,
        Rcpp::Named("points_statistics_Dw") = points_statistics_Dw
    );
}