#include "optimization_utils.h"
#include "matrix_utils.h"
#include "nnls.h"
#include <tuple>
#include "optimization_positivity.h"
#include "optimization_alignment.h"
#include "optimization_error_metrics.h"
#include "optimization_logger.h"


class AlignmentErrorCalculator : public IErrorCalculator {
public:
    std::vector<Metric> calculate(const OptimizationState& state) const override {
        double ratio = static_cast<double>(state.R.n_cols) / static_cast<double>(state.S.n_cols);
        // make sure sigma_fs is a diagonal matrix
        arma::mat sigma_fs = arma::diagmat(arma::diagvec(state.SVRt)) * ratio;
        double alignment_error = arma::norm(state.Omega - sigma_fs * state.X.t(), "fro");
        return {{"alignment_error", alignment_error}};
    }
};



arma::mat subgradient(const arma::mat& Z, const double margin = 0) {
    arma::mat res(Z.n_rows, Z.n_cols, arma::fill::zeros);
    res.elem(arma::find(Z < margin)).fill(-1);
    return res;
}


Rcpp::List optimize_alignment(
    const arma::mat& initial_X,
    const arma::mat& initial_Omega,
    const arma::mat& initial_D_w,
    const arma::mat& SVRt,
    const arma::mat& R,
    const arma::mat& S,
    const double coef_der_X,
    double coef_hinge_W,
    double coef_hinge_H,
    double coef_alignment,
    const int k,
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
    arma::mat points_statistics_X(iterations + 1, k * k, arma::fill::zeros);
    arma::mat points_statistics_Omega(iterations + 1, k * k, arma::fill::zeros);
    arma::mat points_statistics_Dw(iterations + 1, k, arma::fill::zeros);

    CompositeErrorCalculator composite_calc;
    composite_calc.add_calculator(std::make_shared<DeconvErrorCalculator>());
    composite_calc.add_calculator(std::make_shared<HingeProportionsErrorCalculator>());
    composite_calc.add_calculator(std::make_shared<HingeBasisErrorCalculator>());
    composite_calc.add_calculator(std::make_shared<ScaleNormErrorCalculator>());
    composite_calc.add_calculator(std::make_shared<AlignmentErrorCalculator>());
    composite_calc.add_ensemble_metric("total_error", {"deconv_error", "lambda_error", "beta_error"});
    composite_calc.add_ensemble_metric("scaled_total_error", {"deconv_error", "scaled_lambda_error", "scaled_beta_error"});

    RcppLogger rcpp_logger;
    std::map<std::string, std::string> metadata;
    metadata["algorithm"] = "optimize_alignment";
    metadata["iterations"] = std::to_string(iterations);
    metadata["coef_hinge_H"] = std::to_string(coef_hinge_H);
    metadata["coef_hinge_W"] = std::to_string(coef_hinge_W);
    metadata["coef_alignment"] = std::to_string(coef_alignment);
    rcpp_logger.log_metadata(metadata);


    arma::mat X = initial_X;
    arma::mat Omega = initial_Omega;
    arma::mat D_w = initial_D_w;
    arma::mat D_w_sqrt = arma::sqrt(D_w);
    arma::mat D_h = D_w * (N / M);
    arma::vec Sigma = arma::diagvec(SVRt);  // this should be sigma_ss
    arma::vec sqrt_Sigma = arma::sqrt(Sigma);
    arma::vec sigma_fs = Sigma * N / M;

    // calculate X_dtilda and Omega_dtilda in the transformed space for optimization
    arma::mat X_dtilda =  arma::diagmat(D_w_sqrt) * X * arma::diagmat(1 / sqrt_Sigma);
    arma::mat Omega_dtilda =  arma::diagmat(1 / sqrt_Sigma) *  Omega  * arma::diagmat(D_w_sqrt);

    arma::mat align_diff;  // alignment difference
    arma::mat tmp_X, tmp_Omega;


    // derevitives
    arma::mat der_X(k, k, arma::fill::zeros);
    arma::mat hinge_term_H(k, k, arma::fill::zeros);
    arma::mat hinge_term_W(k, k, arma::fill::zeros);
    arma::mat der_dsa(k, k, arma::fill::zeros);

    double shrink_limit = 500; // shouldn't we expose this to user?
    double current_learning_rate = coef_der_X;
    double best_error_value = arma::datum::inf;  // so that all calculated errors will be smaller than this initilzed value
    int best_error_iteration = 0;
    double current_error_value;

    // Start initial inverse search
    tmp_Omega = arma::pinv(X_dtilda);
    if (arma::any(tmp_Omega.row(0) <= 0)) {
       spdl::warn("Couldn't find good initial inverse of X provided. Will try with Omega");
       tmp_X = arma::pinv(Omega_dtilda);
       if (arma::any(tmp_X.col(0) <= 0)) {
            spdl::warn("Couldn't find good initial inverse of Omega provided");
            Rcpp::stop("!!Start with different initialization or ensure X and Omega are inverse!! (try `random_invertible`)");
        }
        else {
            X_dtilda = tmp_X;
        }
    }
    else {
        Omega_dtilda = tmp_Omega;
    }


    // Log initial state (iteration 0)
    OptimizationState init_state{
        X, Omega, D_w, D_h, SVRt, R, S, coef_hinge_H, coef_hinge_W
    };
    std::vector<Metric> init_metrics = composite_calc.calculate(init_state);

    current_error_value = get_metric_value(init_metrics, "total_error");
    best_error_value = current_error_value;
    best_error_iteration = 0;

    init_metrics.push_back({"learning_rate", current_learning_rate});
    init_metrics.push_back({"sum_", accu(D_w) / M});
    init_metrics.push_back({"best_error_value", best_error_value});
    init_metrics.push_back({"best_error_iteration", static_cast<double>(best_error_iteration)});
    rcpp_logger.log_metrics(0, init_metrics);

    points_statistics_X.row(0) = X.as_row();
    points_statistics_Omega.row(0) = Omega.as_row();
    points_statistics_Dw.row(0) = D_w.as_row();

    // here we assume X and Omega are inverse of each other and positive as needed
    int itr_ = 1;
    while ((itr_ < iterations + 1) && (current_learning_rate > convergence_tol)) {
        // -------------------------------------------------------------
        // Main optimization logic
        // -------------------------------------------------------------
        
        // Calculate gradients in the transformed space
        // For data with only pure sample, we could use dual sipmlex alignment to find verticex of sample
        // simplex instead of negativity term.
        align_diff = Omega_dtilda - arma::diagmat(sigma_fs) * X_dtilda.t();  // remember that Omega_dtilda = X_dtilda^-1
        der_dsa = -2 * X_dtilda.t() * align_diff * X_dtilda.t() - 2 * align_diff.t() * arma::diagmat(sigma_fs);
        der_X = coef_alignment * der_dsa;

        // negativity terms for gene and cell simplex. 
        // To save computation time, we only compute it if the coefficient is non-zero.
        if (0 < coef_hinge_W) {  // vertices of cell simplex
            hinge_term_W = (-Omega_dtilda.t()) * arma::diagmat(sqrt_Sigma) * l1_hinge_der_basis_C__(S.t() * arma::diagmat(sqrt_Sigma) * Omega_dtilda, S) * (Omega_dtilda.t());
            der_X += coef_hinge_W * hinge_term_W;
        }
        if (0 < coef_hinge_H) {  // vertices of gene simplex
            hinge_term_H = l1_hinge_der_proportions_C__(X_dtilda * arma::diagmat(sqrt_Sigma) * R, R) * arma::diagmat(sqrt_Sigma);
            der_X +=  coef_hinge_H * hinge_term_H;
        }


        // Update X_dtilda in the transformed space
        tmp_X = (X_dtilda - current_learning_rate * der_X);
        // Ensure if first column of X is all-positive
        if (arma::any(tmp_X.col(0) <= 0)) {
            for (int c=0; c < k; c++) {
                double matrix_value =  tmp_X(c,0);
                if (matrix_value <= 0) {
                int shrink_iteration = 0;
                while((matrix_value <= 0) & (shrink_iteration < shrink_limit)) {
                    der_X.row(c) /=  2;
                    tmp_X = (X_dtilda - current_learning_rate * der_X);
                    matrix_value =  tmp_X(c,0);
                    shrink_iteration++;
                }
                }
            }
            if  (arma::any( tmp_X.col(0) <= 0)) {
                spdl::warn("Any gradient step gives bad X, probably X was bad before");
            }
        }

        // Update Omega_dtilda in the transformed space
        tmp_Omega = arma::pinv(tmp_X);
        // Ensure if first row of Omega is all positive
        if (arma::any( tmp_Omega.row(0) <= 0)) {
            for (int c=0; c < k; c++) {
                double matrix_value =  tmp_Omega(0,c);
                if (matrix_value <= 0) {
                    int shrink_iteration = 0;
                    while((matrix_value <= 0)& (shrink_iteration < shrink_limit)) {
                    der_X /=  2;
                    der_X.row(c) *= 2;
                    tmp_X = (X_dtilda - current_learning_rate * der_X);
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
                        tmp_X = (X_dtilda - current_learning_rate * der_X);
                    }
                }
            }
            Omega_dtilda = tmp_Omega;
            X_dtilda = tmp_X;
        } else {
            Omega_dtilda = tmp_Omega;
            X_dtilda = tmp_X;
        }
 
        // Update Ds, X and Omega in the Sinkhorned space
        std::tie(X_dtilda, Omega_dtilda, D_w_sqrt) = ensure_D_integrity_c(X_dtilda, Omega_dtilda, sqrt_Sigma, N, M);
        X = arma::diagmat(1/D_w_sqrt) * X_dtilda * arma::diagmat(sqrt_Sigma);
        Omega = arma::diagmat(sqrt_Sigma)* Omega_dtilda * arma::diagmat(1/D_w_sqrt);
        D_w = arma::pow(D_w_sqrt, 2);
        D_h = D_w * (N / M);

        // Log errors
        OptimizationState current_state{
            X, Omega, D_w, D_h, SVRt, R, S, coef_hinge_H, coef_hinge_W
        };
        std::vector<Metric> metrics = composite_calc.calculate(current_state);
        metrics.push_back({"learning_rate", current_learning_rate});
        metrics.push_back({"sum_", accu(D_w) / M});

        // -------------------------------------------------------------
        // Optimization Control Flow & Adaptive Step Size
        // -------------------------------------------------------------
        current_error_value = get_metric_value(metrics, "total_error");

        if (current_error_value < best_error_value) {
            best_error_iteration = itr_;
            best_error_value = current_error_value;
        }
        if (itr_ - best_error_iteration > stop_criteria_window) {
            // looks like best solution was not updated for stop_criteria_window iterations. reducing step size.
            current_learning_rate = current_learning_rate / 2;
            best_error_iteration = itr_; // reset iteration counter
        }

        metrics.push_back({"best_error_value", best_error_value});
        metrics.push_back({"best_error_iteration", static_cast<double>(best_error_iteration)});
        rcpp_logger.log_metrics(itr_, metrics);

        // -------------------------------------------------------------
        // Solution Trajectory Recording
        // -------------------------------------------------------------
        points_statistics_X.row(itr_) = X.as_row();
        points_statistics_Omega.row(itr_) = Omega.as_row();
        points_statistics_Dw.row(itr_) = D_w.as_row();

        itr_++;
    }
    
    spdl::info("Optimization completed with number of iterations perfomed: {}", itr_ - 1);
    if (itr_ < iterations + 1) {
        points_statistics_X.resize(itr_, points_statistics_X.n_cols);
        points_statistics_Omega.resize(itr_, points_statistics_Omega.n_cols);
        points_statistics_Dw.resize(itr_, points_statistics_Dw.n_cols);
    }

    return Rcpp::List::create(
        // kept the names here for compatibility. They will be changed after migration
        Rcpp::Named("new_X") = X,
        Rcpp::Named("new_Omega") = Omega,
        Rcpp::Named("new_D_w") = D_w,
        Rcpp::Named("new_D_h") = D_h,
        Rcpp::Named("points_statistics_X") = points_statistics_X,
        Rcpp::Named("points_statistics_Omega") = points_statistics_Omega,
        Rcpp::Named("points_statistics_Dw") = points_statistics_Dw,
        Rcpp::Named("history") = rcpp_logger.to_dataframe(),
        Rcpp::Named("metadata") = rcpp_logger.to_metadata_list()
    );
}



// Rcpp::List optimize_alignment_pgd(
//     const arma::mat& initial_X,
//     const arma::mat& initial_Omega,
//     const arma::mat& initial_D_w,
//     const arma::mat& SVRt,
//     const arma::mat& R,
//     const arma::mat& S,
//     const double coef_der_X,
//     double coef_hinge_W,
//     double coef_hinge_H,
//     const int k,
//     const double N,
//     const double M,
//     const int iterations,
//     // double total_regularization_weight,
//     // const double reg_X,
//     // const double reg_Omega,
//     const double convergence_tol,
//     const int stop_criteria_window,
//     const bool debug_stats
// ) {
//     // Setup local parameters, return variables, metrics, logger and intermediate variables
//     double shrink_limit = 500; // shouldn't we expose this to user?

//     arma::mat points_statistics_X(iterations + 1, k * k, arma::fill::zeros);
//     arma::mat points_statistics_Omega(iterations + 1, k * k, arma::fill::zeros);
//     arma::mat points_statistics_Dw(iterations + 1, k, arma::fill::zeros);

//     CompositeErrorCalculator composite_calc;
//     composite_calc.add_calculator(std::make_shared<DeconvErrorCalculator>());
//     composite_calc.add_calculator(std::make_shared<HingeProportionsErrorCalculator>());
//     composite_calc.add_calculator(std::make_shared<HingeBasisErrorCalculator>());
//     composite_calc.add_calculator(std::make_shared<ScaleNormErrorCalculator>());
//     composite_calc.add_calculator(std::make_shared<AlignmentErrorCalculator>());
//     composite_calc.add_ensemble_metric("total_error", {"deconv_error", "lambda_error", "beta_error"});
//     composite_calc.add_ensemble_metric("scaled_total_error", {"deconv_error", "scaled_lambda_error", "scaled_beta_error"});

//     RcppLogger rcpp_logger;
//     std::map<std::string, std::string> metadata;
//     metadata["algorithm"] = "optimize_alignment";
//     metadata["iterations"] = std::to_string(iterations);
//     metadata["coef_hinge_H"] = std::to_string(coef_hinge_H);
//     metadata["coef_hinge_W"] = std::to_string(coef_hinge_W);
//     rcpp_logger.log_metadata(metadata);

//     double current_learning_rate = coef_der_X;
//     double best_error_value = arma::datum::inf;
//     int best_error_iteration = 0;
//     double current_error_value;

//     arma::mat X = initial_X;
//     arma::mat Omega = initial_Omega;
//     arma::mat D_w = initial_D_w;
//     // arma::mat D_w_sqrt = arma::sqrt(D_w);
//     arma::mat D_h = D_w * (N / M);
//     arma::mat sigma_ss = SVRt;  // this should be sigma_ss
//     // arma::mat sigma_fs = sigma_ss * N / M;

//     arma::mat grad_Q(k, k, arma::fill::zeros);
//     arma::mat grad_hinge_W(k, k, arma::fill::zeros);
//     arma::mat grad_hinge_H(k, k, arma::fill::zeros);

//     // CJLee(TODO): calculate our starting Q
//     arma::mat Q, Q_illegal, Q_cand, U_svd, V_svd, s_svd;

//     // following are JUST for backward compatibility of debug mode
//     double average_gradient_norm = 0;
//     double average_hinge_H_gradient_norm = 0;
//     double average_hinge_W_gradient_norm = 0;
//     double average_hinge_reg_X_gradient_norm = 0;
//     double average_hinge_reg_Omega_gradient_norm = 0;
              
//     // Log initial state (iteration 0)
//     OptimizationState init_state{
//         X, Omega, D_w, D_h, SVRt, R, S, coef_hinge_H, coef_hinge_W
//     };
//     std::vector<Metric> init_metrics = composite_calc.calculate(init_state);

//     current_error_value = get_metric_value(init_metrics, "total_error");
//     best_error_value = current_error_value;
//     best_error_iteration = 0;

//     init_metrics.push_back({"learning_rate", current_learning_rate});
//     init_metrics.push_back({"sum_", accu(D_w) / M});
//     init_metrics.push_back({"best_error_value", best_error_value});
//     init_metrics.push_back({"best_error_iteration", static_cast<double>(best_error_iteration)});

//     // For backward compatibility, we can drop them in the future.
//     init_metrics.push_back({"average_gradient_norm", average_gradient_norm});
//     init_metrics.push_back({"average_hinge_H_gradient_norm", average_hinge_H_gradient_norm});
//     init_metrics.push_back({"average_hinge_W_gradient_norm", average_hinge_W_gradient_norm});
//     init_metrics.push_back({"average_hinge_reg_X_gradient_norm", average_hinge_reg_X_gradient_norm});
//     init_metrics.push_back({"average_hinge_reg_Omega_gradient_norm", average_hinge_reg_Omega_gradient_norm});

//     rcpp_logger.log_metrics(0, init_metrics);
//     points_statistics_X.row(0) = X.as_row();
//     points_statistics_Omega.row(0) = Omega.as_row();
//     points_statistics_Dw.row(0) = D_w.as_row();

//     // -------------------------------------------------------------
//     // Optimization
//     // -------------------------------------------------------------

//     // Start initial inverse search
//     int itr_ = 1;
//     while ((itr_ < iterations + 1) && (current_learning_rate > convergence_tol)) {  // I thought `convergence_tol` is for total error, but clearly no here!
//         // We use projected gradient descent for this optimization problem.
//         // The other option is to use a smoothed hinge loss (c.f. https://mathoverflow.net/questions/51370/smooth-approximation-of-the-hinge-loss-function) paired
//         // with Riemannian gradient descent.

//         // Following are the optimization logic using projected gradient descent:
//         // Staring from a point, Q(t):
//         // 1. Calculate (sub)graident from hinge loss.
//         grad_Q.zeros();
//         if (0 < coef_hinge_W) {
//             grad_Q += coef_hinge_W * subgradient(S.t() * sigma_ss * Q.t()).t() * S.t() * sigma_ss;
//         }
//         if (0 < coef_hinge_H) {
//             grad_Q +=  coef_hinge_H * subgradient(Q*R) * R.t();
//         }

//         // 2. Update to Q(t+1)
//         double current_lr = current_learning_rate;
//         bool step_accepted = false;
//         while (current_lr > 1e-12) {
//             // 2a. Unconstrained Step (Euclidean Space)
//             Q_illegal = Q - current_lr * grad_Q;

//             // 2b. Retract to the closest point on manifold via SVD.
//             // Q_illegal = U_s * Sigma * V_s^T  ->  Q_cand = U_s * V_s^T
//             arma::svd(U_svd, s_svd, V_svd, Q_illegal);
//             Q_cand = U_svd * V_svd.t();
            
//             // 2c. Hard constraint check (positivity of first column)
//             if (all(Q_cand.col(0) > 0)) {
//                 Q = Q_cand;
//                 step_accepted = true;
//                 break; // Exit backtracking loop
//             } else {
//                 // If element in the first column is negative, shrink step
//                 current_lr *= 0.5;
//             }
//         }

//         if (!step_accepted) {
//             // If step size vanishes, we are stuck against the wall or converged
//             spdl::warn("Any gradient step gives bad Q.");
//             break;
//         }

//         // 3. Update X(t+1), Omega(t+1) and D(t+1) from Q(t+1).
//         // CJLee(TODO) implement update logic from Q(t+1)
//         // 

//         // Log errors
//         OptimizationState current_state{
//             X, Omega, D_w, D_h, SVRt, R, S, coef_hinge_H, coef_hinge_W
//         };
//         std::vector<Metric> metrics = composite_calc.calculate(current_state);
//         metrics.push_back({"learning_rate", current_learning_rate});
//         metrics.push_back({"sum_", accu(D_w) / M});

//         // -------------------------------------------------------------
//         // Optimization Control Flow & Adaptive Step Size
//         // -------------------------------------------------------------
//         current_error_value = get_metric_value(metrics, "total_error");

//         if (current_error_value < best_error_value) {
//             best_error_iteration = itr_;
//             best_error_value = current_error_value;
//         }
//         if (itr_ - best_error_iteration > stop_criteria_window) {
//             // looks like best solution was not updated for stop_criteria_window iterations. reducing step size.
//             current_learning_rate = current_learning_rate / 2;
//             best_error_iteration = itr_; // reset iteration counter
//         }

//         metrics.push_back({"best_error_value", best_error_value});
//         metrics.push_back({"best_error_iteration", static_cast<double>(best_error_iteration)});

//         // -------------------------------------------------------------
//         // Additional diagnostics, could be further simplefied if we don't need to support backward compatibility
//         // -------------------------------------------------------------
//         if (debug_stats) {
//             average_gradient_norm = arma::mean(arma::vecnorm(der_X, 2, 1));
//             average_hinge_H_gradient_norm = arma::mean(arma::vecnorm(hinge_term_H, 2, 1));
//             average_hinge_W_gradient_norm = arma::mean(arma::vecnorm(hinge_term_W, 2, 1));
//             // average_hinge_reg_X_gradient_norm = arma::mean(arma::vecnorm(reg_X_term, 2, 1));
//             // average_hinge_reg_Omega_gradient_norm = arma::mean(arma::vecnorm(reg_Omega_term, 2, 1));
//         }
//         metrics.push_back({"average_gradient_norm", average_gradient_norm});
//         metrics.push_back({"average_hinge_H_gradient_norm", average_hinge_H_gradient_norm});
//         metrics.push_back({"average_hinge_W_gradient_norm", average_hinge_W_gradient_norm});
//         metrics.push_back({"average_hinge_reg_X_gradient_norm", average_hinge_reg_X_gradient_norm});
//         metrics.push_back({"average_hinge_reg_Omega_gradient_norm", average_hinge_reg_Omega_gradient_norm});


//         rcpp_logger.log_metrics(itr_, metrics);

//         // -------------------------------------------------------------
//         // Solution Trajectory Recording
//         // -------------------------------------------------------------
//         points_statistics_X.row(itr_) = X.as_row();
//         points_statistics_Omega.row(itr_) = Omega.as_row();
//         points_statistics_Dw.row(itr_) = D_w.as_row();

//         itr_++;
//     }
    
//     spdl::info("Optimization completed with number of iterations perfomed: {}", itr_ - 1);
//     if (itr_ < iterations + 1) {
//         points_statistics_X.resize(itr_, points_statistics_X.n_cols);
//         points_statistics_Omega.resize(itr_, points_statistics_Omega.n_cols);
//         points_statistics_Dw.resize(itr_, points_statistics_Dw.n_cols);
//     }

//     arma::mat errors_statistics = rcpp_logger.to_legacy_matrix();  // For backward compatibility, we can drop them in the future.

//     return Rcpp::List::create(
//         // kept the names here for compatibility. They will be changed after migration
//         Rcpp::Named("new_X") = X,
//         Rcpp::Named("new_Omega") = Omega,
//         Rcpp::Named("new_D_w") = D_w,
//         Rcpp::Named("new_D_h") = D_h,
//         Rcpp::Named("errors_statistics") = errors_statistics,

//         Rcpp::Named("points_statistics_X") = points_statistics_X,
//         Rcpp::Named("points_statistics_Omega") = points_statistics_Omega,
//         Rcpp::Named("points_statistics_Dw") = points_statistics_Dw,
//         Rcpp::Named("history") = rcpp_logger.to_dataframe(),
//         Rcpp::Named("metadata") = rcpp_logger.to_metadata_list()
//     );
// }

