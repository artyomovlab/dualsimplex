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



Rcpp::List optimize_alignment_pgd(
    const arma::mat& initial_X,
    const arma::mat& initial_Omega,
    const arma::mat& initial_D_w,
    const arma::mat& SVRt,
    const arma::mat& R,
    const arma::mat& S,
    const double coef_der_X,
    double coef_hinge_W,
    double coef_hinge_H,
    const int k,
    const double N,
    const double M,
    const int max_iteration,
    const double convergence_tol,
    const int patience,
    const double decade_rate,
    const int max_drop,
    const double epsilon,
    const bool debug_stats
) {
    // Setup local parameters, return variables, metrics, logger and intermediate variables
    arma::mat points_statistics_X(max_iteration + 1, k * k, arma::fill::zeros);
    arma::mat points_statistics_Omega(max_iteration + 1, k * k, arma::fill::zeros);
    arma::mat points_statistics_Dw(max_iteration + 1, k, arma::fill::zeros);

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
    metadata["algorithm"] = "optimize_alignment_pgd";
    metadata["max_iteration"] = std::to_string(max_iteration);
    metadata["coef_hinge_H"] = std::to_string(coef_hinge_H);
    metadata["coef_hinge_W"] = std::to_string(coef_hinge_W);
    rcpp_logger.log_metadata(metadata);

    double current_learning_rate = coef_der_X;
    double best_error_value = arma::datum::inf;
    int best_error_iteration = 0;
    double current_error_value;
    double prev_error_value;
    int consecutive_small_changes = 0;
    int num_drops = 0;

    double c = N / M;
    arma::mat X = initial_X;
    arma::mat Omega = initial_Omega;
    arma::vec d_w = arma::vectorise(initial_D_w);  // let's store diagonal matrix as a vector for convenience
    arma::vec d_h = d_w * c;
    arma::mat sigma_ss = SVRt;
    arma::mat sigma_fs = c * sigma_ss;
    

    arma::mat grad_Q(k, k, arma::fill::zeros);
    arma::mat grad_hinge_W(k, k, arma::fill::zeros);
    arma::mat grad_hinge_H(k, k, arma::fill::zeros);

    // CJLee(TODO): calculate our starting Q
    arma::mat Q, Q_illegal, Q_cand, U_svd, V_svd;
    arma::vec s_svd;
    Q = arma::diagmat(arma::sqrt(c * d_w)) * X;

    // Log initial state (iteration 0)
    OptimizationState init_state{
        X, Omega, d_w, d_h, SVRt, R, S, coef_hinge_H, coef_hinge_W
    };
    std::vector<Metric> init_metrics = composite_calc.calculate(init_state);

    current_error_value = get_metric_value(init_metrics, "total_error");
    best_error_value = current_error_value;
    best_error_iteration = 0;
    prev_error_value = current_error_value;

    init_metrics.push_back({"learning_rate", current_learning_rate});
    init_metrics.push_back({"sum_", accu(d_w) / M});
    init_metrics.push_back({"best_error_value", best_error_value});
    init_metrics.push_back({"best_error_iteration", static_cast<double>(best_error_iteration)});
    rcpp_logger.log_metrics(0, init_metrics);

    points_statistics_X.row(0) = X.as_row();
    points_statistics_Omega.row(0) = Omega.as_row();
    points_statistics_Dw.row(0) = d_w.as_row();


    // -------------------------------------------------------------
    // Optimization
    // -------------------------------------------------------------

    // Start initial inverse search
    int itr_ = 1;
    while (itr_ <= max_iteration) {
        // We use projected gradient descent for this optimization problem.
        // The other option is to use a smoothed hinge loss (c.f. https://mathoverflow.net/questions/51370/smooth-approximation-of-the-hinge-loss-function) paired
        // with Riemannian gradient descent.

        // Following are the optimization logic using projected gradient descent:
        // Staring from a point, Q(t):
        // 1. Calculate (sub)graident from hinge loss.
        grad_Q.zeros();
        if (0 < coef_hinge_W) {
            grad_Q += coef_hinge_W * subgradient(S.t() * sigma_ss * Q.t()).t() * S.t() * sigma_ss;
        }
        if (0 < coef_hinge_H) {
            grad_Q +=  coef_hinge_H * subgradient(Q*R) * R.t();
        }

        // 2. Update to Q(t+1)
        double current_lr = current_learning_rate;
        bool step_accepted = false;
        while (current_lr > 1e-16) {  // since we half the learning rate at each failure, it will terminate at arount 55th failure.
            // 2a. Unconstrained Step (Euclidean Space)
            Q_illegal = Q - current_lr * grad_Q;

            // 2b. Retract to the closest point on manifold via SVD.
            // Q_illegal = U_s * Sigma * V_s^T  ->  Q_cand = U_s * V_s^T
            arma::svd(U_svd, s_svd, V_svd, Q_illegal);
            Q_cand = U_svd * V_svd.t();
            
            // 2c. Hard constraint check (positivity of first column)
            if (arma::all(Q_cand.col(0) > 0)) {
                Q = Q_cand;
                step_accepted = true;
                break; // Exit backtracking loop
            } else {
                // If element in the first column is negative, shrink step
                current_lr *= 0.5;
            }
        }

        if (!step_accepted) {
            // If step size vanishes, we are stuck against the wall or converged
            spdl::warn("Any gradient step gives bad Q.");
            break;
        }

        // 3. Update X(t+1), Omega(t+1) and D(t+1) from Q(t+1).
        d_w = M * arma::square(Q.col(0));
        d_h = d_w * c;
        X = arma::diagmat(1 / arma::sqrt(c * d_w)) * Q;
        Omega = sigma_fs * X.t();


        // -------------------------------------------------------------
        // Logging
        // -------------------------------------------------------------
        OptimizationState current_state{
            X, Omega, d_w, d_h, SVRt, R, S, coef_hinge_H, coef_hinge_W
        };
        std::vector<Metric> metrics = composite_calc.calculate(current_state);
        metrics.push_back({"learning_rate", current_learning_rate});
        metrics.push_back({"sum_", accu(d_w) / M});

        current_error_value = get_metric_value(metrics, "total_error");
        // TODO: Do we still need to track best_error...?
        if (current_error_value < best_error_value) {
            best_error_iteration = itr_;
            best_error_value = current_error_value;
        }
        metrics.push_back({"best_error_value", best_error_value});
        metrics.push_back({"best_error_iteration", static_cast<double>(best_error_iteration)});

        rcpp_logger.log_metrics(itr_, metrics);

        // Log trajectory
        points_statistics_X.row(itr_) = X.as_row();
        points_statistics_Omega.row(itr_) = Omega.as_row();
        points_statistics_Dw.row(itr_) = d_w.as_row();


        // -------------------------------------------------------------
        // Optimization Control Flow & Learning Rate Scheduling
        // -------------------------------------------------------------
        double relative_change = (prev_error_value - current_error_value) / (std::abs(prev_error_value) + epsilon);
        if (relative_change < convergence_tol) {
            consecutive_small_changes++;
        } else {
            consecutive_small_changes = 0;
        }

        if (consecutive_small_changes >= patience) {
            current_learning_rate *= decade_rate;
            num_drops++;
            consecutive_small_changes = 0;
            spdl::info("Plateau detected at iteration {}. Reducing learning rate to {}", itr_, current_learning_rate);
        }

        prev_error_value = current_error_value;


        if (num_drops >= max_drop) {
            spdl::info("Learning rate dropped {} times (max_drop reached). Stopping optimization at iteration {}.", num_drops, itr_);
            itr_++;
            break;
        }

        itr_++;
    }
    
    if (itr_ <= max_iteration) {
        points_statistics_X.resize(itr_, points_statistics_X.n_cols);
        points_statistics_Omega.resize(itr_, points_statistics_Omega.n_cols);
        points_statistics_Dw.resize(itr_, points_statistics_Dw.n_cols);
    } else {
        spdl::info("Reached maximum number of iterations: {}. Stopping optimization.", itr_ - 1);
    }


    return Rcpp::List::create(
        // kept the names here for compatibility. They will be changed after migration
        Rcpp::Named("new_X") = X,
        Rcpp::Named("new_Omega") = Omega,
        Rcpp::Named("new_D_w") = d_w,
        Rcpp::Named("new_D_h") = d_h,
        Rcpp::Named("points_statistics_X") = points_statistics_X,
        Rcpp::Named("points_statistics_Omega") = points_statistics_Omega,
        Rcpp::Named("points_statistics_Dw") = points_statistics_Dw,
        Rcpp::Named("history") = rcpp_logger.to_dataframe(),
        Rcpp::Named("metadata") = rcpp_logger.to_metadata_list()
    );
}

