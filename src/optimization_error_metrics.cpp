#include "optimization_error_metrics.h"
#include <cmath>

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
