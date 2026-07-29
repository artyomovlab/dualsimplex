#pragma once

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <string>
#include <vector>
#include <memory>

struct OptimizationState {
    const arma::mat& X;
    const arma::mat& Omega;
    const arma::mat& D_w;
    const arma::mat& D_h;
    const arma::mat& SVRt;
    const arma::mat& R;
    const arma::mat& S;
    double coef_hinge_H{0.0};
    double coef_hinge_W{0.0};
};

struct Metric {
    std::string name;
    double value;
};

class IErrorCalculator {
public:
    virtual ~IErrorCalculator() = default;
    virtual std::vector<Metric> calculate(const OptimizationState& state) const = 0;
};

class CompositeErrorCalculator : public IErrorCalculator {
private:
    std::vector<std::shared_ptr<IErrorCalculator>> calculators_;
public:
    void add_calculator(std::shared_ptr<IErrorCalculator> calc) {
        calculators_.push_back(calc);
    }

    std::vector<Metric> calculate(const OptimizationState& state) const override {
        std::vector<Metric> combined;
        for (const auto& calc : calculators_) {
            auto metrics = calc->calculate(state);
            combined.insert(combined.end(), metrics.begin(), metrics.end());
        }
        return combined;
    }
};

class DeconvErrorCalculator : public IErrorCalculator {
public:
    std::vector<Metric> calculate(const OptimizationState& state) const override;
};

class HingeProportionsErrorCalculator : public IErrorCalculator {
public:
    std::vector<Metric> calculate(const OptimizationState& state) const override;
};

class HingeBasisErrorCalculator : public IErrorCalculator {
public:
    std::vector<Metric> calculate(const OptimizationState& state) const override;
};

class ScaleNormErrorCalculator : public IErrorCalculator {
public:
    std::vector<Metric> calculate(const OptimizationState& state) const override;
};
