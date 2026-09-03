#pragma once

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(spdl)]]
#include <RcppArmadillo.h>
#include <spdl.h>

#include <string>
#include <vector>
#include <memory>
#include <unordered_map>

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

struct EnsembleSpec {
    std::string name;
    std::vector<std::string> components;
};

class IErrorCalculator {
public:
    virtual ~IErrorCalculator() = default;
    virtual std::vector<Metric> calculate(const OptimizationState& state) const = 0;
};

class CompositeErrorCalculator : public IErrorCalculator {
private:
    std::vector<std::shared_ptr<IErrorCalculator>> calculators_;
    std::vector<EnsembleSpec> ensemble_specs_;

public:
    void add_calculator(std::shared_ptr<IErrorCalculator> calc) {
        calculators_.push_back(calc);
    }

    void add_ensemble_metric(const std::string& ensemble_name, const std::vector<std::string>& component_names) {
        ensemble_specs_.push_back({ensemble_name, component_names});
    }

    std::vector<Metric> calculate(const OptimizationState& state) const override {
        std::vector<Metric> combined;
        std::unordered_map<std::string, double> lookup;

        for (const auto& calc : calculators_) {
            auto metrics = calc->calculate(state);
            for (const auto& m : metrics) {
                lookup[m.name] = m.value;
                combined.push_back(m);
            }
        }

        for (const auto& spec : ensemble_specs_) {
            double ensemble_val = 0.0;
            for (const auto& comp_name : spec.components) {
                auto it = lookup.find(comp_name);
                if (it != lookup.end()) {
                    ensemble_val += it->second;
                } else {
                    spdl::warn("Ensemble metric '{}' requested component '{}' which was not found in calculated metrics.", spec.name, comp_name);
                }
            }
            combined.push_back({spec.name, ensemble_val});
        }

        return combined;
    }
};

inline double get_metric_value(const std::vector<Metric>& metrics, const std::string& name, double default_val = 0.0) {
    for (const auto& m : metrics) {
        if (m.name == name) {
            return m.value;
        }
    }
    return default_val;
}

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
