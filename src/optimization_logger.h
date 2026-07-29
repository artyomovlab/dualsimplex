#pragma once

#include "optimization_error_metrics.h"
#include <vector>
#include <string>
#include <map>
#include <fstream>

class ILogger {
public:
    virtual ~ILogger() = default;
    virtual void log_metrics(int iteration, const std::vector<Metric>& metrics) = 0;
    virtual void log_metadata(const std::map<std::string, std::string>& metadata) = 0;
};

class RcppLogger : public ILogger {
private:
    std::vector<int> iterations_;
    std::vector<std::string> metric_names_;
    std::vector<double> values_;
    std::map<std::string, std::string> metadata_;

public:
    void log_metrics(int iteration, const std::vector<Metric>& metrics) override {
        for (const auto& m : metrics) {
            iterations_.push_back(iteration);
            metric_names_.push_back(m.name);
            values_.push_back(m.value);
        }
    }

    void log_metadata(const std::map<std::string, std::string>& metadata) override {
        metadata_ = metadata;
    }

    Rcpp::DataFrame to_dataframe() const {
        return Rcpp::DataFrame::create(
            Rcpp::Named("iteration") = iterations_,
            Rcpp::Named("metric_name") = metric_names_,
            Rcpp::Named("value") = values_
        );
    }

    Rcpp::List to_metadata_list() const {
        Rcpp::List lst;
        for (const auto& kv : metadata_) {
            lst[kv.first] = kv.second;
        }
        return lst;
    }

    arma::mat to_legacy_matrix() const {
        if (iterations_.empty()) {
            return arma::mat();
        }
        int max_iter = 0;
        for (int it : iterations_) {
            if (it > max_iter) max_iter = it;
        }
        arma::mat mat(max_iter + 1, 21, arma::fill::zeros);

        for (size_t i = 0; i < iterations_.size(); ++i) {
            int r = iterations_[i];
            const std::string& name = metric_names_[i];
            double val = values_[i];

            if (name == "deconv_error") mat(r, 0) = val;
            else if (name == "lambda_error") mat(r, 1) = val;
            else if (name == "beta_error") mat(r, 2) = val;
            else if (name == "D_h_error") mat(r, 3) = val;
            else if (name == "D_w_error") mat(r, 4) = val;
            else if (name == "total_error") mat(r, 5) = val;
            else if (name == "scaled_total_error") mat(r, 6) = val;
            else if (name == "neg_props") mat(r, 7) = val;
            else if (name == "neg_basis") mat(r, 8) = val;
            else if (name == "sum_") mat(r, 9) = val;
            else if (name == "average_norm") mat(r, 10) = val;
            else if (name == "learning_rate") mat(r, 11) = val;
            else if (name == "average_gradient_norm") mat(r, 12) = val;
            else if (name == "average_hinge_H_gradient_norm") mat(r, 13) = val;
            else if (name == "average_hinge_W_gradient_norm") mat(r, 14) = val;
            else if (name == "average_hinge_reg_X_gradient_norm") mat(r, 15) = val;
            else if (name == "average_hinge_reg_Omega_gradient_norm") mat(r, 16) = val;
            else if (name == "best_error_value") mat(r, 17) = val;
            else if (name == "best_error_iteration") mat(r, 18) = val;
            else if (name == "scaled_lambda_error") mat(r, 19) = val;
            else if (name == "scaled_beta_error") mat(r, 20) = val;
        }

        return mat;
    }
};

class CsvLogger : public ILogger {
private:
    std::string filename_;
    std::ofstream out_;

public:
    explicit CsvLogger(const std::string& filename) : filename_(filename) {
        out_.open(filename_);
        if (out_.is_open()) {
            out_ << "iteration,metric_name,value\n";
        }
    }

    ~CsvLogger() override {
        if (out_.is_open()) {
            out_.close();
        }
    }

    void log_metrics(int iteration, const std::vector<Metric>& metrics) override {
        if (!out_.is_open()) return;
        for (const auto& m : metrics) {
            out_ << iteration << "," << m.name << "," << m.value << "\n";
        }
    }

    void log_metadata(const std::map<std::string, std::string>& metadata) override {
        if (!out_.is_open()) return;
        for (const auto& kv : metadata) {
            out_ << "# " << kv.first << ": " << kv.second << "\n";
        }
    }
};
