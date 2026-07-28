# Architectural Blueprint: C++ Error Tracking & R Integration Refactor

This document outlines the architectural plan to replace the rigid `calcErrors` logic and hardcoded `errors_statistics` matrices with a modular, extensible, and R-friendly tracking system.

## Approved Decisions

> [!NOTE]
> Based on the review, the following architectural choices have been finalized:
> - **Serialization**: We will not introduce Apache Arrow. We will start with memory-mapped `Rcpp::DataFrame` serialization for typical runs, and fallback to on-disk `CSV` streaming for massive runs.
> - **Data Shape**: The log structure will natively adopt the "tidy" long-format (`iteration`, `metric_name`, `value`) to flexibly support any arbitrary errors without relying on sparse matrices.

## Proposed Architecture

### Milestone 1: Core C++ Abstraction & Interface Design

To eliminate the `if/else` explosion and rigid 21-column matrices, we will implement the **Strategy** and **Composite** design patterns. This decouples the mathematical calculations from the main optimization loop.

1. **State Encapsulation**
   Instead of passing individual `arma::mat` references (X, Omega, D_w) to various functions, we wrap the current iteration's state in an immutable `OptimizationState` struct. This makes extending calculators trivial since the signature never changes.

2. **The `IErrorCalculator` Interface (Strategy Pattern)**
   We will define a standard interface for any error component.
   ```cpp
   struct Metric {
       std::string name;
       double value;
   };

   class IErrorCalculator {
   public:
       virtual ~IErrorCalculator() = default;
       virtual std::vector<Metric> calculate(const OptimizationState& state) const = 0;
   };
   ```

3. **Concrete Strategies & The Composite Pattern**
   Specific error terms (e.g., `PositivityErrorW`, `AlignmentError`, `HingeError`) will be implemented as isolated classes inheriting from `IErrorCalculator`. 
   A `CompositeErrorCalculator` will hold a collection of these strategies. When a new algorithm is initialized, it dynamically registers the specific error strategies it needs into the composite calculator.

### Milestone 2: Logging & Serialization Strategy

To handle the dynamic nature of these logs, we will abstract the storage mechanism away from `arma::mat`.

1. **The `ILogger` Interface (Observer Pattern)**
   The optimization loop will report metrics to an `ILogger` interface, completely unaware of how the data is stored.
   ```cpp
   class ILogger {
   public:
       virtual void log_metrics(int iteration, const std::vector<Metric>& metrics) = 0;
       virtual void log_metadata(const std::map<std::string, std::string>& metadata) = 0;
       virtual ~ILogger() = default;
   };
   ```

2. **Long-Format Data Model (Tidy Data)**
   Logs will be stored in a flat, long format: `[iteration, metric_name, value]`. This means `optimize_alignment` can log its unique `coef_alignment` error without forcing other algorithms to have a blank column for it.

3. **Concrete Loggers**
   - `CsvLogger`: Streams long-format data to a CSV file incrementally (O(1) memory overhead). Highly robust and R-friendly, this will be the primary choice for extremely long optimization runs.
   - `RcppLogger`: Stores data in standard C++ vectors to build an `Rcpp::DataFrame` in memory for immediate use within R.
### Milestone 3: R Interoperability Layer

We need to return both the optimization results and the rich logs seamlessly back to R using `Rcpp`.

1. **In-Memory Tracking (`RcppLogger`)**
   If the number of iterations fits comfortably in RAM, we implement an `RcppLogger` that buffers the `[iteration, metric_name, value]` triplets in standard C++ vectors. At the end of the run, it binds them into an `Rcpp::DataFrame`.

2. **Rich Metadata Dictionary**
   We will capture hyperparameters (`coef_hinge_W`, `coef_alignment`), algorithm name, and start/end timestamps in a `std::map<std::string, std::string>`. This maps perfectly to a named `Rcpp::List` (an R `list()`).

3. **Unified Return Structure**
   The `Rcpp::export` functions will return a standardized `Rcpp::List`:
   ```R
   List(
     _["parameters"] = List(_["X"] = final_X, _["Omega"] = final_Omega, ...),
     _["history"]    = RcppDataFrame,  // The long-format metrics
     _["metadata"]   = RcppList        // Hyperparameters and run details
   )
   ```

### Milestone 4: Back Compatibility Layer [NEW]

- Build a dedicated C++ or R backward compatibility module.
- Design it to intercept and seamlessly convert the newly refactored data/output formats back into the legacy matrix structure.
- The outcome must ensure that all surrounding legacy code blocks and external dependencies continue to operate with zero downstream friction.

### Milestone 5: Verification, Migration, and Test Strategy

To prevent regressions, the migration will be iterative and safely scoped.

1. **Phase 1: Shadow Logging (Dual Run)**
   Implement `OptimizationState`, `IErrorCalculator`, and `RcppLogger`. Run the new calculators side-by-side with the legacy `arma::mat errors_statistics` inside `optimize_alignment`.
2. **Phase 2: R-Side Validation**
   Compare the final legacy matrix to the new long-format DataFrame in R (using `tidyr::pivot_wider` to match shapes). Ensure absolute numerical equivalence.
3. **Phase 3: Cleanup & Scaling**
   Remove the hardcoded 21-column matrices and old `calcErrors` functions. Extend the pattern to other algorithms in the repository.
4. **Phase 4: C++ Unit Testing**
   Leverage the newly decoupled architecture to write isolated C++ unit tests for individual error math without needing to run the full optimization loop.
