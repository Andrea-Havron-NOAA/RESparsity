#include "core/laplace/model_analysis_report.hpp"
#include "core/laplace/laplace_objective.hpp"
#include "core/model/parameter_partition.hpp"
#include "core/model/quadra_model.hpp"
#include "include/quadra/stats.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

DECLARE_ADGRAPH();

class ScalingAr1Model : public quadra::QuadraModel<ScalingAr1Model> {
public:
  explicit ScalingAr1Model(int n) : y_m(static_cast<std::size_t>(n)) {
    for (int i = 0; i < n; ++i) {
      const double time = static_cast<double>(i + 1);
      y_m[static_cast<std::size_t>(i)] =
          0.5 * std::sin(0.13 * time) + 0.25 * std::cos(0.07 * time);
    }
    parameters_m.add("unconstrained_phi", 0.4,
                     quadra::ParameterTransform::Identity, false);
    for (int i = 0; i < n; ++i) {
      parameters_m.add("x_" + std::to_string(i), 0.0,
                       quadra::ParameterTransform::Identity, true);
    }
  }

  std::vector<std::string> parameter_names_impl() const {
    return parameters_m.names();
  }
  const quadra::ParameterSet &parameters() const { return parameters_m; }
  int size() const { return static_cast<int>(y_m.size()); }

  template <typename T>
  T evaluate_impl(const std::vector<T> &parameters,
                  quadra::ModelReportContext &) const {
    const T phi =
        quadra::stats::correlation_from_unconstrained(parameters[0]);
    std::vector<T> states(parameters.begin() + 1, parameters.end());
    T nll = quadra::stats::ar1_stationary_nll(states, T(0.0), phi, T(0.8));
    for (std::size_t i = 0; i < y_m.size(); ++i) {
      nll += quadra::stats::normal_nll(T(y_m[i]), states[i], T(0.5));
    }
    return nll;
  }

private:
  std::vector<double> y_m;
  quadra::ParameterSet parameters_m;
};

template <class Function> double elapsed_ms(Function &&function) {
  const auto start = std::chrono::steady_clock::now();
  function();
  const auto end = std::chrono::steady_clock::now();
  return std::chrono::duration<double, std::milli>(end - start).count();
}

double mean(const std::vector<double> &values) {
  double total = 0.0;
  for (double value : values) {
    total += value;
  }
  return total / static_cast<double>(values.size());
}

int main(int argc, char **argv) {
  std::cout << std::unitbuf;
  std::cout << std::setprecision(17);
  std::cout << "engine,n,setup_ms,first_laplace_ms,warm_laplace_mean_ms,"
               "gradient_setup_ms,gradient_ms,gradient,gradient_method,"
               "hessian_nnz,hessian_density,"
               "hessian_unique_nnz,hessian_bandwidth,objective,"
               "objective_phase_ms,factorization_phase_ms,"
               "sensitivity_phase_ms,hdot_phase_ms,"
               "hdot_validation_ms,hdot_direction_setup_ms,"
               "hdot_reverse_ms,hdot_contraction_ms,"
               "trace_phase_ms,gradient_internal_total_ms,"
               "optimization_ms,estimate,optimization_converged,repetitions,"
               "gradient_repetitions,success,detected_structure,"
               "selected_backend,symbolic_reuse,active_directions,"
               "hdot_workers,objective_tape_rebuilds,hdot_tape_rebuilds\n";

  std::vector<int> sizes = {30, 100, 300, 1000};
  if (argc > 1)
    sizes = {std::stoi(argv[1])};

  for (int n : sizes) {
    ScalingAr1Model *model_pointer = nullptr;
    const double setup_ms = elapsed_ms([&]() {
      model_pointer = new ScalingAr1Model(n);
    });
    ScalingAr1Model &model = *model_pointer;

    quadra::LaplaceObjectiveOptions objective_options;
    objective_options.include_constant_m = true;
    objective_options.newton_m.gradient_tolerance_m = 1e-9;
    objective_options.newton_m.step_tolerance_m = 1e-11;

    quadra::LaplaceObjectiveResult result;
    std::vector<double> random_start(static_cast<std::size_t>(n), 0.0);
    quadra::stats::LaplaceEvaluator<ScalingAr1Model> automatic_evaluator(
        model, random_start, model.parameters(), objective_options);
    const double first_ms = elapsed_ms([&]() {
      result = automatic_evaluator.evaluate({0.4});
    });

    const int repetitions = n <= 30 ? 1000 : (n <= 100 ? 300 :
                            (n <= 300 ? 100 : 30));
    std::vector<double> warm_times;
    warm_times.reserve(static_cast<std::size_t>(repetitions));
    for (int repetition = 0; repetition < repetitions; ++repetition) {
      random_start = result.u_hat_m;
      warm_times.push_back(elapsed_ms([&]() {
        result = automatic_evaluator.evaluate({0.4});
      }));
    }

    const int gradient_repetitions =
        n <= 30 ? 20 : (n <= 100 ? 10 : (n <= 300 ? 5 : 3));
    std::unique_ptr<quadra::stats::ExactLaplaceEvaluator<ScalingAr1Model>>
        gradient_evaluator;
    const double gradient_setup_ms = elapsed_ms([&]() {
      gradient_evaluator.reset(
          new quadra::stats::ExactLaplaceEvaluator<ScalingAr1Model>(
              model, {0.4}, result.u_hat_m, model.parameters(),
              objective_options));
    });
    quadra::stats::ExactLaplaceResult gradient;
    const double gradient_total_ms = elapsed_ms([&]() {
      for (int repetition = 0; repetition < gradient_repetitions;
           ++repetition) {
        gradient = gradient_evaluator->evaluate({0.4}, result.u_hat_m);
      }
    });
    const double gradient_ms =
        gradient_total_ms / static_cast<double>(gradient_repetitions);

    quadra::stats::LaplaceOptimizerOptions optimizer_options;
    optimizer_options.max_iterations = 30;
    optimizer_options.gradient_tolerance = n <= 100 ? 1e-6 : 1e-4;
    optimizer_options.step_tolerance = 1e-10;
    optimizer_options.initial_step_scale = 0.5;
    quadra::stats::LaplaceOptimizerResult optimization;
    double optimization_ms = std::nan("");
    double estimate = std::nan("");
    int optimization_converged = -1;
    optimization_ms = elapsed_ms([&]() {
      optimization = quadra::stats::optimize_laplace(
          *gradient_evaluator, {0.4}, optimizer_options);
    });
    estimate = optimization.fixed[0];
    optimization_converged = optimization.converged ? 1 : 0;
    if (!optimization.converged) {
      std::cerr << "n=" << n << " optimization: " << optimization.message
                << " gradient_norm=" << optimization.gradient_norm << "\n";
    }

    const double density =
        static_cast<double>(result.hessian_random_m.nonZeros()) /
        static_cast<double>(n * n);
    int unique_nonzeros = 0;
    int bandwidth = 0;
    for (int outer = 0; outer < result.hessian_random_m.outerSize(); ++outer) {
      for (Eigen::SparseMatrix<double>::InnerIterator entry(
               result.hessian_random_m, outer);
           entry; ++entry) {
        if (entry.row() <= entry.col()) {
          ++unique_nonzeros;
        }
        bandwidth = std::max(
            bandwidth,
            static_cast<int>(std::abs(entry.row() - entry.col())));
      }
    }
    const bool success = result.converged_m && result.logdet_ok_m &&
                         gradient.success && optimization.converged;
    const auto diagnostics =
        quadra::laplace::analyze_hessian_structure(result.hessian_random_m);

    std::cout << "quadra," << n << "," << setup_ms << "," << first_ms << ","
              << mean(warm_times) << "," << gradient_setup_ms << ","
              << gradient_ms << "," << gradient.gradient[0]
              << ",exact_marginal_autodiff,"
              << result.hessian_random_m.nonZeros() << "," << density << ","
              << unique_nonzeros << "," << bandwidth << ","
              << result.laplace_objective_m << ","
              << gradient.timings.objective_ms << ","
              << gradient.timings.factorization_ms << ","
              << gradient.timings.mode_sensitivity_ms << ","
              << gradient.timings.hdot_ms << ","
              << gradient.timings.hdot_validation_ms << ","
              << gradient.timings.hdot_direction_setup_ms << ","
              << gradient.timings.hdot_reverse_ms << ","
              << gradient.timings.hdot_contraction_ms << ","
              << gradient.timings.trace_ms << ","
              << gradient.timings.total_ms << "," << optimization_ms << ","
              << estimate << "," << optimization_converged << ","
              << repetitions << ","
              << gradient_repetitions << ","
              << (success ? 1 : 0) << ","
              << quadra::laplace::ToString(diagnostics.structure) << ","
              << quadra::laplace::ToString(diagnostics.backend) << ","
              << (diagnostics.supports_symbolic_reuse ? 1 : 0) << ","
              << gradient.active_directions.size() << ","
              << gradient_evaluator->hdot_worker_count() << ","
              << gradient_evaluator->objective_tape_rebuild_count() << ","
              << gradient_evaluator->hdot_tape_rebuild_count() << "\n";

    delete model_pointer;
  }
  return 0;
}
