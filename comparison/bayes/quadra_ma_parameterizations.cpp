#include "core/laplace/model_analysis_report.hpp"
#include "core/model/quadra_model.hpp"
#include "include/quadra/stats.hpp"
#include "../quadra_diagnostics.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

DECLARE_ADGRAPH();

class BayesianMaParityModel
    : public quadra::QuadraModel<BayesianMaParityModel> {
public:
  explicit BayesianMaParityModel(bool innovations) : innovations_(innovations) {
    parameters_.add("unconstrained_theta", -0.6,
                    quadra::ParameterTransform::Identity, false);
    for (int i = 0; i < 4; ++i)
      parameters_.add("latent_" + std::to_string(i), 0.0,
                      quadra::ParameterTransform::Identity, true);
  }

  std::vector<std::string> parameter_names_impl() const {
    return parameters_.names();
  }
  const quadra::ParameterSet &parameters() const { return parameters_; }

  template <typename T>
  T evaluate_impl(const std::vector<T> &parameters,
                  quadra::ModelReportContext &) const {
    const T theta =
        quadra::stats::correlation_from_unconstrained(parameters[0]);
    std::vector<T> latent(parameters.begin() + 1, parameters.end());
    std::vector<T> innovations(4), states(4);
    T nll = T(0.0);

    if (innovations_) {
      innovations = latent;
      for (int i = 0; i < 4; ++i) {
        nll += quadra::stats::normal_nll(innovations[i], T(0.0), T(0.8));
        states[i] = innovations[i] +
                    (i == 0 ? T(0.0) : theta * innovations[i - 1]);
      }
    } else {
      states = latent;
      innovations[0] = states[0];
      nll += quadra::stats::normal_nll(innovations[0], T(0.0), T(0.8));
      for (int i = 1; i < 4; ++i) {
        innovations[i] = states[i] - theta * innovations[i - 1];
        nll += quadra::stats::normal_nll(innovations[i], T(0.0), T(0.8));
      }
    }

    const double observations[4] = {0.7, 1.0, 0.1, 0.9};
    for (int i = 0; i < 4; ++i)
      nll += quadra::stats::normal_nll(T(observations[i]), states[i], T(0.5));
    return nll;
  }

  std::vector<double> states_from_mode(const double unconstrained_theta,
                                       const std::vector<double> &mode) const {
    if (!innovations_)
      return mode;
    const double theta =
        2.0 / (1.0 + std::exp(-unconstrained_theta)) - 1.0;
    std::vector<double> states(4);
    for (int i = 0; i < 4; ++i)
      states[i] = mode[i] + (i == 0 ? 0.0 : theta * mode[i - 1]);
    return states;
  }

private:
  bool innovations_;
  quadra::ParameterSet parameters_;
};

void run_case(const char *label, const bool innovations,
              const std::string &diagnostics_dir) {
  BayesianMaParityModel model(innovations);
  quadra::LaplaceObjectiveOptions options;
  options.include_constant_m = true;
  options.newton_m.gradient_tolerance_m = 1e-11;
  options.newton_m.step_tolerance_m = 1e-13;
  quadra::stats::ExactLaplaceEvaluator<BayesianMaParityModel> evaluator(
      model, {-0.6}, std::vector<double>(4, 0.0), model.parameters(), options);
  const auto result = evaluator.evaluate({-0.6}, std::vector<double>(4, 0.0));
  if (!result.success)
    throw std::runtime_error(std::string("Laplace failure for ") + label);
  const auto diagnostics = quadra::laplace::analyze_hessian_structure(
      result.objective.hessian_random_m);
  const auto states = model.states_from_mode(-0.6, result.objective.u_hat_m);
  if (!diagnostics_dir.empty()) {
    std::vector<std::string> names;
    for (int i = 0; i < 4; ++i)
      names.push_back("latent_" + std::to_string(i));
    comparison::write_quadra_diagnostics(
        diagnostics_dir + "/bayes_ma_" + label,
        std::string("Bayesian MA Diagnostics: ") + label,
        "Quadra diagnostics at the parameterization parity evaluation point.",
        result.objective.laplace_objective_m, result.gradient,
        result.objective.newton_iterations_m, result.objective.converged_m,
        result.objective.message_m, result.objective.hessian_random_m,
        result.objective.u_hat_m, names, 1);
  }
  std::cout << "quadra," << label << "," << result.objective.laplace_objective_m
            << "," << result.gradient[0] << ","
            << quadra::laplace::ToString(diagnostics.structure) << ","
            << quadra::laplace::ToString(diagnostics.backend);
  for (double state : states)
    std::cout << "," << state;
  std::cout << "\n";
}

int main(int argc, char **argv) {
  const std::string diagnostics_dir = argc > 1 ? argv[1] : "";
  std::cout << std::setprecision(17);
  std::cout << "engine,parameterization,objective,gradient,structure,backend,"
               "state0,state1,state2,state3\n";
  run_case("innovation_good", true, diagnostics_dir);
  run_case("inverse_filter_bad", false, diagnostics_dir);
  return 0;
}
