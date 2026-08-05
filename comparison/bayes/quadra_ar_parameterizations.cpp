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

enum class Parameterization { Centered, CenteredManual, Scaled, Standard };

class BayesianArParityModel
    : public quadra::QuadraModel<BayesianArParityModel> {
public:
  explicit BayesianArParityModel(const Parameterization parameterization)
      : parameterization_(parameterization) {
    parameters_.add("unconstrained_phi", 0.4,
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
    const T phi =
        quadra::stats::correlation_from_unconstrained(parameters[0]);
    const T innovation_sd = T(0.8) * sqrt(T(1.0) - phi * phi);
    std::vector<T> latent(parameters.begin() + 1, parameters.end());
    std::vector<T> states(4);
    T nll = T(0.0);

    if (parameterization_ == Parameterization::Centered ||
        parameterization_ == Parameterization::CenteredManual) {
      states = latent;
      T previous = T(-0.2);
      for (int i = 0; i < 4; ++i) {
        nll += quadra::stats::normal_nll(states[i], phi * previous,
                                         innovation_sd);
        previous = states[i];
      }
    } else {
      T previous = T(-0.2);
      for (int i = 0; i < 4; ++i) {
        if (parameterization_ == Parameterization::Standard) {
          nll += quadra::stats::normal_nll(latent[i], T(0.0), T(1.0));
          states[i] = phi * previous + innovation_sd * latent[i];
        } else {
          nll += quadra::stats::normal_nll(latent[i], T(0.0), innovation_sd);
          states[i] = phi * previous + latent[i];
        }
        previous = states[i];
      }
    }

    const double observations[4] = {-0.1, 0.6, 1.0, 0.5};
    for (int i = 0; i < 4; ++i)
      nll += quadra::stats::normal_nll(T(observations[i]), states[i], T(0.5));
    return nll;
  }

  std::vector<double> states_from_mode(const double unconstrained_phi,
                                       const std::vector<double> &mode) const {
    if (parameterization_ == Parameterization::Centered ||
        parameterization_ == Parameterization::CenteredManual)
      return mode;
    const double phi = 2.0 / (1.0 + std::exp(-unconstrained_phi)) - 1.0;
    const double scale = 0.8 * std::sqrt(1.0 - phi * phi);
    std::vector<double> states(4);
    double previous = -0.2;
    for (int i = 0; i < 4; ++i) {
      const double innovation = parameterization_ == Parameterization::Standard
                                    ? scale * mode[i]
                                    : mode[i];
      states[i] = phi * previous + innovation;
      previous = states[i];
    }
    return states;
  }

private:
  Parameterization parameterization_;
  quadra::ParameterSet parameters_;
};

void run_case(const char *label, const Parameterization parameterization,
              const std::string &diagnostics_dir) {
  BayesianArParityModel model(parameterization);
  quadra::LaplaceObjectiveOptions options;
  options.include_constant_m = true;
  options.newton_m.gradient_tolerance_m = 1e-11;
  options.newton_m.step_tolerance_m = 1e-13;

  quadra::stats::ExactLaplaceEvaluator<BayesianArParityModel> evaluator(
      model, {0.4}, std::vector<double>(4, 0.0), model.parameters(), options);
  const auto result = evaluator.evaluate({0.4}, std::vector<double>(4, 0.0));
  if (!result.success)
    throw std::runtime_error(std::string("Laplace failure for ") + label);

  const auto diagnostics = quadra::laplace::analyze_hessian_structure(
      result.objective.hessian_random_m);
  const auto states = model.states_from_mode(0.4, result.objective.u_hat_m);
  if (!diagnostics_dir.empty()) {
    std::vector<std::string> names;
    for (int i = 0; i < 4; ++i)
      names.push_back("latent_" + std::to_string(i));
    comparison::write_quadra_diagnostics(
        diagnostics_dir + "/bayes_ar_" + label,
        std::string("Bayesian AR Diagnostics: ") + label,
        "Quadra diagnostics at the parameterization parity evaluation point. "
        "The latent mode converged; the fixed effects were held fixed and were "
        "not optimized.",
        result.objective.laplace_objective_m, result.gradient,
        result.objective.newton_iterations_m, result.objective.converged_m,
        result.objective.message_m, result.objective.hessian_random_m,
        result.objective.u_hat_m, names, 1, true);
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
  run_case("centered", Parameterization::Centered, diagnostics_dir);
  run_case("centered_manual", Parameterization::CenteredManual,
           diagnostics_dir);
  run_case("scaled_innovations", Parameterization::Scaled, diagnostics_dir);
  run_case("standard_innovations", Parameterization::Standard,
           diagnostics_dir);
  return 0;
}
