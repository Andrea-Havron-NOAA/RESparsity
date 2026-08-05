#include "core/model/parameter_partition.hpp"
#include "core/model/quadra_model.hpp"
#include "include/quadra/stats.hpp"
#include "../quadra_diagnostics.hpp"

#include <iomanip>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

DECLARE_ADGRAPH();

class LatentAr1Model : public quadra::QuadraModel<LatentAr1Model> {
public:
  explicit LatentAr1Model(std::vector<double> y) : y_m(std::move(y)) {
    parameters_m.add("unconstrained_phi", 0.4,
                     quadra::ParameterTransform::Identity, false);
    for (std::size_t i = 0; i < y_m.size(); ++i) {
      parameters_m.add("x_" + std::to_string(i), 0.0,
                       quadra::ParameterTransform::Identity, true);
    }
  }

  std::vector<std::string> parameter_names_impl() const {
    return parameters_m.names();
  }

  const quadra::ParameterSet &parameters() const { return parameters_m; }

  template <typename T>
  T evaluate_impl(const std::vector<T> &parameters,
                  quadra::ModelReportContext &context) const {
    const T phi =
        quadra::stats::correlation_from_unconstrained(parameters[0]);
    std::vector<T> states(parameters.begin() + 1, parameters.end());

    T nll = quadra::stats::ar1_stationary_nll(states, T(0.0), phi, T(0.8));
    for (std::size_t i = 0; i < y_m.size(); ++i) {
      nll += quadra::stats::normal_nll(T(y_m[i]), states[i], T(0.5));
    }
    context.report("phi", phi);
    return nll;
  }

private:
  std::vector<double> y_m;
  quadra::ParameterSet parameters_m;
};

int main(int argc, char **argv) {
  LatentAr1Model model({-0.1, 0.6, 1.0, 0.5});

  quadra::LaplaceObjectiveOptions options;
  options.include_constant_m = true;
  options.newton_m.gradient_tolerance_m = 1e-11;
  options.newton_m.step_tolerance_m = 1e-13;

  quadra::stats::LaplaceEvaluator<LatentAr1Model> objective_evaluator(
      model, std::vector<double>(4, 0.0), model.parameters(), options);
  const auto discovery = objective_evaluator.evaluate({0.4});

  quadra::stats::ExactLaplaceEvaluator<LatentAr1Model> exact_evaluator(
      model, {0.4}, discovery.u_hat_m, model.parameters(), options);
  const auto exact = exact_evaluator.evaluate({0.4}, discovery.u_hat_m);

  if (!discovery.converged_m || !discovery.logdet_ok_m || !exact.success) {
    std::cerr << "Quadra Laplace evaluation failed: " << discovery.message_m
              << "\n";
    return 1;
  }
  if (argc > 1) {
    std::vector<std::string> names;
    for (std::size_t i = 0; i < exact.objective.u_hat_m.size(); ++i)
      names.push_back("x_" + std::to_string(i));
    comparison::write_quadra_diagnostics(
        argv[1], "Latent AR(1) Laplace Diagnostics",
        "Quadra diagnostics at the objective/gradient parity evaluation point.",
        exact.objective.laplace_objective_m, exact.gradient,
        exact.objective.newton_iterations_m, exact.objective.converged_m,
        exact.objective.message_m, exact.objective.hessian_random_m,
        exact.objective.u_hat_m, names, 1);
  }

  std::cout << std::setprecision(17);
  std::cout << "engine,objective,joint,logdet,gradient,"
               "u0,u1,u2,u3\n";
  std::cout << "quadra," << exact.objective.laplace_objective_m << ","
            << exact.objective.joint_objective_m << ","
            << exact.objective.log_det_hessian_m << ","
            << exact.gradient[0];
  for (double value : exact.objective.u_hat_m) {
    std::cout << "," << value;
  }
  std::cout << "\n";
  return 0;
}
