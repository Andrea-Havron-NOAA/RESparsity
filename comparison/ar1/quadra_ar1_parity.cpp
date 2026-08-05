#include "core/autodiff/model_gradient.hpp"
#include "core/model/quadra_model.hpp"
#include "include/quadra/stats.hpp"
#include "../quadra_diagnostics.hpp"

#include <iomanip>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

DECLARE_ADGRAPH();

class Ar1ParityModel : public quadra::QuadraModel<Ar1ParityModel> {
public:
  Ar1ParityModel(std::vector<double> observations, double mean,
                 double innovation_sd)
      : observations_m(std::move(observations)), mean_m(mean),
        innovation_sd_m(innovation_sd) {}

  std::vector<std::string> parameter_names_impl() const {
    return {"unconstrained_phi"};
  }

  template <typename T>
  T evaluate_impl(const std::vector<T> &parameters,
                  quadra::ModelReportContext &context) const {
    std::vector<T> observations;
    observations.reserve(observations_m.size());
    for (double value : observations_m) {
      observations.push_back(T(value));
    }

    const T phi =
        quadra::stats::correlation_from_unconstrained(parameters[0]);
    context.report("phi", phi);
    return quadra::stats::ar1_stationary_nll(
        observations, T(mean_m), phi, T(innovation_sd_m));
  }

private:
  std::vector<double> observations_m;
  double mean_m;
  double innovation_sd_m;
};

int main(int argc, char **argv) {
  Ar1ParityModel model({-0.2, 0.4, 1.1, 0.7}, 0.3, 0.8);
  const auto result = quadra::evaluate_gradient(model, {0.4});
  if (argc > 1) {
    const auto curvature = had::evaluate_value_gradient_hessian(
        [&](const std::vector<had::AReal> &parameters) {
          quadra::ModelReportContext context;
          model.initialize(context);
          return model.evaluate<had::AReal>(parameters, context);
        },
        {0.4});
    Eigen::MatrixXd hessian(1, 1);
    hessian(0, 0) = curvature.hessian[0][0];
    comparison::write_quadra_fixed_diagnostics(
        argv[1], "Observed AR(1) Fixed-Effect Diagnostics",
        result.objective_value_m, result.gradient_m, hessian,
        {"unconstrained_phi"});
  }

  std::cout << std::setprecision(17);
  std::cout << "engine,objective,gradient,phi\n";
  std::cout << "quadra," << result.objective_value_m << ","
            << result.gradient_m[0] << "," << result.reports_m[0].value_m
            << "\n";
  return 0;
}
