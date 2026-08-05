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

class Ma1ParityModel : public quadra::QuadraModel<Ma1ParityModel> {
public:
  Ma1ParityModel(std::vector<double> innovations, std::vector<double> y,
                 double mean, double innovation_sd, double observation_sd)
      : innovations_m(std::move(innovations)), y_m(std::move(y)),
        mean_m(mean), innovation_sd_m(innovation_sd),
        observation_sd_m(observation_sd) {}

  std::vector<std::string> parameter_names_impl() const {
    return {"unconstrained_theta"};
  }

  template <typename T>
  T evaluate_impl(const std::vector<T> &parameters,
                  quadra::ModelReportContext &context) const {
    std::vector<T> innovations;
    innovations.reserve(innovations_m.size());
    for (double value : innovations_m) {
      innovations.push_back(T(value));
    }

    const T theta =
        quadra::stats::correlation_from_unconstrained(parameters[0]);
    const std::vector<T> process =
        quadra::stats::ma1_from_innovations(innovations, T(mean_m), theta);

    T nll = quadra::stats::ma1_innovations_nll(
        innovations, T(innovation_sd_m));
    for (std::size_t i = 0; i < y_m.size(); ++i) {
      nll += quadra::stats::normal_nll(T(y_m[i]), process[i],
                                      T(observation_sd_m));
    }

    context.report("theta", theta);
    return nll;
  }

private:
  std::vector<double> innovations_m;
  std::vector<double> y_m;
  double mean_m;
  double innovation_sd_m;
  double observation_sd_m;
};

int main(int argc, char **argv) {
  Ma1ParityModel model({-0.3, 0.2, 0.8, -0.1, 0.4},
                       {0.7, 1.0, 0.1, 0.9}, 0.3, 0.8, 0.5);
  const auto result = quadra::evaluate_gradient(model, {-0.6});
  if (argc > 1) {
    const auto curvature = had::evaluate_value_gradient_hessian(
        [&](const std::vector<had::AReal> &parameters) {
          quadra::ModelReportContext context;
          model.initialize(context);
          return model.evaluate<had::AReal>(parameters, context);
        },
        {-0.6});
    Eigen::MatrixXd hessian(1, 1);
    hessian(0, 0) = curvature.hessian[0][0];
    comparison::write_quadra_fixed_diagnostics(
        argv[1], "Observed MA(1) Fixed-Effect Diagnostics",
        result.objective_value_m, result.gradient_m, hessian,
        {"unconstrained_theta"});
  }

  std::cout << std::setprecision(17);
  std::cout << "engine,objective,gradient,theta\n";
  std::cout << "quadra," << result.objective_value_m << ","
            << result.gradient_m[0] << "," << result.reports_m[0].value_m
            << "\n";
  return 0;
}
