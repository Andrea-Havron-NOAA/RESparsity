#include "core/laplace/model_analysis_report.hpp"
#include "core/model/quadra_model.hpp"
#include "include/quadra/stats.hpp"
#include "../quadra_diagnostics.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

DECLARE_ADGRAPH();

struct Observation {
  int year = 0;
  int fleet = 0;
  int age = 0;
  double value = 0.0;
};

std::vector<Observation> read_observations(const std::string &path) {
  std::ifstream input(path);
  if (!input)
    throw std::runtime_error("cannot open fisheries observations: " + path);
  std::string line;
  std::getline(input, line);
  std::vector<Observation> observations;
  while (std::getline(input, line)) {
    std::stringstream fields(line);
    std::string field;
    Observation observation;
    std::getline(fields, field, ','); observation.year = std::stoi(field);
    std::getline(fields, field, ','); observation.fleet = std::stoi(field);
    std::getline(fields, field, ','); observation.age = std::stoi(field);
    std::getline(fields, field, ','); observation.value = std::stod(field);
    observations.push_back(observation);
  }
  return observations;
}

class SimpleFsaArGateModel
    : public quadra::QuadraModel<SimpleFsaArGateModel> {
public:
  explicit SimpleFsaArGateModel(std::vector<Observation> observations,
                                bool expanded_nuisance = false,
                                bool full_fixed = false,
                                std::string process = "ar")
      : observations_(std::move(observations)),
        expanded_nuisance_(expanded_nuisance), full_fixed_(full_fixed),
        process_(std::move(process)) {
    parameters_.add(process_ == "ar" ? "tphiR" : "tthetaR", 0.0,
                    quadra::ParameterTransform::Identity, false);
    parameters_.add("logSdR", 0.0, quadra::ParameterTransform::Identity, false);
    if (expanded_nuisance_) {
      parameters_.add("logMuR", 0.0, quadra::ParameterTransform::Identity,
                      false);
      parameters_.add("logSdCatch", 0.0,
                      quadra::ParameterTransform::Identity, false);
      parameters_.add("logSdSurvey", 0.0,
                      quadra::ParameterTransform::Identity, false);
      for (int age = 1; age <= 5; ++age)
        parameters_.add("logQ_" + std::to_string(age), 0.0,
                        quadra::ParameterTransform::Identity, false);
    }
    if (full_fixed_) {
      for (int age = 1; age <= 7; ++age)
        parameters_.add("logN1Y_" + std::to_string(age), 0.0,
                        quadra::ParameterTransform::Identity, false);
      for (int year = 1; year <= 45; ++year)
        parameters_.add("logFY_" + std::to_string(year), 0.0,
                        quadra::ParameterTransform::Identity, false);
      for (int age = 1; age <= 4; ++age)
        parameters_.add("logFA_" + std::to_string(age), 0.0,
                        quadra::ParameterTransform::Identity, false);
    }
    const int random_count = process_ == "ar" ? 44 : 45;
    const std::string random_prefix = process_ == "ma_good" ? "eps_" : "logN1A_";
    for (int index = 0; index < random_count; ++index)
      parameters_.add(random_prefix + std::to_string(index + 1), 0.0,
                      quadra::ParameterTransform::Identity, true);
  }

  std::vector<std::string> parameter_names_impl() const {
    return parameters_.names();
  }
  const quadra::ParameterSet &parameters() const { return parameters_; }
  std::size_t fixed_size() const {
    return full_fixed_ ? 66u : (expanded_nuisance_ ? 10u : 2u);
  }
  std::size_t random_size() const { return process_ == "ar" ? 44u : 45u; }

  template <typename T>
  T evaluate_impl(const std::vector<T> &parameters,
                  quadra::ModelReportContext &) const {
    const T phi = quadra::stats::correlation_from_unconstrained(parameters[0]);
    const T recruitment_sd = exp(parameters[1]);
    const T conditional_sd = recruitment_sd * sqrt(T(1.0) - phi * phi);
    const T recruitment_mean = expanded_nuisance_ ? parameters[2] : T(0.0);
    const T catch_sd = expanded_nuisance_ ? exp(parameters[3]) : T(1.0);
    const T survey_sd = expanded_nuisance_ ? exp(parameters[4]) : T(1.0);
    std::vector<T> latent(parameters.begin() + fixed_size(), parameters.end());
    std::vector<T> recruits(44, T(0.0));
    T nll = T(0.0);
    if (process_ == "ar") {
      nll = quadra::stats::normal_nll(T(0.0), recruitment_mean,
                                      recruitment_sd);
      T previous = T(0.0);
      recruits = latent;
      for (const T &recruit : recruits) {
        const T conditional_mean = recruitment_mean +
                                   phi * (previous - recruitment_mean);
        nll += quadra::stats::normal_nll(recruit, conditional_mean,
                                         conditional_sd);
        previous = recruit;
      }
    } else if (process_ == "ma_good") {
      for (const T &innovation : latent)
        nll += quadra::stats::normal_nll(innovation, T(0.0), recruitment_sd);
      for (int year = 0; year < 44; ++year)
        recruits[year] = recruitment_mean + latent[year + 1] + phi * latent[year];
    } else {
      std::vector<T> innovations(45, T(0.0));
      innovations[0] = latent[0];
      nll += quadra::stats::normal_nll(innovations[0], T(0.0), recruitment_sd);
      for (int i = 1; i < 45; ++i) {
        innovations[i] = latent[i] - phi * innovations[i - 1];
        nll += quadra::stats::normal_nll(innovations[i], T(0.0), recruitment_sd);
      }
      for (int year = 0; year < 44; ++year)
        recruits[year] = latent[year + 1] + recruitment_mean;
    }

    std::vector<T> log_n(7 * 45, T(0.0));
    const double mortality[7] = {0.8, 0.35, 0.25, 0.2, 0.2, 0.2, 0.2};
    if (full_fixed_) {
      for (int age = 0; age < 7; ++age)
        log_n[age] = parameters[10 + age];
    }
    for (int year = 1; year < 45; ++year) {
      log_n[year * 7] = recruits[year - 1];
      for (int age = 1; age < 7; ++age) {
        const T fishing = full_fixed_
                              ? exp(parameters[17 + year - 1] +
                                    (age - 1 < 4
                                         ? parameters[62 + age - 1]
                                         : T(0.0)))
                              : T(1.0);
        log_n[year * 7 + age] =
            log_n[(year - 1) * 7 + age - 1] -
            fishing - T(mortality[age - 1]);
      }
    }

    for (const Observation &observation : observations_) {
      const int year = observation.year - 1963;
      const int age = observation.age - 1;
      T prediction = log_n[year * 7 + age];
      const T fishing = full_fixed_
                            ? exp(parameters[17 + year] +
                                  (age < 4 ? parameters[62 + age] : T(0.0)))
                            : T(1.0);
      if (observation.fleet == 1) {
        const T total_mortality = fishing + T(mortality[age]);
        prediction += log(fishing) - log(total_mortality) +
                      log(T(1.0) - exp(-total_mortality));
      } else {
        const T log_q = expanded_nuisance_ && age < 5
                            ? parameters[5 + age]
                            : T(0.0);
        prediction += log_q - (fishing + T(mortality[age])) * T(0.125);
      }
      nll += quadra::stats::normal_nll(T(std::log(observation.value)),
                                       prediction,
                                       observation.fleet == 1 ? catch_sd
                                                              : survey_sd);
    }
    return nll;
  }

private:
  std::vector<Observation> observations_;
  bool expanded_nuisance_ = false;
  bool full_fixed_ = false;
  std::string process_ = "ar";
  quadra::ParameterSet parameters_;
};

int main(int argc, char **argv) {
  if (argc < 2 || argc > 5) {
    std::cerr << "usage: quadra_simplefsa_ar_gate observations.csv "
                 "[diagnostics_base]\n";
    return 1;
  }
  const std::string gate = argc >= 4 ? std::string(argv[3]) : "core";
  const bool ma_model = gate == "ma_good" || gate == "ma_bad";
  const bool full_fixed = gate == "full" || ma_model;
  const bool expanded_nuisance = gate == "nuisance" || full_fixed;
  const std::string process = gate == "ma_good" ? "ma_good" :
                              (gate == "ma_bad" ? "ma_bad" : "ar");
  SimpleFsaArGateModel model(read_observations(argv[1]), expanded_nuisance,
                             full_fixed, process);
  quadra::LaplaceObjectiveOptions objective_options;
  objective_options.include_constant_m = true;
  objective_options.newton_m.gradient_tolerance_m = 1e-8;
  objective_options.newton_m.step_tolerance_m = 1e-11;
  std::vector<double> random_start(model.random_size(), 0.0);
  const std::vector<double> zero_fixed(model.fixed_size(), 0.0);
  std::vector<double> fixed_start = zero_fixed;
  if (argc == 5) {
    std::ifstream start_input(argv[4]);
    std::string value;
    std::size_t index = 0;
    while (std::getline(start_input, value) && index < fixed_start.size())
      fixed_start[index++] = std::stod(value);
    if (index != fixed_start.size())
      throw std::runtime_error("fixed-start file has wrong length");
  }
  quadra::stats::ExactLaplaceEvaluator<SimpleFsaArGateModel> evaluator(
      model, zero_fixed, random_start, model.parameters(), objective_options);
  const auto initial = evaluator.evaluate(zero_fixed, random_start);
  if (!initial.success)
    throw std::runtime_error("initial Quadra fisheries evaluation failed");
  quadra::stats::ExactLaplaceResult ma_probe;
  if (ma_model) {
    std::vector<double> probe_fixed = zero_fixed;
    probe_fixed[0] = 0.75;
    ma_probe = evaluator.evaluate(probe_fixed, random_start);
    if (!ma_probe.success)
      throw std::runtime_error("off-zero Quadra fisheries MA probe failed");
  }

  quadra::stats::LaplaceOptimizerOptions optimizer_options;
  optimizer_options.max_iterations = full_fixed ? 500 :
                                     (expanded_nuisance ? 200 : 100);
  optimizer_options.gradient_tolerance = 1e-4;
  quadra::stats::LaplaceOptimizerResult fit;
  if (ma_model && argc != 5) {
    bool have_fit = false;
    for (double dependence_start : {0.0, 0.75, -0.75}) {
      std::vector<double> candidate_start = fixed_start;
      candidate_start[0] = dependence_start;
      (void)evaluator.evaluate(candidate_start, random_start);
      auto candidate = quadra::stats::optimize_laplace(
          evaluator, candidate_start, optimizer_options);
      if (candidate.converged &&
          (!have_fit || candidate.objective < fit.objective)) {
        fit = std::move(candidate);
        have_fit = true;
      }
    }
    if (!have_fit)
      throw std::runtime_error("Quadra fisheries MA multistart failed");
  } else {
    fit = quadra::stats::optimize_laplace(evaluator, fixed_start,
                                           optimizer_options);
  }
  if (!fit.converged) {
    for (const auto &iteration : fit.history) {
      std::cerr << "iteration=" << iteration.iteration
                << " objective=" << iteration.objective
                << " gradient_norm=" << iteration.gradient_norm
                << " step_scale=" << iteration.step_scale
                << " step_norm=" << iteration.step_norm << "\n";
    }
    throw std::runtime_error("Quadra fisheries optimization failed: " +
                             fit.message);
  }
  const auto final = evaluator.evaluate(fit.fixed, fit.random_mode);
  const auto diagnostics = quadra::laplace::analyze_hessian_structure(
      final.objective.hessian_random_m);
  if (argc >= 3) {
    std::vector<std::string> random_names;
    for (std::size_t index = 0; index < model.random_size(); ++index)
      random_names.push_back((process == "ma_good" ? "eps_" : "logN1A_") +
                             std::to_string(index + 1));
    comparison::write_quadra_diagnostics(
        argv[2],
        ma_model
            ? (process == "ma_good"
                   ? "Simple FSA Complete MA Innovation Diagnostics"
                   : "Simple FSA Complete MA Inverse-Filter Diagnostics")
            : full_fixed
            ? "Simple FSA Complete AR Fixed-Block Diagnostics"
            : expanded_nuisance
            ? "Simple FSA AR Observation/Recruitment Nuisance Diagnostics"
            : "Simple FSA AR Correlation/Scale Diagnostics",
        ma_model
            ? "Quadra diagnostics for all 66 mapped fixed effects and 45 "
              "latent MA effects using all 440 observations."
            : full_fixed
            ? "Quadra diagnostics for all 66 mapped fixed effects and 44 "
              "latent recruits using all 440 observations."
            : expanded_nuisance
            ? "Quadra diagnostics for the 10-fixed-effect nuisance block and "
              "44 latent recruits using all 440 observations."
            : "Quadra diagnostics at the optimized Laplace solution using all "
              "440 observations.",
        final.objective.laplace_objective_m, final.gradient, fit.iterations,
        fit.converged, fit.message, final.objective.hessian_random_m,
        final.objective.u_hat_m, random_names, model.fixed_size());
  }

  std::cout << std::setprecision(17);
  if (full_fixed) {
    double initial_gradient_norm = 0.0;
    double fit_gradient_norm = 0.0;
    for (double value : initial.gradient)
      initial_gradient_norm += value * value;
    for (double value : final.gradient)
      fit_gradient_norm += value * value;
    std::cout << "engine,objective,gradient_norm,fit_objective,"
                 "fit_gradient_norm,optimization_converged,structure,backend,"
                 "fixed_effects,random_effects,active_directions,hdot_workers\n";
    std::cout << "quadra," << initial.objective.laplace_objective_m << ","
              << std::sqrt(initial_gradient_norm) << ","
              << final.objective.laplace_objective_m << ","
              << std::sqrt(fit_gradient_norm) << ",1,"
              << quadra::laplace::ToString(diagnostics.structure) << ","
              << quadra::laplace::ToString(diagnostics.backend) << ",66,"
              << model.random_size() << ","
              << initial.active_directions.size() << ","
              << evaluator.hdot_worker_count() << "\n";
    std::cout << "parameter,name,initial_gradient,estimate\n";
    const auto names = model.parameters().names();
    for (std::size_t i = 0; i < model.fixed_size(); ++i)
      std::cout << "parameter," << names[i] << "," << initial.gradient[i]
                << "," << fit.fixed[i] << "\n";
    if (ma_model)
      std::cout << "probe,0.75," << ma_probe.objective.laplace_objective_m
                << "," << ma_probe.gradient[0] << "\n";
    return 0;
  }
  if (expanded_nuisance) {
    std::cout << "engine,objective,gradient_phi,gradient_log_sd,gradient_mu,"
                 "gradient_catch_sd,gradient_survey_sd,gradient_q1,gradient_q2,"
                 "gradient_q3,gradient_q4,gradient_q5,estimate_phi,"
                 "estimate_log_sd,estimate_mu,estimate_catch_sd,"
                 "estimate_survey_sd,estimate_q1,estimate_q2,estimate_q3,"
                 "estimate_q4,estimate_q5,fit_objective,fit_gradient_norm,"
                 "optimization_converged,structure,backend,random_effects,"
                 "active_directions,hdot_workers\n";
    double fit_gradient_norm = 0.0;
    for (double value : final.gradient)
      fit_gradient_norm += value * value;
    std::cout << "quadra," << initial.objective.laplace_objective_m;
    for (int i = 0; i < 10; ++i)
      std::cout << "," << initial.gradient[static_cast<std::size_t>(i)];
    for (int i = 0; i < 10; ++i)
      std::cout << "," << fit.fixed[static_cast<std::size_t>(i)];
    std::cout << "," << final.objective.laplace_objective_m << ","
              << std::sqrt(fit_gradient_norm) << ",1,"
              << quadra::laplace::ToString(diagnostics.structure) << ","
              << quadra::laplace::ToString(diagnostics.backend) << ",44,"
              << initial.active_directions.size() << ","
              << evaluator.hdot_worker_count() << "\n";
    return 0;
  }
  std::cout << "engine,objective,gradient_phi,gradient_log_sd,estimate_phi,"
               "estimate_log_sd,fit_objective,fit_gradient_phi,"
               "fit_gradient_log_sd,"
               "optimization_converged,structure,backend,random_effects,"
               "active_directions,hdot_workers\n";
  std::cout << "quadra," << initial.objective.laplace_objective_m << ","
            << initial.gradient[0] << "," << initial.gradient[1] << ","
            << fit.fixed[0] << "," << fit.fixed[1] << ","
            << final.objective.laplace_objective_m << ","
            << final.gradient[0] << "," << final.gradient[1] << ",1,"
            << quadra::laplace::ToString(diagnostics.structure) << ","
            << quadra::laplace::ToString(diagnostics.backend) << ",44,"
            << initial.active_directions.size() << ","
            << evaluator.hdot_worker_count() << "\n";
  return 0;
}
