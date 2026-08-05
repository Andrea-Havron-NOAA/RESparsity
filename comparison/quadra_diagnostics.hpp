#pragma once

#include "core/diagnostics/functional_analysis_markdown.hpp"
#include "core/laplace/functional_analysis_report.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <string>
#include <vector>

namespace comparison {

inline void replace_all(std::string &text, const std::string &from,
                        const std::string &to) {
  std::size_t position = 0;
  while ((position = text.find(from, position)) != std::string::npos) {
    text.replace(position, from.size(), to);
    position += to.size();
  }
}

inline void clarify_parity_point_markdown(const std::string &path) {
  std::ifstream input(path);
  const std::string original((std::istreambuf_iterator<char>(input)),
                             std::istreambuf_iterator<char>());
  std::string clarified = original;
  replace_all(clarified, "**Optimization quality:**",
              "**Parity-point quality:**");
  replace_all(clarified, "**Optimization:** converged",
              "**Latent mode:** converged");
  replace_all(clarified, "| Optimization |", "| Latent mode |");
  replace_all(clarified, "gradient norm =",
              "marginal fixed-effect gradient norm =");
  replace_all(clarified, "## Optimization", "## Parity-Point Evaluation");
  replace_all(clarified, "- Gradient norm:",
              "- Marginal fixed-effect gradient norm:");
  replace_all(clarified, "- Converged:", "- Latent mode converged:");
  replace_all(clarified, "- Max gradient parameter:",
              "- Maximum marginal-gradient parameter:");
  std::ofstream output(path);
  output << clarified;
}

inline void write_quadra_fixed_diagnostics(
    const std::string &path, const std::string &title, double objective,
    const std::vector<double> &gradient, const Eigen::MatrixXd &hessian,
    const std::vector<std::string> &parameter_names) {
  const auto curvature =
      quadra::summarize_laplace_hessian_structure(hessian);
  const auto spectrum = quadra::summarize_spectral_structure(hessian);
  const auto backend = quadra::laplace::analyze_hessian_structure(
      hessian.sparseView(1.0e-8, 1.0));
  double gradient_norm = 0.0;
  for (double value : gradient)
    gradient_norm += value * value;
  gradient_norm = std::sqrt(gradient_norm);

  std::ofstream md(path);
  md << std::setprecision(15);
  md << "# " << title << "\n\n";
  md << "Quadra fixed-effect curvature diagnostics at the parity evaluation "
        "point. This model has no latent Hessian; the spectrum below is the "
        "complete fixed-effect objective Hessian spectrum.\n\n";
  md << "## Evaluation\n\n";
  md << "- Objective: `" << objective << "`\n";
  md << "- Gradient norm: `" << gradient_norm << "`\n";
  md << "- Parameter count: `" << parameter_names.size() << "`\n\n";
  md << "## Curvature\n\n";
  md << "- Positive definite: `"
     << (curvature.positive_definite ? "yes" : "no") << "`\n";
  md << "- Minimum eigenvalue: `" << curvature.min_eigenvalue << "`\n";
  md << "- Maximum eigenvalue: `" << curvature.max_eigenvalue << "`\n";
  md << "- Condition number: `" << curvature.condition_number_abs << "`\n";
  md << "- Hessian inertia (positive / near-zero / negative): `"
     << spectrum.positive_eigen_count << " / " << spectrum.near_zero_eigen_count
     << " / " << spectrum.negative_eigen_count << "`\n";
  md << "- Normalized spectral entropy: `"
     << spectrum.normalized_spectral_entropy << "`\n";
  md << "- Participation ratio: `" << spectrum.participation_ratio << "`\n";
  md << "- Stable rank: `" << spectrum.stable_rank << "`\n\n";
  md << "## Quadra Backend Selection\n\n";
  md << "| Quantity | Selection |\n";
  md << "|---|---|\n";
  md << "| Detected Hessian structure | `"
     << quadra::laplace::ToString(backend.structure) << "` |\n";
  md << "| Factorization backend | `"
     << quadra::laplace::ToString(backend.backend) << "` |\n";
  md << "| Solver recommendation | `"
     << quadra::laplace::ToString(backend.solver) << "` |\n";
  md << "| Bandwidth | `" << backend.bandwidth << "` |\n";
  md << "| Expected complexity | `" << backend.complexity << "` |\n";
  md << "| Symbolic reuse supported | `"
     << (backend.supports_symbolic_reuse ? "yes" : "no") << "` |\n";
  md << "| Selection reason | " << backend.backend_reason << " |\n\n";
  md << "## Full Eigenvalue Spectrum\n\n";
  md << "| Rank | Eigenvalue | Cumulative positive-curvature share |\n";
  md << "|---:|---:|---:|\n";
  for (std::size_t i = 0; i < spectrum.signed_eigenvalues_desc.size(); ++i) {
    md << "| " << (i + 1) << " | `" << spectrum.signed_eigenvalues_desc[i]
       << "` | `"
       << (i < spectrum.cumulative_share.size()
               ? std::to_string(spectrum.cumulative_share[i])
               : std::string("n/a"))
       << "` |\n";
  }
}

inline void write_quadra_diagnostics(
    const std::string &base_path, const std::string &title,
    const std::string &subtitle, double objective,
    const std::vector<double> &gradient, int iterations, bool converged,
    const std::string &message, const Eigen::SparseMatrix<double> &hessian,
    const std::vector<double> &latent_states,
    const std::vector<std::string> &random_effect_names,
    std::size_t fixed_effect_count, bool parity_point = false) {
  quadra::FunctionalOptimizationSummary optimization;
  optimization.objective_value = objective;
  optimization.iterations = iterations;
  optimization.converged = converged;
  optimization.message = message;
  double squared_norm = 0.0;
  for (std::size_t i = 0; i < gradient.size(); ++i) {
    squared_norm += gradient[i] * gradient[i];
    if (i == 0 || std::abs(gradient[i]) > optimization.max_abs_gradient) {
      optimization.max_abs_gradient = std::abs(gradient[i]);
      optimization.max_gradient_value = gradient[i];
      optimization.max_gradient_parameter = "fixed_" + std::to_string(i);
    }
  }
  optimization.gradient_norm = std::sqrt(squared_norm);

  const Eigen::MatrixXd dense_hessian(hessian);
  const auto report = quadra::make_functional_analysis_report(
      optimization, dense_hessian, latent_states, 1.0e-8,
      random_effect_names);
  const std::string csv_path = base_path + ".csv";
  const std::string text_path = base_path + ".txt";
  const std::string markdown_path = base_path + ".md";
  quadra::write_functional_analysis_report_csv(report, csv_path);
  quadra::write_functional_analysis_report_text(report, text_path);

  std::string effective_entries_95 = "unknown";
  for (const auto &row : report.laplace_structure.effective_sparsity) {
    if (row.label == "95%")
      effective_entries_95 = std::to_string(row.entries_required);
  }
  std::string effective_bandwidth_95 = "unknown";
  for (const auto &row : report.laplace_structure.effective_bandwidth) {
    if (row.label == "95%")
      effective_bandwidth_95 = std::to_string(row.bandwidth);
  }

  quadra::diagnostics::MarkdownReportConfig config;
  config.title = title;
  config.subtitle = subtitle;
  config.output_path = markdown_path;
  config.functional_csv_path = csv_path;
  config.structure_txt_path = text_path;
  config.fixed_effects = std::to_string(fixed_effect_count);
  config.total_estimated =
      std::to_string(fixed_effect_count + latent_states.size());
  config.effective_entries_95 = effective_entries_95;
  config.effective_bandwidth_95 = effective_bandwidth_95;
  quadra::diagnostics::write_markdown_report(config);
  if (parity_point)
    clarify_parity_point_markdown(markdown_path);
}

} // namespace comparison
