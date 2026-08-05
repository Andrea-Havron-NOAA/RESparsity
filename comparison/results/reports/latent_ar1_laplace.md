# Latent AR(1) Laplace Diagnostics

Quadra diagnostics at the objective/gradient parity evaluation point.

## Executive Summary

- **Overall status:** `REVIEW`.
- **Confidence:** `LOW`.
- **Optimization quality:** `REVIEW`.
- **Uncertainty structure:** `LOCAL`.
- **Optimization:** converged = `yes`, gradient norm = `0.231236105000721`.
- **Curvature health:** positive definite = `yes`, condition number = `1.19586075231367`.
- **Quadra factorization:** structure = `tridiagonal`, backend = `tridiagonal`.
- **Latent structure:** `4` random effects were estimated.
- **Symbolic vs numerical structure:** structural density = `0.625`, but 95% of curvature is retained by `7` entries.
- **Spectral complexity:** entropy effective rank = `3.99083541706421`, with 90% curvature requiring `4` eigen-directions.

## Model Health Assessment

| Check | Status | Evidence |
|---|---:|---|
| Optimization | `PASS` | converged = `yes` |
| Gradient quality | `CHECK` | gradient norm = `0.231236105000721` |
| Curvature | `PASS` | positive definite = `yes` |
| Conditioning | `EXCELLENT` | condition number = `1.19586075231367` |
| Overall status | `REVIEW` | rule-based v1 diagnostic |
| Confidence | `LOW` | based on convergence, gradient, PD status, and conditioning |

**Interpretation:** the rule-based health check is intentionally simple. It flags obvious numerical issues quickly, but it does not replace scientific review or model-specific diagnostics.

## Model Complexity

| Quantity | Value |
|---|---:|
| Fixed effects | `1` |
| Random effects | `4` |
| Total estimated quantities | `5` |
| Structural nonzeros | `10` |
| Structural density | `0.625` |
| Entries for 95% curvature | `7` |
| Effective bandwidth for 95% curvature | `0` |
| 95% curvature compression | `1.42857x` |

## Optimization

- Quality: `REVIEW`
- Objective value: `4.22118997734713`
- Gradient norm: `0.231236105000721`
- Converged: `yes`
- Max gradient parameter: `fixed_0`

## Curvature

- Positive definite: `yes`
- Condition number: `1.19586075231367`
- Minimum eigenvalue: `5.10642935602523`
- Maximum eigenvalue: `6.10657845133293`

## Quadra Backend Selection

| Quantity | Selection |
|---|---|
| Detected Hessian structure | `tridiagonal` |
| Factorization backend | `tridiagonal` |
| Random-effect solver | `Newton` |
| Bandwidth | `1` |
| Expected complexity | `O(n)` |
| Symbolic reuse supported | `no` |
| Selection reason | unit bandwidth |

## Spectral Structure

- Largest eigenvalue share: `0.272959468785324`
- Hessian inertia (positive / near-zero / negative): `4 / 0 / 0`
- Entropy effective rank: `3.99083541706421`
- Normalized spectral entropy: `0.998345391506612`
- Participation ratio: `3.98172240278072`
- Stable rank: `3.370797729766`
- Leading spectral-gap ratio: `1.05814138997149`
- Eigenvectors needed for 90% curvature: `4`
- Eigenvectors needed for 95% curvature: `4`

**Interpretation:** curvature is distributed across many latent-state directions rather than being dominated by one or two modes. That is a good sign for numerical stability.

### Full Eigenvalue Spectrum

| Rank | Eigenvalue | Cumulative curvature share |
|---:|---:|---:|
| `rank_1` | `6.10657845133293` | `0.272959468785324` |
| `rank_2` | `5.7710420452388` | `0.530920712313164` |
| `rank_3` | `5.38769082563393` | `0.771746444343774` |
| `rank_4` | `5.10642935602523` | `1` |

## Effective Structure

- Structural density: `0.625`
- Structural nonzeros: `10`
- Entries for 95% curvature: `7`
- Effective bandwidth for 95% curvature: `0`
- 95% curvature compression: `1.42857x`

**Interpretation:** symbolic density alone overstates practical complexity. The detailed Laplace report below shows that large amounts of curvature can be retained with far fewer entries or a narrow effective bandwidth.

## Correlation Graph

- Classification: `LOCAL`
- Average degree: `0`
- Maximum degree: `0`
- Connected components: `4`
- Largest component size: `1`
- Graph diameter: `0`

**Interpretation:** a LOCAL graph means the strongest uncertainty relationships are neighborhood-like rather than globally tangled.

## Latent State Summary

- Count: `4`
- Mean: `0.395080290790782`
- Standard deviation: `0.288113933246675`

## Key Takeaway

This report demonstrates why Quadra's functional analysis diagnostics are useful: a model can look dense from a symbolic Hessian pattern, while numerical curvature, graph structure, and effective bandwidth reveal a simpler local-dependence structure.

## Full Laplace Structure Report

```text
Functional Analysis Report
==========================

Optimization
------------
objective_value:            4.22118997734713
gradient_norm:              0.231236105000721
max_gradient_parameter:     fixed_0
max_gradient_value:         -0.231236105000721
max_abs_gradient:           0.231236105000721
iterations:                 0
converged:                  yes
message:                    Converged: gradient norm below tolerance.

Curvature
---------
positive_definite:          yes
min_eigenvalue:             5.10642935602523
max_eigenvalue:             6.10657845133293
condition_number_abs:       1.19586075231367

Backend Selection
-----------------
detected_structure:         tridiagonal
selected_backend:           tridiagonal
solver_recommendation:      Newton
bandwidth:                  1
expected_complexity:        O(n)
symbolic_reuse:             not required
selection_reason:           unit bandwidth

Spectral Structure
------------------
available:                  yes
eigen_count:                4
positive_eigen_count:       4
near_zero_eigen_count:      0
negative_eigen_count:       0
largest_eigen_share:        0.272959468785324
effective_rank_entropy:     3.99083541706421
normalized_entropy:         0.998345391506612
participation_ratio:        3.98172240278072
stable_rank:                3.370797729766
leading_gap_ratio:          1.05814138997149
eigen_count_for_50%:        2
eigen_count_for_90%:        4
eigen_count_for_95%:        4
eigen_count_for_99%:        4

rank,eigenvalue,cumulative_curvature_share
1,6.10657845133293,0.272959468785324
2,5.7710420452388,0.530920712313164
3,5.38769082563393,0.771746444343774
4,5.10642935602523,1

Huu Structure
-------------
random_effects:             4
total_entries:              16
structural_nonzeros:        10
structural_density:         0.625

Effective Sparsity
------------------
curvature_retained,entries_required,entry_share,compression_vs_structural
90%,4,0.25,2.5
95%,7,0.4375,1.42857142857143
97%,8,0.5,1.25
98%,9,0.5625,1.11111111111111
99%,10,0.625,1
99.5%,10,0.625,1
99.9%,10,0.625,1
100%,10,0.625,1

Effective Bandwidth
-------------------
curvature_retained,bandwidth,entry_count_if_banded,entry_share_if_banded
90%,0,4,0.25
95%,0,4,0.25
97%,1,10,0.625
98%,1,10,0.625
99%,1,10,0.625
99.5%,1,10,0.625
99.9%,1,10,0.625
100%,3,16,1

Uncertainty
-----------
covariance_available:       yes
correlation_available:      yes
covariance_size:            4 x 4
min_variance:               0.178913056082264
max_variance:               0.180325236631532
min_variance_index:         2
max_variance_index:         0
max_abs_correlation:        0.0552249856038601
max_abs_correlation_pair:   0,1
count_abs_corr_gt_0_5:      0
count_abs_corr_gt_0_8:      0
count_abs_corr_gt_0_9:      0

Parameter Influence
-------------------
available:                  yes
Top parameter importance
index,name,variance,sd,variance_share,correlation_centrality,correlation_centrality_share,curvature_column_norm,curvature_diagonal,importance_score,importance_share
1,x_1,0.178913056082265,0.422981153341688,0.249017239686087,0.113272520375858,0.329849803043014,5.64025830796596,5.62337033911544,2.65595981732747,0.257364514793428
2,x_2,0.178913056082264,0.422981153341688,0.249017239686087,0.113272520375858,0.329849803043014,5.64025830796596,5.62337033911544,2.65595981732747,0.257364514793428
0,x_0,0.180325236631532,0.424647190773154,0.250982760313913,0.0584306598759833,0.170150196956986,5.57104264522072,5.5625,2.5039586342494,0.242635485206572
3,x_3,0.180325236631532,0.424647190773154,0.250982760313913,0.0584306598759833,0.170150196956986,5.57104264522072,5.5625,2.5039586342494,0.242635485206572

Top correlation pairs
i,j,name_i,name_j,correlation,abs_correlation
0,1,x_0,x_1,0.0552249856038601,0.0552249856038601
2,3,x_2,x_3,0.0552249856038601,0.0552249856038601
1,2,x_1,x_2,0.0550096288127401,0.0550096288127401
1,3,x_1,x_3,0.00303790595925726,0.00303790595925726
0,2,x_0,x_2,0.00303790595925726,0.00303790595925726
0,3,x_0,x_3,0.000167768312865863,0.000167768312865863

Correlation Graph
-----------------
available:                  yes
abs_correlation_threshold:  0.5
node_count:                 4
edge_count:                 0
average_degree:             0
maximum_degree:             0
maximum_degree_parameter:   x_0
connected_components:       4
largest_component_size:     1
graph_diameter:             0

Parameter Geometry
------------------
available:                  no
dominant_parameter:
dominant_parameter_index:   0
dominant_curvature_norm:    nan
index,name,gradient,abs_gradient,curvature_column_norm,curvature_diagonal,curvature_share

Gradient Volatility
-------------------
available:                  no
perturbation_scale:         0
samples:                    0
baseline_gradient_norm:     nan
mean_gradient_norm:         nan
sd_gradient_norm:           nan
max_gradient_norm:          nan
gradient_norm_cv:           nan
most_volatile_parameter:
most_volatile_parameter_sd: nan
most_sign_flips_parameter:
most_sign_flips:            0

Latent States
-------------
count:                      4
mean:                       0.395080290790782
sd:                         0.288113933246675
min_value:                  -0.0460804099873798
max_value:                  0.758893393773622
min_index:                  0
max_index:                  2
l2_norm:                    0.97795311687667
```
