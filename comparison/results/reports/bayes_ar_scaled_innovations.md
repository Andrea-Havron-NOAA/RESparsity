# Bayesian AR Diagnostics: scaled_innovations

Quadra diagnostics at the parameterization parity evaluation point.

## Executive Summary

- **Overall status:** `REVIEW`.
- **Confidence:** `LOW`.
- **Optimization quality:** `REVIEW`.
- **Uncertainty structure:** `LOCAL`.
- **Optimization:** converged = `yes`, gradient norm = `0.472367612872141`.
- **Curvature health:** positive definite = `yes`, condition number = `1.57854369363476`.
- **Quadra factorization:** structure = `dense`, backend = `dense_ldlt`.
- **Latent structure:** `4` random effects were estimated.
- **Symbolic vs numerical structure:** structural density = `1`, but 95% of curvature is retained by `10` entries.
- **Spectral complexity:** entropy effective rank = `3.93970979377318`, with 90% curvature requiring `4` eigen-directions.

## Model Health Assessment

| Check | Status | Evidence |
|---|---:|---|
| Optimization | `PASS` | converged = `yes` |
| Gradient quality | `CHECK` | gradient norm = `0.472367612872141` |
| Curvature | `PASS` | positive definite = `yes` |
| Conditioning | `EXCELLENT` | condition number = `1.57854369363476` |
| Overall status | `REVIEW` | rule-based v1 diagnostic |
| Confidence | `LOW` | based on convergence, gradient, PD status, and conditioning |

**Interpretation:** the rule-based health check is intentionally simple. It flags obvious numerical issues quickly, but it does not replace scientific review or model-specific diagnostics.

## Model Complexity

| Quantity | Value |
|---|---:|
| Fixed effects | `1` |
| Random effects | `4` |
| Total estimated quantities | `5` |
| Structural nonzeros | `16` |
| Structural density | `1` |
| Entries for 95% curvature | `10` |
| Effective bandwidth for 95% curvature | `1` |
| 95% curvature compression | `1.6x` |

## Optimization

- Quality: `REVIEW`
- Objective value: `4.17136656947898`
- Gradient norm: `0.472367612872141`
- Converged: `yes`
- Max gradient parameter: `fixed_0`

## Curvature

- Positive definite: `yes`
- Condition number: `1.57854369363476`
- Minimum eigenvalue: `4.58117348135282`
- Maximum eigenvalue: `7.2315825084363`

## Quadra Backend Selection

| Quantity | Selection |
|---|---|
| Detected Hessian structure | `dense` |
| Factorization backend | `dense_ldlt` |
| Random-effect solver | `LBFGS` |
| Bandwidth | `3` |
| Expected complexity | `O(n^3)` |
| Symbolic reuse supported | `no` |
| Selection reason | small matrix; dense LDLT preferred |

## Spectral Structure

- Largest eigenvalue share: `0.314646280832335`
- Hessian inertia (positive / near-zero / negative): `4 / 0 / 0`
- Entropy effective rank: `3.93970979377318`
- Normalized spectral entropy: `0.989044681020492`
- Participation ratio: `3.88057844431002`
- Stable rank: `2.60290506008368`
- Leading spectral-gap ratio: `1.19357371363784`
- Eigenvectors needed for 90% curvature: `4`
- Eigenvectors needed for 95% curvature: `4`

**Interpretation:** curvature is distributed across many latent-state directions rather than being dominated by one or two modes. That is a good sign for numerical stability.

### Full Eigenvalue Spectrum

| Rank | Eigenvalue | Cumulative curvature share |
|---:|---:|---:|
| `rank_1` | `7.2315825084363` | `0.314646280832335` |
| `rank_2` | `6.05876488884417` | `0.578263246619339` |
| `rank_3` | `5.11169217351714` | `0.800673062075445` |
| `rank_4` | `4.58117348135282` | `1` |

## Effective Structure

- Structural density: `1`
- Structural nonzeros: `16`
- Entries for 95% curvature: `10`
- Effective bandwidth for 95% curvature: `1`
- 95% curvature compression: `1.6x`

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
- Mean: `0.341641021182953`
- Standard deviation: `0.252535811932874`

## Key Takeaway

This report demonstrates why Quadra's functional analysis diagnostics are useful: a model can look dense from a symbolic Hessian pattern, while numerical curvature, graph structure, and effective bandwidth reveal a simpler local-dependence structure.

## Full Laplace Structure Report

```text
Functional Analysis Report
==========================

Optimization
------------
objective_value:            4.17136656947898
gradient_norm:              0.472367612872141
max_gradient_parameter:     fixed_0
max_gradient_value:         -0.472367612872141
max_abs_gradient:           0.472367612872141
iterations:                 1
converged:                  yes
message:                    Converged: gradient norm below tolerance.

Curvature
---------
positive_definite:          yes
min_eigenvalue:             4.58117348135282
max_eigenvalue:             7.2315825084363
condition_number_abs:       1.57854369363476

Backend Selection
-----------------
detected_structure:         dense
selected_backend:           dense_ldlt
solver_recommendation:      LBFGS
bandwidth:                  3
expected_complexity:        O(n^3)
symbolic_reuse:             not required
selection_reason:           small matrix; dense LDLT preferred

Spectral Structure
------------------
available:                  yes
eigen_count:                4
positive_eigen_count:       4
near_zero_eigen_count:      0
negative_eigen_count:       0
largest_eigen_share:        0.314646280832335
effective_rank_entropy:     3.93970979377318
normalized_entropy:         0.989044681020492
participation_ratio:        3.88057844431002
stable_rank:                2.60290506008368
leading_gap_ratio:          1.19357371363784
eigen_count_for_50%:        2
eigen_count_for_90%:        4
eigen_count_for_95%:        4
eigen_count_for_99%:        4

rank,eigenvalue,cumulative_curvature_share
1,7.2315825084363,0.314646280832335
2,6.05876488884417,0.578263246619339
3,5.11169217351714,0.800673062075445
4,4.58117348135282,1

Huu Structure
-------------
random_effects:             4
total_entries:              16
structural_nonzeros:        16
structural_density:         1

Effective Sparsity
------------------
curvature_retained,entries_required,entry_share,compression_vs_structural
90%,8,0.5,2
95%,10,0.625,1.6
97%,10,0.625,1.6
98%,11,0.6875,1.45454545454545
99%,13,0.8125,1.23076923076923
99.5%,14,0.875,1.14285714285714
99.9%,16,1,1
100%,16,1,1

Effective Bandwidth
-------------------
curvature_retained,bandwidth,entry_count_if_banded,entry_share_if_banded
90%,1,10,0.625
95%,1,10,0.625
97%,1,10,0.625
98%,1,10,0.625
99%,2,14,0.875
99.5%,2,14,0.875
99.9%,3,16,1
100%,3,16,1

Uncertainty
-----------
covariance_available:       yes
correlation_available:      yes
covariance_size:            4 x 4
min_variance:               0.176335211134886
max_variance:               0.181235272782077
min_variance_index:         0
max_variance_index:         3
max_abs_correlation:        0.139414365876299
max_abs_correlation_pair:   0,1
count_abs_corr_gt_0_5:      0
count_abs_corr_gt_0_8:      0
count_abs_corr_gt_0_9:      0

Parameter Influence
-------------------
available:                  yes
Top parameter importance
index,name,variance,sd,variance_share,correlation_centrality,correlation_centrality_share,curvature_column_norm,curvature_diagonal,importance_score,importance_share
1,latent_1,0.179830476460299,0.42406423624293,0.250723191988205,0.284054858340175,0.331559680796725,5.90506467136517,5.78773645533904,3.21543610617514,0.267860467759505
2,latent_2,0.179846116106827,0.424082676027714,0.250744997091413,0.280809155453506,0.327771172410165,5.89491361825207,5.78166585863433,3.20193418258453,0.266735696049164
0,latent_0,0.176335211134886,0.419922863315259,0.245850024233312,0.147746089853997,0.17245487958629,5.84829702428651,5.78797294767827,2.81867344841001,0.23480833187835
3,latent_3,0.181235272782077,0.425717362556517,0.25268178668707,0.144113059004676,0.168214267206821,5.68318497620148,5.62583779049879,2.76810205213014,0.230595504312981

Top correlation pairs
i,j,name_i,name_j,correlation,abs_correlation
0,1,latent_0,latent_1,-0.139414365876299,0.139414365876299
1,2,latent_1,latent_2,-0.136945676810699,0.136945676810699
2,3,latent_2,latent_3,-0.135974999008304,0.135974999008304
0,2,latent_0,latent_2,-0.00788847963450345,0.00788847963450345
1,3,latent_1,latent_3,-0.00769481565317716,0.00769481565317716
0,3,latent_0,latent_3,-0.00044324434319495,0.00044324434319495

Correlation Graph
-----------------
available:                  yes
abs_correlation_threshold:  0.5
node_count:                 4
edge_count:                 0
average_degree:             0
maximum_degree:             0
maximum_degree_parameter:   latent_0
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
mean:                       0.341641021182953
sd:                         0.252535811932874
min_value:                  -0.0161058074976641
max_value:                  0.660557114494651
min_index:                  0
max_index:                  2
l2_norm:                    0.849689175318897
```
