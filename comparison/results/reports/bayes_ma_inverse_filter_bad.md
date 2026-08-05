# Bayesian MA Diagnostics: inverse_filter_bad

Quadra diagnostics at the parameterization parity evaluation point. The latent mode converged; the fixed effects were held fixed and were not optimized.

## Executive Summary

- **Overall status:** `REVIEW`.
- **Confidence:** `LOW`.
- **Parity-point quality:** `REVIEW`.
- **Uncertainty structure:** `LOCAL`.
- **Latent mode:** converged = `yes`, marginal fixed-effect gradient norm = `0.655796536381041`.
- **Curvature health:** positive definite = `yes`, condition number = `1.31843011968793`.
- **Quadra factorization:** structure = `dense`, backend = `dense_ldlt`.
- **Latent structure:** `4` random effects were estimated.
- **Symbolic vs numerical structure:** structural density = `1`, but 95% of curvature is retained by `9` entries.
- **Spectral complexity:** entropy effective rank = `3.9772256974195`, with 90% curvature requiring `4` eigen-directions.

## Model Health Assessment

| Check | Status | Evidence |
|---|---:|---|
| Latent mode | `PASS` | converged = `yes` |
| Gradient quality | `CHECK` | marginal fixed-effect gradient norm = `0.655796536381041` |
| Curvature | `PASS` | positive definite = `yes` |
| Conditioning | `EXCELLENT` | condition number = `1.31843011968793` |
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
| Entries for 95% curvature | `9` |
| Effective bandwidth for 95% curvature | `1` |
| 95% curvature compression | `1.77778x` |

## Parity-Point Evaluation

- Quality: `REVIEW`
- Objective value: `5.0511811613296`
- Marginal fixed-effect gradient norm: `0.655796536381041`
- Latent mode converged: `yes`
- Maximum marginal-gradient parameter: `fixed_0`

## Curvature

- Positive definite: `yes`
- Condition number: `1.31843011968793`
- Minimum eigenvalue: `5.01054878016581`
- Maximum eigenvalue: `6.60605842793621`

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

- Largest eigenvalue share: `0.291384759215163`
- Hessian inertia (positive / near-zero / negative): `4 / 0 / 0`
- Entropy effective rank: `3.9772256974195`
- Normalized spectral entropy: `0.995881216865985`
- Participation ratio: `3.95412929790155`
- Stable rank: `2.97862252080156`
- Leading spectral-gap ratio: `1.14077833377841`
- Eigenvectors needed for 90% curvature: `4`
- Eigenvectors needed for 95% curvature: `4`

**Interpretation:** curvature is distributed across many latent-state directions rather than being dominated by one or two modes. That is a good sign for numerical stability.

### Full Eigenvalue Spectrum

| Rank | Eigenvalue | Cumulative curvature share |
|---:|---:|---:|
| `rank_1` | `6.60605842793621` | `0.291384759215163` |
| `rank_2` | `5.79083440869362` | `0.546811033178534` |
| `rank_3` | `5.26381423656281` | `0.778991123713` |
| `rank_4` | `5.01054878016581` | `1` |

## Effective Structure

- Structural density: `1`
- Structural nonzeros: `16`
- Entries for 95% curvature: `9`
- Effective bandwidth for 95% curvature: `1`
- 95% curvature compression: `1.77778x`

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
- Mean: `0.417073989107842`
- Standard deviation: `0.282109175625996`

## Key Takeaway

This report demonstrates why Quadra's functional analysis diagnostics are useful: a model can look dense from a symbolic Hessian pattern, while numerical curvature, graph structure, and effective bandwidth reveal a simpler local-dependence structure.

## Full Laplace Structure Report

```text
Functional Analysis Report
==========================

Optimization
------------
objective_value:            5.0511811613296
gradient_norm:              0.655796536381041
max_gradient_parameter:     fixed_0
max_gradient_value:         -0.655796536381041
max_abs_gradient:           0.655796536381041
iterations:                 1
converged:                  yes
message:                    Converged: gradient norm below tolerance.

Curvature
---------
positive_definite:          yes
min_eigenvalue:             5.01054878016581
max_eigenvalue:             6.60605842793621
condition_number_abs:       1.31843011968793

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
largest_eigen_share:        0.291384759215163
effective_rank_entropy:     3.9772256974195
normalized_entropy:         0.995881216865985
participation_ratio:        3.95412929790155
stable_rank:                2.97862252080156
leading_gap_ratio:          1.14077833377841
eigen_count_for_50%:        2
eigen_count_for_90%:        4
eigen_count_for_95%:        4
eigen_count_for_99%:        4

rank,eigenvalue,cumulative_curvature_share
1,6.60605842793621,0.291384759215163
2,5.79083440869362,0.546811033178534
3,5.26381423656281,0.778991123713
4,5.01054878016581,1

Huu Structure
-------------
random_effects:             4
total_entries:              16
structural_nonzeros:        16
structural_density:         1

Effective Sparsity
------------------
curvature_retained,entries_required,entry_share,compression_vs_structural
90%,6,0.375,2.66666666666667
95%,9,0.5625,1.77777777777778
97%,10,0.625,1.6
98%,11,0.6875,1.45454545454545
99%,13,0.8125,1.23076923076923
99.5%,14,0.875,1.14285714285714
99.9%,16,1,1
100%,16,1,1

Effective Bandwidth
-------------------
curvature_retained,bandwidth,entry_count_if_banded,entry_share_if_banded
90%,0,4,0.25
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
min_variance:               0.176611672346553
max_variance:               0.181011065489628
min_variance_index:         0
max_variance_index:         3
max_abs_correlation:        0.0851972289113325
max_abs_correlation_pair:   0,1
count_abs_corr_gt_0_5:      0
count_abs_corr_gt_0_8:      0
count_abs_corr_gt_0_9:      0

Parameter Influence
-------------------
available:                  yes
Top parameter importance
index,name,variance,sd,variance_share,correlation_centrality,correlation_centrality_share,curvature_column_norm,curvature_diagonal,importance_score,importance_share
1,latent_1,0.177908472035967,0.421791977206735,0.249304866261558,0.184699276271753,0.32437641548594,5.75073548238375,5.70635120847091,2.87362325635419,0.261156172182867
2,latent_2,0.178086918421084,0.422003457830719,0.249554924910633,0.179780997300457,0.315738733004058,5.73636358518653,5.69509849714589,2.85597286249234,0.259552096461289
0,latent_0,0.176611672346553,0.420251915339541,0.247487648287407,0.106117624169834,0.186368107352239,5.73084800636272,5.70730614774164,2.66397352145339,0.242103110114003
3,latent_3,0.181011065489628,0.425453952255268,0.253652560540403,0.0988000839054795,0.173516744157762,5.58280089702708,5.5625,2.60989710657068,0.237188621241841

Top correlation pairs
i,j,name_i,name_j,correlation,abs_correlation
0,1,latent_0,latent_1,-0.0851972289113325,0.0851972289113325
1,2,latent_1,latent_2,-0.0832072736311268,0.0832072736311268
2,3,latent_2,latent_3,-0.0790793192935074,0.0790793192935074
0,2,latent_0,latent_2,-0.0174944043758226,0.0174944043758226
1,3,latent_1,latent_3,-0.0162947737292935,0.0162947737292935
0,3,latent_0,latent_3,-0.00342599088267868,0.00342599088267868

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
mean:                       0.417073989107842
sd:                         0.282109175625996
min_value:                  -0.0478086072981184
max_value:                  0.652897173713833
min_index:                  2
max_index:                  1
l2_norm:                    1.00704776324206
```
