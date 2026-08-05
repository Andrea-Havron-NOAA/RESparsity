# Bayesian AR Diagnostics: centered

Quadra diagnostics at the parameterization parity evaluation point. The latent mode converged; the fixed effects were held fixed and were not optimized.

## Executive Summary

- **Overall status:** `REVIEW`.
- **Confidence:** `LOW`.
- **Parity-point quality:** `REVIEW`.
- **Uncertainty structure:** `LOCAL`.
- **Latent mode:** converged = `yes`, marginal fixed-effect gradient norm = `0.47236761287214`.
- **Curvature health:** positive definite = `yes`, condition number = `1.20162758653246`.
- **Quadra factorization:** structure = `tridiagonal`, backend = `tridiagonal`.
- **Latent structure:** `4` random effects were estimated.
- **Symbolic vs numerical structure:** structural density = `0.625`, but 95% of curvature is retained by `7` entries.
- **Spectral complexity:** entropy effective rank = `3.99036361958929`, with 90% curvature requiring `4` eigen-directions.

## Model Health Assessment

| Check | Status | Evidence |
|---|---:|---|
| Latent mode | `PASS` | converged = `yes` |
| Gradient quality | `CHECK` | marginal fixed-effect gradient norm = `0.47236761287214` |
| Curvature | `PASS` | positive definite = `yes` |
| Conditioning | `EXCELLENT` | condition number = `1.20162758653246` |
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

## Parity-Point Evaluation

- Quality: `REVIEW`
- Objective value: `4.17136656947898`
- Marginal fixed-effect gradient norm: `0.47236761287214`
- Latent mode converged: `yes`
- Maximum marginal-gradient parameter: `fixed_0`

## Curvature

- Positive definite: `yes`
- Condition number: `1.20162758653246`
- Minimum eigenvalue: `5.16012260443924`
- Maximum eigenvalue: `6.2005456713839`

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

- Largest eigenvalue share: `0.273231660392753`
- Hessian inertia (positive / near-zero / negative): `4 / 0 / 0`
- Entropy effective rank: `3.99036361958929`
- Normalized spectral entropy: `0.998260108596734`
- Participation ratio: `3.98080297243514`
- Stable rank: `3.36486213075471`
- Leading spectral-gap ratio: `1.05709609323106`
- Eigenvectors needed for 90% curvature: `4`
- Eigenvectors needed for 95% curvature: `4`

**Interpretation:** curvature is distributed across many latent-state directions rather than being dominated by one or two modes. That is a good sign for numerical stability.

### Full Eigenvalue Spectrum

| Rank | Eigenvalue | Cumulative curvature share |
|---:|---:|---:|
| `rank_1` | `6.2005456713839` | `0.273231660392753` |
| `rank_2` | `5.8656405137509` | `0.531705475727372` |
| `rank_3` | `5.46705574391751` | `0.772615356492257` |
| `rank_4` | `5.16012260443924` | `1` |

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
- Mean: `0.38886797426605`
- Standard deviation: `0.289135047153467`

## Key Takeaway

This report demonstrates why Quadra's functional analysis diagnostics are useful: a model can look dense from a symbolic Hessian pattern, while numerical curvature, graph structure, and effective bandwidth reveal a simpler local-dependence structure.

## Full Laplace Structure Report

```text
Functional Analysis Report
==========================

Parity-Point and Latent-Mode Status
-----------------------------------
objective_value:            4.17136656947898
marginal_fixed_gradient_norm: 0.47236761287214
max_marginal_gradient_parameter: fixed_0
max_marginal_gradient_value:  -0.47236761287214
max_abs_marginal_gradient:    0.47236761287214
latent_mode_iterations:       1
latent_mode_converged:        yes
latent_mode_message:          Converged: latent gradient norm below tolerance.

Curvature
---------
positive_definite:          yes
min_eigenvalue:             5.16012260443924
max_eigenvalue:             6.2005456713839
condition_number_abs:       1.20162758653246

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
largest_eigen_share:        0.273231660392753
effective_rank_entropy:     3.99036361958929
normalized_entropy:         0.998260108596734
participation_ratio:        3.98080297243514
stable_rank:                3.36486213075471
leading_gap_ratio:          1.05709609323106
eigen_count_for_50%:        2
eigen_count_for_90%:        4
eigen_count_for_95%:        4
eigen_count_for_99%:        4

rank,eigenvalue,cumulative_curvature_share
1,6.2005456713839,0.273231660392753
2,5.8656405137509,0.531705475727372
3,5.46705574391751,0.772615356492257
4,5.16012260443924,1

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
min_variance:               0.176335211134886
max_variance:               0.17832688748533
min_variance_index:         0
max_variance_index:         3
max_abs_correlation:        0.0568127710393671
max_abs_correlation_pair:   2,3
count_abs_corr_gt_0_5:      0
count_abs_corr_gt_0_8:      0
count_abs_corr_gt_0_9:      0

Parameter Influence
-------------------
available:                  yes
Top parameter importance
index,name,variance,sd,variance_share,correlation_centrality,correlation_centrality_share,curvature_column_norm,curvature_diagonal,importance_score,importance_share
2,latent_2,0.17690620387186,0.420602191948473,0.249702409347893,0.116596111391438,0.330322666851352,5.7072473871507,5.68917558099758,2.68036748323358,0.257259995053755
1,latent_1,0.176899846989772,0.420594634998798,0.249693436633885,0.116296925549559,0.329475058264803,5.7072473871507,5.68917558099758,2.67960114818044,0.257186442694531
0,latent_0,0.176335211134886,0.419922863315259,0.248896455350677,0.0598741517691259,0.169626493129277,5.69821864838292,5.68917558099758,2.53607989686929,0.24341136273502
3,latent_3,0.17832688748533,0.422287683321844,0.251707698667545,0.0602092282667579,0.170575781754568,5.63498250382057,5.62583779049879,2.52285660569801,0.242142199516693

Top correlation pairs
i,j,name_i,name_j,correlation,abs_correlation
2,3,latent_2,latent_3,0.0568127710393671,0.0568127710393671
1,2,latent_1,latent_2,0.0565864525095903,0.0565864525095903
0,1,latent_0,latent_1,0.0564956398696118,0.0564956398696118
1,3,latent_1,latent_3,0.00321483317035737,0.00321483317035737
0,2,latent_0,latent_2,0.0031968878424807,0.0031968878424807
0,3,latent_0,latent_3,0.000181624057033392,0.000181624057033392

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
mean:                       0.38886797426605
sd:                         0.289135047153467
min_value:                  -0.0555808715426449
max_value:                  0.751569092823876
min_index:                  0
max_index:                  2
l2_norm:                    0.969159175578953
```
