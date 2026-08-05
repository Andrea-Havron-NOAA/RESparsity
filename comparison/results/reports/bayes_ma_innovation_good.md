# Bayesian MA Diagnostics: innovation_good

Quadra diagnostics at the parameterization parity evaluation point. The latent mode converged; the fixed effects were held fixed and were not optimized.

## Executive Summary

- **Overall status:** `REVIEW`.
- **Confidence:** `LOW`.
- **Parity-point quality:** `REVIEW`.
- **Uncertainty structure:** `LOCAL`.
- **Latent mode:** converged = `yes`, marginal fixed-effect gradient norm = `0.655796536381041`.
- **Curvature health:** positive definite = `yes`, condition number = `1.95600405252386`.
- **Quadra factorization:** structure = `tridiagonal`, backend = `tridiagonal`.
- **Latent structure:** `4` random effects were estimated.
- **Symbolic vs numerical structure:** structural density = `0.625`, but 95% of curvature is retained by `9` entries.
- **Spectral complexity:** entropy effective rank = `3.87882789690221`, with 90% curvature requiring `4` eigen-directions.

## Model Health Assessment

| Check | Status | Evidence |
|---|---:|---|
| Latent mode | `PASS` | converged = `yes` |
| Gradient quality | `CHECK` | marginal fixed-effect gradient norm = `0.655796536381041` |
| Curvature | `PASS` | positive definite = `yes` |
| Conditioning | `EXCELLENT` | condition number = `1.95600405252386` |
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
| Entries for 95% curvature | `9` |
| Effective bandwidth for 95% curvature | `1` |
| 95% curvature compression | `1.11111x` |

## Parity-Point Evaluation

- Quality: `REVIEW`
- Objective value: `5.0511811613296`
- Marginal fixed-effect gradient norm: `0.655796536381041`
- Latent mode converged: `yes`
- Maximum marginal-gradient parameter: `fixed_0`

## Curvature

- Positive definite: `yes`
- Condition number: `1.95600405252386`
- Minimum eigenvalue: `3.96075781839799`
- Maximum eigenvalue: `7.74725834385204`

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

- Largest eigenvalue share: `0.332952538259815`
- Hessian inertia (positive / near-zero / negative): `4 / 0 / 0`
- Entropy effective rank: `3.87882789690221`
- Normalized spectral entropy: `0.977810382547452`
- Participation ratio: `3.77064076961005`
- Stable rank: `2.39232501186755`
- Leading spectral-gap ratio: `1.19044874449018`
- Eigenvectors needed for 90% curvature: `4`
- Eigenvectors needed for 95% curvature: `4`

**Interpretation:** curvature is distributed across many latent-state directions rather than being dominated by one or two modes. That is a good sign for numerical stability.

### Full Eigenvalue Spectrum

| Rank | Eigenvalue | Cumulative curvature share |
|---:|---:|---:|
| `rank_1` | `7.74725834385204` | `0.332952538259815` |
| `rank_2` | `6.50784704483004` | `0.612639118468192` |
| `rank_3` | `5.05249325100038` | `0.829779218590983` |
| `rank_4` | `3.96075781839799` | `1` |

## Effective Structure

- Structural density: `0.625`
- Structural nonzeros: `10`
- Entries for 95% curvature: `9`
- Effective bandwidth for 95% curvature: `1`
- 95% curvature compression: `1.11111x`

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
- Mean: `0.518155648016863`
- Standard deviation: `0.233536717714028`

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
min_eigenvalue:             3.96075781839799
max_eigenvalue:             7.74725834385204
condition_number_abs:       1.95600405252386

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
largest_eigen_share:        0.332952538259815
effective_rank_entropy:     3.87882789690221
normalized_entropy:         0.977810382547452
participation_ratio:        3.77064076961005
stable_rank:                2.39232501186755
leading_gap_ratio:          1.19044874449018
eigen_count_for_50%:        2
eigen_count_for_90%:        4
eigen_count_for_95%:        4
eigen_count_for_99%:        4

rank,eigenvalue,cumulative_curvature_share
1,7.74725834385204,0.332952538259815
2,6.50784704483004,0.612639118468192
3,5.05249325100038,0.829779218590983
4,3.96075781839799,1

Huu Structure
-------------
random_effects:             4
total_entries:              16
structural_nonzeros:        10
structural_density:         0.625

Effective Sparsity
------------------
curvature_retained,entries_required,entry_share,compression_vs_structural
90%,8,0.5,1.25
95%,9,0.5625,1.11111111111111
97%,10,0.625,1
98%,10,0.625,1
99%,10,0.625,1
99.5%,10,0.625,1
99.9%,10,0.625,1
100%,10,0.625,1

Effective Bandwidth
-------------------
curvature_retained,bandwidth,entry_count_if_banded,entry_share_if_banded
90%,1,10,0.625
95%,1,10,0.625
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
min_variance:               0.176611672346553
max_variance:               0.187874118792824
min_variance_index:         0
max_variance_index:         3
max_abs_correlation:        0.207624153507196
max_abs_correlation_pair:   2,3
count_abs_corr_gt_0_5:      0
count_abs_corr_gt_0_8:      0
count_abs_corr_gt_0_9:      0

Parameter Influence
-------------------
available:                  yes
Top parameter importance
index,name,variance,sd,variance_share,correlation_centrality,correlation_centrality_share,curvature_column_norm,curvature_diagonal,importance_score,importance_share
2,latent_2,0.184554311587924,0.42959784867702,0.251732162843913,0.454785543908904,0.321322480721309,6.1276958503365,5.90195215269348,3.82964286516317,0.27536437577046
1,latent_1,0.184097488842026,0.429065832759993,0.251109056415937,0.449981114355111,0.317927976996748,6.1276958503365,5.90195215269348,3.81226848438784,0.274115098570198
0,latent_0,0.176611672346553,0.420251915339541,0.240898399396632,0.25164802480372,0.177798456176843,6.01588296291219,5.90195215269348,3.16439943582078,0.227531105644249
3,latent_3,0.187874118792824,0.433444481788411,0.256260381343518,0.258940827969013,0.182951086105101,5.68323982080415,5.5625,3.10123573166752,0.222989420015093

Top correlation pairs
i,j,name_i,name_j,correlation,abs_correlation
2,3,latent_2,latent_3,0.207624153507196,0.207624153507196
1,2,latent_1,latent_2,0.205697754956321,0.205697754956321
0,1,latent_0,latent_1,0.201575537147654,0.201575537147654
1,3,latent_1,latent_3,0.0427078222511369,0.0427078222511369
0,2,latent_0,latent_2,0.0414636354453868,0.0414636354453868
0,3,latent_0,latent_3,0.00860885221067942,0.00860885221067942

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
mean:                       0.518155648016863
sd:                         0.233536717714028
min_value:                  0.17893556839791
max_value:                  0.778353445763383
min_index:                  2
max_index:                  1
l2_norm:                    1.13670519325358
```
