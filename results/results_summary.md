# Reduced statistical set: re-run results

N layouts = 300; N errors = 2708; categories = 13. Seed = 2026.

## RQ1. Error shares and layout prevalence (Wilson 95% CI)

| Domain | Errors | Share % | Layouts affected | Prevalence % [95% CI] |
|---|---|---|---|---|
| Graphic Drafting | 1080 | 39.88 | 269 | 89.7 [85.7, 92.6] |
| Spatial Circulation | 879 | 32.46 | 243 | 81.0 [76.2, 85.0] |
| Ergonomic Layout | 672 | 24.82 | 251 | 83.7 [79.1, 87.4] |
| Structural Logic | 48 | 1.77 | 44 | 14.7 [11.1, 19.1] |
| Regulatory Compliance | 29 | 1.07 | 23 | 7.7 [5.2, 11.2] |

Layouts with at least one error: 300/300. Category-level table: rq1_shares_prevalence.csv

## RQ2. Category distribution by model (layout-level permutation, 10,000 permutations)

Omnibus: chi2 = 274.0, Cramer V = 0.184, permutation p = 0.0001 (min attainable 0.0001); permutation-null median chi2 = 55.3 versus naive df = 36, i.e. error-level chi2 is inflated by clustering.

Planned contrasts (Bonferroni, M = 2):

| Contrast | chi2 | V | perm p | Bonferroni p |
|---|---|---|---|---|
| GPT-Image-1 vs GPT-Image-2 | 109.2 | 0.274 | 0.0001 | 0.0002 |
| Nano Banana Pro vs Nano Banana 2 | 7.8 | 0.079 | 0.9508 | 1.0000 |

Errors per layout by model:

| model           |   count |   mean |   std |
|:----------------|--------:|-------:|------:|
| GPT-Image-1     |      79 |   9.91 |  4.87 |
| GPT-Image-2     |      79 |   8.53 |  4.89 |
| Nano Banana 2   |      66 |   8.77 |  5.54 |
| Nano Banana Pro |      76 |   8.84 |  6.14 |

## RQ3. Negative binomial GLM with log(rooms) as a covariate (reference: Rectangular, Nano Banana Pro; robust HC1 SE)

| Term | IRR | 95% CI | p |
|---|---|---|---|
| Asymmetrical/Irregular | 1.115 | [0.942, 1.319] | 0.205 |
| Composite/Sprawling | 1.258 | [0.932, 1.698] | 0.134 |
| L-Shape | 0.899 | [0.738, 1.095] | 0.290 |
| GPT-Image-1 | 1.152 | [0.952, 1.394] | 0.145 |
| GPT-Image-2 | 0.983 | [0.806, 1.199] | 0.869 |
| Nano Banana 2 | 0.994 | [0.808, 1.224] | 0.958 |
| bedroom_count | 1.031 | [0.930, 1.144] | 0.562 |
| floor_area | 0.999 | [0.997, 1.001] | 0.349 |
| log_rooms | 1.164 | [0.869, 1.560] | 0.307 |

log(rooms) coefficient (elasticity) = 0.152 (95% CI -0.140 to 0.445); the offset model assumes exactly 1.0.
MLE alpha = 0.217. AIC: covariate model 1802.9 vs offset model 1834.3.
Spearman rho (prompted rooms vs errors per layout) = 0.153, p = 0.0078.
Predictor collinearity (Spearman): rooms-bedrooms 0.66, rooms-area 0.76, bedrooms-area 0.68.

## RQ4a. Pairwise co-occurrence: logistic regression of A on B + log(rooms), all 78 pairs, BH-FDR

Nominal p < 0.05: 11 pairs (about 3.9 expected by chance). BH q < 0.05: 0. Bonferroni < 0.05: 0. Smallest q = 0.067.

Top 8 by p:

| A | B | both/nA/nB | phi | OR [95% CI] | p | q (BH) |
|---|---|---|---|---|---|---|
| Dimension & Measurement Errors | Scale Errors | 32/97/66 | 0.183 | 2.48 [1.41, 4.37] | 0.0016 | 0.067 |
| Line & Layer Errors | Missing Openings & Elements | 84/117/182 | 0.182 | 2.21 [1.35, 3.64] | 0.0017 | 0.067 |
| Line & Layer Errors | Unclear Wall Types | 16/117/23 | 0.181 | 3.96 [1.58, 9.96] | 0.0034 | 0.089 |
| Redundant / Unnecessary Elements | Space & Clearance Errors | 77/199/99 | 0.170 | 2.15 [1.23, 3.75] | 0.0070 | 0.119 |
| Missing Openings & Elements | Symbol Errors | 106/182/156 | 0.155 | 1.90 [1.19, 3.04] | 0.0076 | 0.119 |
| Dimension & Measurement Errors | Missing Openings & Elements | 69/97/182 | 0.148 | 1.98 [1.18, 3.34] | 0.0100 | 0.131 |
| Missing Required Spaces | Space & Clearance Errors | 5/35/99 | -0.145 | 0.30 [0.11, 0.80] | 0.0165 | 0.167 |
| Poor Room Layout & Placement | Space & Clearance Errors | 70/182/99 | 0.144 | 1.88 [1.12, 3.15] | 0.0171 | 0.167 |

## RQ4b. Co-occurrence network (phi > 0.05; 34 edges) and bootstrap stability (2000 resamples of layouts)

| Category | Degree | Betweenness (dist = 1/phi) | P(top betweenness) | P(top degree) | Eigenvector (suppl.) | P(top eigenvector) |
|---|---|---|---|---|---|---|
| Symbol Errors | 0.750 | 0.288 | 0.40 | 0.54 | 0.381 | 0.25 |
| Missing Openings & Elements | 0.583 | 0.197 | 0.19 | 0.13 | 0.428 | 0.30 |
| Line & Layer Errors | 0.583 | 0.121 | 0.22 | 0.12 | 0.415 | 0.24 |
| Redundant / Unnecessary Elements | 0.583 | 0.106 | 0.06 | 0.10 | 0.333 | 0.08 |
| Dimension & Measurement Errors | 0.417 | 0.091 | 0.01 | 0.01 | 0.333 | 0.04 |
| Unclear Wall Types | 0.417 | 0.045 | 0.03 | 0.01 | 0.252 | 0.01 |

Full table (13 categories, all variants): rq4_network_centrality_bootstrap.csv

Betweenness of Symbol Errors under alternative definitions: distance = 1/phi (used in the paper) 0.288; unweighted 0.250; phi as length 0.303.
