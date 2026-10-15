# Count-model sensitivity benchmark

Date: 2026-10-15. This retrospective analysis follows the frozen plan in RETROSPECTIVE_ANALYSIS_PLAN.md. Existing manuscript outcomes were already known. No findings were withheld.

## Main finding

Count-aware inference supports the original aggregate condition-associated structure, while narrowing local candidate claims. Exact DM-score tests identify 20 positive contrast–cluster tests (19 distinct clusters): 18 flight–ground and 2 flight–vivarium. Sixteen overlap the original 43 positives. All three full-depth global comparisons and all 30 original A-half global comparisons survive their respective Holm corrections. However, the count-model sensitivity does not establish local read-depth stability, and it does not test independent-mission transfer. Each of the ten original A-half DM families has zero local BY discoveries. Complementary B-half draws were not included in this frozen benchmark and cannot be assumed to have the same result. The parent complementary-split audit retains additional B-half outcomes; no universal failure at 50% depth is claimed. The stricter 95% gate also has zero DM-score BY discoveries, with a discrete correction boundary explained below.

## Statistical model and calibration

For animal i and local alternative j, x_ij are deposited STAR counts and N_i their selected sum. A Dirichlet–multinomial model has alpha_j > 0, A=sum_j alpha_j, pi_j=alpha_j/A and rho=1/(A+1). Its log likelihood, omitting the fixed multinomial coefficient, is sum_i{log Gamma(A)-log Gamma(N_i+A)+sum_j[log Gamma(x_ij+alpha_j)-log Gamma(alpha_j)]}. Two fixed starting concentrations (10 and 1,000) initialize bounded optimization in log alpha. Fits use all 26 animals without group labels. The higher likelihood among converged solutions is selected; complete diagnostic records are retained. A pooled model intentionally absorbs both within-arm and between-arm heterogeneity, so rho must not be interpreted as a within-group dispersion or biological effect.

The mean-logit score for each animal is S_ij=alpha_j{d_ij-sum_l pi_l d_il}, where d_ij=digamma(x_ij+alpha_j)-digamma(alpha_j). For each comparison, its pooled empirical covariance is eigen-whitened with the prespecified relative eigenvalue tolerance 10^-10. The test is the squared distance between group-average whitened score vectors. The multinomial-limit sensitivity uses S_ij=x_ij-N_i pi_j and the identical calibration. Exact tests enumerate 48,620 flight–ground assignments or 24,310 assignments for each vivarium comparison. Inclusive ties use the stated numerical tolerance. No asymptotic likelihood-ratio or chi-square p values are used. All 923 clusters × 3 contrasts remain in each BY/BH correction family.

This is an independently specified, permutation-calibrated DM-score sensitivity analysis using an established count likelihood. It is not the DRIMSeq package, its dispersion shrinkage or regression likelihood-ratio workflow; it is not LeafCutter/rMATS. Standard likelihood scores provide a scientifically interpretable check, while the empirical whitening/permutation wrapper and lack of a complete independent software benchmark are explicit limits. Whitening uses at most K-1 dimensions and is noisy with small cohorts, despite full retained numerical rank in this dataset. Permutation validity is conditional on exchangeability of complete animal count vectors; depth differences, unmeasured batch, and collection date remain relevant. Exact arithmetic does not remove these assumptions.

## Complete model agreement

Native animal-average composition effects are kept identical across model comparisons. Their units do not change when a different test statistic is used. Evidence-rank concordance and overlap, rather than the union of favorable p values, quantify model agreement.

| model | contrast | n_tests | n_valid | n_original_BY | n_model_BY | n_shared_BY | spearman_negative_log_p | median_covariance_condition |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DM_score | FLT_vs_GC | 923 | 923 | 40 | 18 | 15 | 0.943475 | 5.9761493 |
| DM_score | FLT_vs_VIV | 923 | 923 | 3 | 2 | 1 | 0.9207111 | 5.5481392 |
| DM_score | GC_vs_VIV | 923 | 923 | 0 | 0 | 0 | 0.9208067 | 5.7255117 |
| multinomial_score | FLT_vs_GC | 923 | 923 | 40 | 21 | 19 | 0.8944183 | 7.1647746 |
| multinomial_score | FLT_vs_VIV | 923 | 923 | 3 | 2 | 0 | 0.8503727 | 6.3938388 |
| multinomial_score | GC_vs_VIV | 923 | 923 | 0 | 0 | 0 | 0.8429891 | 6.3788732 |

The DM score yields 20 BY-positive tests versus 23 for the multinomial-limit score. Their differences show sensitivity to both the score transformation and covariance geometry. The original method has 43 positive tests; 16 are shared with DM. All 20 DM findings and their original annotation are retained in all20_DM_BY_tests_with_original_annotation.csv, without substituting candidates into the manuscript based on this analysis. Complete results include all nonsignificant tests, all covariance ranks, condition numbers, raw exact fractions and both BH/BY values.

## Fitting and numerical diagnostics

All 10,153 fits (923 full depth plus 9,230 original A-half fits) report convergence. Full-depth median rho is 0.017242 (interquartile range 0.010241–0.027525). One full-depth fit and 44 original A-half fits touch an alpha bound. At full depth the two starting points agree in likelihood within 2.85e-07; three thinning fits have starting-solution likelihood gaps greater than 0.01 (maximum 2.749300), and the prespecified better converged solution is retained. This is an informative numerical limitation, not proof of a unique global maximum. Maximum gradients and all messages are disclosed. The exact label calibration depends on a fixed label-blind transformation and does not require that a particular optimizer recover a biological true parameter.

| replicate | n_fits | n_converged | median_rho | n_bound_hits | median_iterations | maximum_gradient | maximum_start_loglik_gap |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 923 | 923 | 0.0163074 | 5 | 13.0 | 0.00165271 | 0.00010757 |
| 10 | 923 | 923 | 0.01602475 | 4 | 13.0 | 0.01269258 | 2.74929993 |
| 2 | 923 | 923 | 0.01553157 | 4 | 13.0 | 0.00574008 | 0.00028275 |
| 3 | 923 | 923 | 0.01576959 | 5 | 13.0 | 0.00687347 | 0.0004373 |
| 4 | 923 | 923 | 0.01590182 | 9 | 13.0 | 0.02440307 | 0.57999384 |
| 5 | 923 | 923 | 0.01618111 | 3 | 13.0 | 0.04293304 | 0.00087183 |
| 6 | 923 | 923 | 0.0159011 | 2 | 13.0 | 0.0134538 | 0.00048454 |
| 7 | 923 | 923 | 0.01593372 | 6 | 13.0 | 0.00327533 | 0.30418979 |
| 8 | 923 | 923 | 0.01619643 | 3 | 13.0 | 0.01476863 | 0.00118165 |
| 9 | 923 | 923 | 0.01619839 | 3 | 13.0 | 0.00769658 | 0.00053904 |
| full | 923 | 923 | 0.01724208 | 1 | 12.0 | 0.0118951 | 2.9e-07 |

## Global association survives count-model sensitivity

Global scores average the whitened local statistic divided by its covariance rank, equally across valid frozen clusters. Because they have a different normalization and geometry from Hellinger distance, their numeric magnitudes are not directly interchangeable with previous effect sizes.

| replicate | model | contrast | n_valid_clusters | statistic | n_ge_observed | n_assignments | exact_p | Holm_p_3 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| full | DM_score | FLT_vs_GC | 923 | 0.463372073 | 2 | 48620 | 4.1135e-05 | 0.000123406 |
| full | DM_score | FLT_vs_VIV | 923 | 0.311208378 | 48 | 24310 | 0.001974496 | 0.00370218 |
| full | DM_score | GC_vs_VIV | 923 | 0.319802324 | 45 | 24310 | 0.00185109 | 0.00370218 |
| full | multinomial_score | FLT_vs_GC | 923 | 0.454490313 | 2 | 48620 | 4.1135e-05 | 0.000123406 |
| full | multinomial_score | FLT_vs_VIV | 923 | 0.313238004 | 36 | 24310 | 0.001480872 | 0.002961744 |
| full | multinomial_score | GC_vs_VIV | 923 | 0.319302233 | 47 | 24310 | 0.001933361 | 0.002961744 |

The ten original A-half read-thinnings use the original all-animal count and dominance gates (580–600 valid clusters). All 30 global model-score tests survive Holm within replicate; there is no claim of a simultaneous error rate across this retrospective family of sensitivity analyses.

| contrast | n_thinnings | n_Holm_positive | minimum_Holm | maximum_Holm |
| --- | --- | --- | --- | --- |
| FLT_vs_GC | 10 | 10 | 0.000123406 | 0.000123406 |
| FLT_vs_VIV | 10 | 10 | 0.001151789 | 0.005841218 |
| GC_vs_VIV | 10 | 10 | 0.001151789 | 0.005841218 |

## Discrete BY correction explains the stricter-gate boundary

The smallest attainable two-sided exact p is 2/48,620=1/24,310=4.113533525e-05. With m=2,769, BY uses H_m=sum_{i=1}^m 1/i. Thus m H_m p_min=0.968596171, and at least ceil(m H_m p_min/0.05)=20 tests at this floor are needed for that tied minimum-p block to cross BY 0.05 in isolation. All 20 full-depth DM discoveries lie at this floor and share q=0.048429809. The 95% dominance gate removes one of the floor tests, leaving 19 and a minimum q=0.050978746; every retained raw p stays unchanged. No local test passes. This is a correction-family resolution boundary rather than evidence that the measured Tcaf1 effect vanished. It is not a rationale to loosen BY or shrink the denominator. Other, higher-p ranks can in principle create their own step-up rejection block; the complete ranked family was checked here.

| family | n_tests | smallest_attainable_two_sided_p | n_at_attainable_floor | minimum_observed_BY_q | n_BY | BY_floor_multiplier | minimum_floor_ties_for_BY05 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| original_hellinger | 2769 | 4.1135e-05 | 30 | 0.032286539 | 43 | 0.968596171 | 20 |
| DM_full | 2769 | 4.1135e-05 | 20 | 0.048429809 | 20 | 0.968596171 | 20 |
| multinomial_full | 2769 | 4.1135e-05 | 23 | 0.042112877 | 23 | 0.968596171 | 20 |
| DM_dominance95 | 2769 | 4.1135e-05 | 19 | 0.050978746 | 0 | 0.968596171 | 20 |
| DM_thinning_1 | 2769 | 4.1135e-05 | 9 | 0.107621797 | 0 | 0.968596171 | 20 |
| DM_thinning_2 | 2769 | 4.1135e-05 | 7 | 0.138370882 | 0 | 0.968596171 | 20 |
| DM_thinning_3 | 2769 | 4.1135e-05 | 8 | 0.121074521 | 0 | 0.968596171 | 20 |
| DM_thinning_4 | 2769 | 4.1135e-05 | 8 | 0.121074521 | 0 | 0.968596171 | 20 |
| DM_thinning_5 | 2769 | 4.1135e-05 | 11 | 0.088054197 | 0 | 0.968596171 | 20 |
| DM_thinning_6 | 2769 | 4.1135e-05 | 6 | 0.161432695 | 0 | 0.968596171 | 20 |
| DM_thinning_7 | 2769 | 4.1135e-05 | 11 | 0.088054197 | 0 | 0.968596171 | 20 |
| DM_thinning_8 | 2769 | 4.1135e-05 | 13 | 0.074507398 | 0 | 0.968596171 | 20 |
| DM_thinning_9 | 2769 | 4.1135e-05 | 8 | 0.121074521 | 0 | 0.968596171 | 20 |
| DM_thinning_10 | 2769 | 4.1135e-05 | 9 | 0.107621797 | 0 | 0.968596171 | 20 |

## Tcaf1 remains specific to the flight–ground comparison

The previously discussed cluster 6:2:42663893:start has pooled rho=0.015649. Its native flight–ground maximum alternative difference remains 16.41 percentage points and flight–vivarium 12.27 percentage points, by definition using unchanged animal means. With DM scores, flight–ground p=2/48,620 and BY q=0.048429809; flight–vivarium p=88/24,310=0.00361991 and q=0.546388. The multinomial score similarly retains flight–ground (q=0.042113) but not flight–vivarium (q=0.910893). Under the stricter gate, flight–ground raw p remains identical while q rises to 0.050979. All original A-half Tcaf1 results are disclosed. These results narrow its interpretation to a comparator- and method-sensitive short-read junction candidate, not a validated functional isoform.

| model | contrast | statistic | n_ge_observed | n_assignments | exact_p | BY_q_full_2769 | native_mean_total_variation |
| --- | --- | --- | --- | --- | --- | --- | --- |
| DM_score | FLT_vs_GC | 3.080825521 | 2 | 48620 | 4.1135e-05 | 0.048429809 | 0.16408843 |
| DM_score | FLT_vs_VIV | 1.949105667 | 88 | 24310 | 0.00361991 | 0.546387584 | 0.122661485 |
| DM_score | GC_vs_VIV | 1.007244594 | 2815 | 24310 | 0.115795969 | 1.0 | 0.041426944 |
| multinomial_score | FLT_vs_GC | 3.043417611 | 2 | 48620 | 4.1135e-05 | 0.042112877 | 0.16408843 |
| multinomial_score | FLT_vs_VIV | 1.702014442 | 221 | 24310 | 0.009090909 | 0.91089257 | 0.122661485 |
| multinomial_score | GC_vs_VIV | 1.201962853 | 1767 | 24310 | 0.072686137 | 1.0 | 0.041426944 |

## Independent checks and reproducibility

190/190 checks pass. The verifier does not import analysis code. It compares the scalar likelihood to scipy.stats.dirichlet_multinomial.logpmf, finite-difference mean-logit derivatives to stored scores, direct covariance-pseudoinverse statistics to whitened statistics, and every assignment for six fixed cluster indices (including Tcaf1) across all three contrasts and both models. Exhaustive reference labels show super-uniform empirical rejection frequencies at 0.01, 0.05 and 0.10 with inclusive ties. It also checks every input hash, all count-to-composition mappings, all model ranks, retained native effects, every BY family, failed-gate p=1 behavior and every saved full-depth global null. These checks validate implementation and conditional reference arithmetic, not causal identification or independent replication.

Run scripts/analyse_benchmark.py, scripts/verify_independently.py, then scripts/build_report_figure.py from any working directory using the package Python dependencies. The extension requires the parent package frozen counts and saved thinning matrices but no live data downloads. Figure outputs are vector PDF/SVG and RGB 600-dpi PNG at 183×160 mm with minimum 8-pt Arial/Helvetica text.

## Implications for the manuscript

The defensible stronger result is that broad comparator-associated cerebellar junction structure is concordant across a square-root composition test and an explicitly count-based exact score analysis. Local discoveries depend on the model, read depth and the discrete correction family. The full 20-test count-model result is a sensitivity result and does not replace the primary 43-test result or expand a pooled discovery list. Independent-mission nontransfer and the perfect RR23 condition–collection-date alias remain unchanged. Neither additional likelihood assumptions nor repeated computations create causal identification, independent biology, protein evidence, or a guaranteed npj publication. The study gains credibility by precisely locating which claims survive these checks.

## Verified references

1. Mosimann, J. E. On the compound multinomial distribution, the multivariate β-distribution, and correlations among proportions. Biometrika 49, 65–82 (1962). https://doi.org/10.1093/biomet/49.1-2.65
2. Nowicka, M. & Robinson, M. D. DRIMSeq: a Dirichlet-multinomial framework for multivariate count outcomes in genomics. F1000Research 5, 1356 (2016). https://doi.org/10.12688/f1000research.8900.2
3. Winkler, A. M. et al. Permutation inference for the general linear model. NeuroImage 92, 381–397 (2014). https://doi.org/10.1016/j.neuroimage.2014.01.060

Crossref metadata was retrieved on 2026-10-15 and retained in inputs/verified_reference_metadata.json. The mathematical distribution definition was independently checked against SciPy documentation (saved HTML and URL) and its evaluated log PMF. Incorrect initial DOI lookup candidates are kept in an explicitly named audit file and are not cited as sources.
