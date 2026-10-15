# RR23 multiway robustness: additional retrospective analysis

Author: Isarar Siddique

Affiliation:

Email:

Date: 15 October 2026. Plan frozen before new outcomes. This family reuses the original 26 animals, 923 clusters and 2,769 local tests; it is not independent replication or prospective preregistration.

## Main conclusion

Aggregate global multiway associations withstand the stronger gate and all ten 50% depth thinnings, but the local findings are depth sensitive. The 95% gate retains 25 of the original 43 significant contrast–cluster tests (24 distinct clusters): 23/40 flight–ground and 2/3 flight–vivarium. Every thinning has zero BY-positive local tests when the full family remains 2,769. None of the original 41 unique local BY clusters is recovered in any thinning at this threshold. Tcaf1 remains measurable in every thinning but loses local BY significance in both contrasts; its original flight–vivarium result also fails the 95% gate sensitivity after full-family correction. Original local findings must therefore be described as retrospective, depth-sensitive associations. Global stability cannot be used to certify individual junction choices, and no result is replaced by a favorable new candidate.

## Locked methods and provenance

The plan SHA256 is 56d56dcd96b34d12421b15689a56bbf7cdad833ca1f7b936dfac4041eb2a7a0a. Inputs are linked exclusively to the NEW npj package in frozen_input_links.csv with file hashes. The original full count matrix has 805,562 junctions×26 collapsed animals. Structural definitions and all alternative lists are inherited unchanged. For the 95% gate, selected/all canonical boundary counts must reach 95% in every animal; failed clusters receive p=1 for every comparison. Passing clusters retain their original exact p values, but full-family BH/BY ranks change when other measurements are lost.

For each of ten replicates, every cell of the ENTIRE count matrix is independently Binomial(original count,0.5), with PCG64 seeds 2026101501–2026101510. Each junction is sampled once and reused at all overlapping boundaries. All ten complete matrices are retained in original row order; source hash, seed and file hash are recorded. Alternatives are not reselected and pooled composition gates are not tuned. A thinned measurement must retain ≥20 reads and ≥90% selected/all canonical boundary dominance in every animal. Failed measurements are p=1 in all 2,769 slots. Global secondary comparisons equally weight the valid frozen-cluster subset and report its size. Exact sqrt-centroid tests enumerate every 48,620/24,310/24,310 assignment, including observed/complement ties, with 10⁻¹² absolute tolerance. Holm3 is within each explicitly secondary global family/replicate; local BH/BY is within all 2,769 tests in each sensitivity. This is not a combined correction over all retrospective analyses ever performed.

## Local stability, including null results

| contrast | original_positive_tests | survive_measurement95 | BY_positive95 | total_valid_thinning_tests | BY_positive_in_any_thinning |
| --- | --- | --- | --- | --- | --- |
| FLT_vs_GC | 40 | 33 | 23 | 345 | 0 |
| FLT_vs_VIV | 3 | 3 | 2 | 29 | 0 |

The gate retains 607/923 measurements; it does not discard failed hypotheses from the correction family. Of the original positive tests, 33/40 flight–ground and all 3 flight–vivarium remain measurable at 95%, but only 23 and 2 survive BY. Across the ten thinnings, 345/400 original flight–ground test instances and 29/30 original flight–vivarium instances remain measurable (374/430 total), while zero survive BY. Thus loss of significance is not explained only by measurement exclusion: many passing original events fail the fixed-family correction after binomial read noise and changes in the rank distribution. Point estimates, raw p values and both corrections for every event/replicate are supplied; no threshold or thinning fraction was changed after observing these losses.

| replicate | contrast | n_tests | n_measurement_valid | min_exact_p | min_BH_q | min_BY_q | BH_positive | BY_positive |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | FLT_vs_GC | 923 | 587 | 4.1135e-05 | 0.006327986 | 0.053810898 | 68 | 0 |
| 1 | FLT_vs_VIV | 923 | 587 | 0.000123406 | 0.011390374 | 0.096859617 | 7 | 0 |
| 1 | GC_vs_VIV | 923 | 587 | 8.2271e-05 | 0.008761826 | 0.074507398 | 7 | 0 |
| 2 | FLT_vs_GC | 923 | 591 | 4.1135e-05 | 0.008761826 | 0.074507398 | 78 | 0 |
| 2 | FLT_vs_VIV | 923 | 591 | 8.2271e-05 | 0.011390374 | 0.096859617 | 10 | 0 |
| 2 | GC_vs_VIV | 923 | 591 | 0.00061703 | 0.02669619 | 0.227014728 | 5 | 0 |
| 3 | FLT_vs_GC | 923 | 599 | 4.1135e-05 | 0.007593583 | 0.064573078 | 82 | 0 |
| 3 | FLT_vs_VIV | 923 | 599 | 4.1135e-05 | 0.007593583 | 0.064573078 | 16 | 0 |
| 3 | GC_vs_VIV | 923 | 599 | 0.000452489 | 0.020540019 | 0.174664883 | 10 | 0 |
| 4 | FLT_vs_GC | 923 | 598 | 4.1135e-05 | 0.008761826 | 0.074507398 | 83 | 0 |
| 4 | FLT_vs_VIV | 923 | 598 | 8.2271e-05 | 0.010354886 | 0.088054197 | 13 | 0 |
| 4 | GC_vs_VIV | 923 | 598 | 4.1135e-05 | 0.008761826 | 0.074507398 | 7 | 0 |
| 5 | FLT_vs_GC | 923 | 590 | 4.1135e-05 | 0.00670022 | 0.056976245 | 77 | 0 |
| 5 | FLT_vs_VIV | 923 | 590 | 8.2271e-05 | 0.008135982 | 0.069185441 | 14 | 0 |
| 5 | GC_vs_VIV | 923 | 590 | 0.000411353 | 0.021093286 | 0.179369661 | 7 | 0 |
| 6 | FLT_vs_GC | 923 | 584 | 4.1135e-05 | 0.005994934 | 0.050978746 | 79 | 0 |
| 6 | FLT_vs_VIV | 923 | 584 | 4.1135e-05 | 0.005994934 | 0.050978746 | 13 | 0 |
| 6 | GC_vs_VIV | 923 | 584 | 0.00061703 | 0.025887214 | 0.220135493 | 13 | 0 |
| 7 | FLT_vs_GC | 923 | 600 | 4.1135e-05 | 0.007593583 | 0.064573078 | 80 | 0 |
| 7 | FLT_vs_VIV | 923 | 600 | 0.000164541 | 0.013017571 | 0.110696705 | 9 | 0 |
| 7 | GC_vs_VIV | 923 | 600 | 0.000287947 | 0.016610963 | 0.141253608 | 10 | 0 |
| 8 | FLT_vs_GC | 923 | 589 | 4.1135e-05 | 0.00670022 | 0.056976245 | 88 | 0 |
| 8 | FLT_vs_VIV | 923 | 589 | 4.1135e-05 | 0.00670022 | 0.056976245 | 16 | 0 |
| 8 | GC_vs_VIV | 923 | 589 | 0.000493624 | 0.022407294 | 0.190543509 | 8 | 0 |
| 9 | FLT_vs_GC | 923 | 580 | 4.1135e-05 | 0.007593583 | 0.064573078 | 76 | 0 |
| 9 | FLT_vs_VIV | 923 | 580 | 8.2271e-05 | 0.011390374 | 0.096859617 | 9 | 0 |
| 9 | GC_vs_VIV | 923 | 580 | 0.000904977 | 0.035294118 | 0.300128391 | 5 | 0 |
| 10 | FLT_vs_GC | 923 | 597 | 4.1135e-05 | 0.006327986 | 0.053810898 | 76 | 0 |
| 10 | FLT_vs_VIV | 923 | 597 | 4.1135e-05 | 0.006327986 | 0.053810898 | 11 | 0 |
| 10 | GC_vs_VIV | 923 | 597 | 0.000123406 | 0.009235439 | 0.078534825 | 9 | 0 |

## The three original flight–vivarium findings

All three original findings, rather than a selected substitute, are tracked. The chromosome-11 and chromosome-2 clusters have no single gene consistent with every retained intron and remain unnamed. Only the chromosome-6 cluster maps all alternatives to Tcaf1 (ENSMUSG00000036667). The two unnamed events retain 95%-gate BY significance, but none survives BY in any thinning; chromosome 11 is invalid in one thinning. Tcaf1 is valid in all ten, so its thinning null is not a denominator/coverage exclusion.

| cluster_id | contrast | n_replicates | n_measurement_valid | n_BY_positive | min_BY_q | max_BY_q | valid_dominance95 | BY_q_dominance95 | BY_positive_dominance95 | gene_id | gene_symbol | BY_q_full_3E |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 11:2:109886437:end | FLT_vs_VIV | 10 | 9 | 0 | 0.123650575 | 1.0 | True | 0.038743847 | True | unassigned | unassigned | 0.032286539 |
| 2:1:55311420:start | FLT_vs_VIV | 10 | 10 | 0 | 0.050978746 | 0.138370882 | True | 0.038743847 | True | unassigned | unassigned | 0.032286539 |
| 6:2:42663893:start | FLT_vs_VIV | 10 | 10 | 0 | 0.069185441 | 0.939244772 | True | 0.053810898 | False | ENSMUSG00000036667 | Tcaf1 | 0.045050985 |

Tcaf1 flight–ground retains q=0.038744 under the 95% gate, while flight–vivarium has q=0.053811. Across half-depth replicates, flight–ground q ranges 0.056976–0.167642 and flight–vivarium q 0.069185–0.939245. All 20 contrast–replicate results are in Tcaf1_thinning_results_all20.csv. These modest differences from the 0.05 threshold are reported quantitatively; they are not recategorized as robust significance.

## Whole-animal uncertainty

Two thousand resamples use the original group sizes (9/9/8), sampling animals with replacement within arms using seed 2026101520. Each resample applies identical animal weights across all clusters and comparisons, preserving junction/gene overlap. All 2,769 local distance intervals and 10,128 signed alternative-composition intervals are saved. The bootstrap NPZ retains every effect and resampling weight. These are pointwise approximate uncertainty summaries, without simultaneous correction or a new discovery family.

Nonnegative distance percentile intervals show upward bias: global percentile lower endpoints exceed the corresponding observed distances. They therefore cannot be used as evidence of significance. The preferred displayed distance interval projects a 95th-percentile centered-bootstrap centroid-error radius using the reverse triangle inequality, clipped to the [0,1] distance range. This is a bootstrap confidence-ball-derived approximation, not a finite-sample confidence guarantee. It includes zero for the two vivarium global comparisons; the flight–ground interval excludes zero. The exact association tests and uncertainty intervals address different reference constructions.

| contrast | observed_distance | percentile_low | percentile_high | centered_error_radius95 | confidence_ball_distance_low | confidence_ball_distance_high |
| --- | --- | --- | --- | --- | --- | --- |
| FLT_vs_GC | 0.083129639 | 0.08452631 | 0.110766844 | 0.061647427 | 0.021482211 | 0.144777066 |
| FLT_vs_VIV | 0.066147417 | 0.073153895 | 0.099519585 | 0.067605407 | 0.0 | 0.133752824 |
| GC_vs_VIV | 0.065514911 | 0.069246273 | 0.100642973 | 0.067728743 | 0.0 | 0.133243655 |

Tcaf1 intervals for all three comparisons, including the null ground–vivarium comparison, are:

| cluster_id | contrast | observed_distance | percentile_low | percentile_high | centered_error_radius95 | confidence_ball_distance_low | confidence_ball_distance_high | bootstrap_replicates |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 6:2:42663893:start | FLT_vs_GC | 0.127635615 | 0.096473648 | 0.160620583 | 0.03923117 | 0.088404445 | 0.166866785 | 2000 |
| 6:2:42663893:start | FLT_vs_VIV | 0.101115875 | 0.057486877 | 0.155188833 | 0.050949222 | 0.050166653 | 0.152065096 | 2000 |
| 6:2:42663893:start | GC_vs_VIV | 0.039801261 | 0.016828524 | 0.08238968 | 0.055029907 | 0.0 | 0.094831168 | 2000 |

Positive lower endpoints in the pointwise Tcaf1 bootstrap intervals do not restore simultaneous local discovery or override its failed depth sensitivity. Signed native-composition intervals are descriptive and not multiplicity-controlled assertions. Existing single-animal omission effects are independently crosschecked for every one of 52 omissions and retained without another hypothesis-testing family.

## Aggregate global sensitivity

All three 95%-gate comparisons survive Holm3:

| replicate | contrast | n_clusters | statistic | hellinger_centroid_distance | n_ge_observed | n_assignments | exact_p | Holm_p_3 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| dominance95 | FLT_vs_GC | 607 | 0.011512863 | 0.07587115 | 2 | 48620 | 4.1135e-05 | 0.000123406 |
| dominance95 | FLT_vs_VIV | 607 | 0.006733537 | 0.058023861 | 18 | 24310 | 0.000740436 | 0.001480872 |
| dominance95 | GC_vs_VIV | 607 | 0.007124159 | 0.05968316 | 73 | 24310 | 0.003002879 | 0.003002879 |

Every thinning retains all three aggregate associations. Each uses 580–600 valid frozen clusters, so these results compare the well-measured subset after depth perturbation; their effect magnitudes must not be presented as estimates on a universally complete set of 923 thinned measurements. All 30 exact results are supplied, with ranges below:

| contrast | n_thinnings | minimum_effect | maximum_effect | minimum_Holm_p | maximum_Holm_p | n_Holm_positive |
| --- | --- | --- | --- | --- | --- | --- |
| FLT_vs_GC | 10 | 0.080886491 | 0.084981264 | 0.000123406 | 0.000123406 | 10 |
| FLT_vs_VIV | 10 | 0.062644724 | 0.066346813 | 0.001892225 | 0.005347594 | 10 |
| GC_vs_VIV | 10 | 0.063062761 | 0.065241866 | 0.002303579 | 0.008391608 | 10 |

| replicate | seed | probability | n_clusters_valid | original_total_counts | thinned_total_counts |
| --- | --- | --- | --- | --- | --- |
| 1 | 2026101501 | 0.5 | 587 | 595410122 | 297686882 |
| 2 | 2026101502 | 0.5 | 591 | 595410122 | 297725590 |
| 3 | 2026101503 | 0.5 | 599 | 595410122 | 297715948 |
| 4 | 2026101504 | 0.5 | 598 | 595410122 | 297719993 |
| 5 | 2026101505 | 0.5 | 590 | 595410122 | 297713146 |
| 6 | 2026101506 | 0.5 | 584 | 595410122 | 297691562 |
| 7 | 2026101507 | 0.5 | 600 | 595410122 | 297699974 |
| 8 | 2026101508 | 0.5 | 589 | 595410122 | 297730360 |
| 9 | 2026101509 | 0.5 | 580 | 595410122 | 297699583 |
| 10 | 2026101510 | 0.5 | 597 | 595410122 | 297728337 |

## Dispersion and studentization

Within-arm dispersion is the equally weighted mean-cluster squared sqrt-space radius divided by two, using the sample-size-unbiased n/(n−1) adjustment. Flight and ground hardware have similar point estimates (0.011828 and 0.011740); vivarium is higher (0.014182). The flight/ground variance ratio is 1.0075; the ratios versus vivarium are 0.8340 and 0.8278. The reference dispersion label enumerations recompute group centroids for every assignment, and their Holm values are 0.9588, 0.5672 and 0.5672. This is no evidence that variances are equal: small samples, differences in location, and the reference construction limit inference.

| group | n_animals | unbiased_mean_squared_hellinger_radius | bootstrap_percentile_low | bootstrap_percentile_high |
| --- | --- | --- | --- | --- |
| FLT | 9 | 0.011827912 | 0.008265427 | 0.012321501 |
| GC | 9 | 0.011739777 | 0.008222575 | 0.012024637 |
| VIV | 8 | 0.014181734 | 0.009404851 | 0.014515334 |

A Welch-type centroid ratio T/(V_A/n_A+V_B/n_B) recomputes numerator and variance denominator for every assignment. All three studentized centroid references survive Holm (0.0001234, 0.0016454, 0.0029206). Studentization provides a sensitivity calculation, not a guarantee of small-sample robustness to covariance differences. Exact label enumeration is exact arithmetic under exchangeability of complete distributions. When the groups have different locations, the dispersion permutation is not an exact test of equal dispersion alone; when covariances differ, label exchangeability for the centroid permutation is also not assured. Nonsignificant dispersion cannot establish independence from dispersion, and these diagnostics do not remove date confounding. All cluster/arm descriptive dispersions are retained.

| contrast | variance_A_squared_hellinger | variance_B_squared_hellinger | variance_ratio_A_to_B | absolute_variance_difference | dispersion_exact_reference_p | studentized_centroid_statistic | studentized_exact_reference_p | n_assignments | dispersion_Holm_3 | studentized_Holm_3 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| FLT_vs_GC | 0.011827912 | 0.011739777 | 1.007507335 | 8.8134e-05 | 0.958782394 | 2.638987317 | 4.1135e-05 | 48620 | 0.958782394 | 0.000123406 |
| FLT_vs_VIV | 0.011827912 | 0.014181734 | 0.834024337 | 0.002353823 | 0.192554504 | 1.417421821 | 0.000822707 | 24310 | 0.567174002 | 0.001645413 |
| GC_vs_VIV | 0.011739777 | 0.014181734 | 0.827809692 | 0.002441957 | 0.189058001 | 1.394869438 | 0.002920609 | 24310 | 0.567174002 | 0.002920609 |

## Independent verification and reproducibility

200 checks pass, including independent reconstruction of every canonical denominator and alternative, every full seeded thinning matrix, reuse of every overlapping junction, all 30,459 new local test slots (including invalid p=1), every global exact reference, and SciPy BH/BY over the complete families. Independent code enumerates group-A sums and their pooled complements for local tests rather than importing the analysis routine; global and dispersion checks use direct within/across kernel sums. Bootstrap weights, distances/errors, all local intervals and all signed alternative intervals are checked; all inherited omission effects are crosschecked. The verification file lists detailed numerical assertions. No synthetic animal or library replicate is introduced.

All generated files are in this extension. README supplies package-relative offline commands. Full thinned matrices, bootstrap replicates, complete null/positive results, original-hit tracking, report and editable figure are retained and hashed in MANIFEST.csv. Figure specification is fixed 183×157 mm, Arial/Helvetica at 8 pt, minimum 1 pt lines, lowercase a)–d) panels, editable vector PDF/SVG and 600 dpi PNG.

## Interpretation limits

These analyses reuse the same 26 bulk cerebellar animals and previously selected local boundaries; ten count perturbations and 2,000 resamples are not new cohorts. The count thinning assesses read-depth sensitivity conditional on the deposited matrix and does not reproduce alignment, mapping, batch or library preparation uncertainty. Exact local arithmetic and bootstrap intervals do not establish complete isoforms, protein changes, cell-intrinsic changes, injury or a radiation mechanism. The flight condition is perfectly aliased with collection on 14 January 2021 versus both ground arms on 17 January 2021. No statistical robustness calculation can separate this alias or support a causal flight claim. The stronger conclusion is aggregate condition/date-associated local composition variation with depth-sensitive individual events. No outreach was performed.
