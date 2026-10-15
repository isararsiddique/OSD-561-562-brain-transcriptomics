# Identification and comparator sensitivity

Author: Isarar Siddique

Affiliation:

Email:

15 October 2026. Retrospective analysis lock plus a transparently logged U-statistic uncertainty amendment. This extension quantifies existing design limitations; it does not claim a causal microgravity effect or biological replication.

## Main findings

The flight indicator and early collection-date indicator are identical across all 26 animals. The design [intercept, flight, ground hardware, early date] has rank 3 for 4 columns, smallest singular value 4.32×10⁻¹⁶ and exact null vector (0,1,0,−1). Adding any constant to the flight coefficient while subtracting it from the date coefficient leaves fitted values exactly unchanged. No increase in algorithmic complexity or precision resolves this lack of identification.

Nevertheless, the observed aggregate profile is insensitive to the two measured quality covariates. Simultaneously adjusting for RIN and log2 summed read depth changes the multiway FLT–GC norm from 0.08313 to 0.08337 (+0.291%) and FLT–VIV from 0.06615 to 0.06662 (+0.707%). The corresponding adjusted/unadjusted vector cosines are 0.99859 and 0.99130. Binary FLT–GC changes by +0.009% and FLT–VIV by −0.053%. Thus measured RIN and library depth do not materially explain the observed centroid pattern under this linear sensitivity model. Adjustment can control measured covariates but cannot separate flight from date; RIN and read depth may also be downstream of treatment or technical processing.

An additive unmeasured-bias vector with norm equal to the observed FLT–GC or FLT–VIV difference can completely erase that difference. Relative to the same-date GC–VIV norm, the plug-in tipping multipliers are 1.269 and 1.010 in the multiway frame; binary values are 1.361 and 0.936. These are quantitative assumption thresholds, not estimates of actual date effects. Under the conservative joint bootstrap sampling allowance, a positive lower norm bound persists only below 0.094× the GC–VIV benchmark for multiway FLT–GC, and it is already zero at λ=0 for FLT–VIV. The latter is uncertainty, not a contradiction of the exact global permutation test, which targets a different null under exchangeability.

Naive shared-control geometry is materially affected by finite-sample noise. In multiway Hellinger coordinates, the FLT–VIV / GC–VIV inner product is 0.000879 before and −0.000894 after distinct-animal moment correction (approximate centered basic interval −0.003323 to 0.000545). The corrected binary inner product is −5.545 percentage-points squared (approximate interval −13.817 to −0.467). These exploratory moment estimates do not establish an antiparallel biological axis: only 45.47% of multiway and 37.39% of binary centered bootstrap Gram matrices admit the required nonnegative two-dimensional geometry. Corrected cosines are ratios of unbiased moment estimates, not unbiased or stable estimators, and are not used as findings.

## Mathematics and estimation

For cluster e and animal i, let p_ie be the frozen vector of junction proportions. Set X_i=concat_e sqrt(p_ie)/sqrt(2E), E=923. Then ||mean_A X−mean_B X|| is the root mean squared Hellinger centroid distance across clusters. For the binary sensitivity frame set X_i=100 f_i/sqrt(5249), giving the RMS difference in percentage points. Samples are 9 FLT, 9 GC and 8 VIV animals; technical replicates are collapsed first. The canonical raw junction matrix was checked independently against every selected multiway junction count.

For group g, the plug-in centroid Gram diagonal is G_gg=||mean_g X||². The distinct-animal estimator replaces this with G_U,gg=[n_g(n_g−1)]⁻¹ sum_(i≠j in g) X_iᵀX_j = G_gg−tr(S_g)/n_g. Cross-group entries remain G_gh=mean_g Xᵀmean_h X. Contrast squared norms and inner products are linear combinations of G_U. These linear moments are unbiased for population-centroid moments under independent animals and a fixed feature frame; the matrix need not be positive semidefinite. The feature frame is conditional on the deposited cohort and its pre-existing measurement gates.

The additive sensitivity model is delta=theta+b, where theta is the residual contrast after an assumed bias vector, not an identified microgravity effect. Under ||b||≤B, the ambient Euclidean bounds are max(0,||delta||−B)≤||theta||≤||delta||+B. They are attained by choosing b parallel or antiparallel to delta. Set B=λ||c|| with c=mean_GC X−mean_VIV X solely as a transparent scenario benchmark. GC and VIV were collected on the same date; their contrast provides no estimate of the three-day flight-versus-ground batch/date bias. Assuming B is small does not remove radiation, habitat, re-entry or other co-exposures from an interpreted flight effect.

Joint within-group resampling uses 10,000 PCG64 count-weight bootstrap draws, seed 2026101531, shared across frames and contrasts. Let r be the 95th percentile of max_k ||delta*_k−delta_k|| across all three contrasts. An approximate conservative lower bound with sampling allowance is max(0,||delta||−r−λ(||c||+r)). It is conditional on the chosen additive bound and the animal-resampling assumptions and is not a causal confidence interval. Existing global permutation results are retained unchanged.

For U-statistic moment uncertainty, ordinary empirical bootstrap percentiles require care: E*(G_U*)=G, the plug-in Gram, exactly. We use centered errors G_U*−G and basic intervals [G_U,contrast−q97.5(error), G_U,contrast−q2.5(error)] for the four linear moments. This was documented in the amendment after the initial diagnostic run. No intervals for corrected nonlinear cosine, position or orthogonal-fraction ratios are interpreted. Raw plug-in bootstrap percentiles are saved only as descriptive resampling ranges; upward norm bias makes them unsuitable for claiming nonzero signal. The correction addresses animal sampling noise and cannot remove batch confounding.

Measured-quality sensitivity fits X_i=alpha+beta_F I_FLT+beta_G I_GC+gamma_R RIN_z+gamma_D log2(depth)_z+epsilon_i by multivariate OLS, alongside group-only and each single-covariate model. RIN is the mean of library RIN values after technical-replicate collapse; depth is their sum. Covariates are standardized in the full sample. The full design has rank 5, 21 residual degrees of freedom and maximum leverage 0.419. All 10,000 within-arm bootstrap refits retain full rank. Output tables include all three contrasts under all four models; no covariate-adjusted significance claim is introduced.

## Complete summary tables

### Bias tipping values

| frame | contrast | observed_norm | control_benchmark_norm | plugin_tipping_lambda | joint_error_radius95 | conservative_tipping_lambda |
| --- | --- | --- | --- | --- | --- | --- |
| multiway_Hellinger | FLT_vs_GC | 0.0831296 | 0.0655149 | 1.26887 | 0.0703452 | 0.0941 |
| multiway_Hellinger | FLT_vs_VIV | 0.0661474 | 0.0655149 | 1.00965 | 0.0703452 | 0 |
| binary_RMS_pp | FLT_vs_GC | 5.54417 | 4.07456 | 1.36068 | 3.97671 | 0.194685 |
| binary_RMS_pp | FLT_vs_VIV | 3.81445 | 4.07456 | 0.936163 | 3.97671 | 0 |

### All measured-QC sensitivity models

| frame | model | contrast | norm_ratio | cosine_with_unadjusted | ratio_bootstrap_low | ratio_bootstrap_high |
| --- | --- | --- | --- | --- | --- | --- |
| multiway_Hellinger | group_only | FLT_vs_GC | 1 | 1 | 1 | 1 |
| multiway_Hellinger | group_only | FLT_vs_VIV | 1 | 1 | 1 | 1 |
| multiway_Hellinger | group_only | GC_vs_VIV | 1 | 1 | 1 | 1 |
| multiway_Hellinger | group_RIN | FLT_vs_GC | 1.00118 | 0.99886 | 0.978363 | 1.05208 |
| multiway_Hellinger | group_RIN | FLT_vs_VIV | 1.00712 | 0.991288 | 0.961153 | 1.13263 |
| multiway_Hellinger | group_RIN | GC_vs_VIV | 1.00169 | 0.99732 | 0.952836 | 1.08613 |
| multiway_Hellinger | group_log2_depth | FLT_vs_GC | 1.00018 | 0.999998 | 0.983124 | 1.06795 |
| multiway_Hellinger | group_log2_depth | FLT_vs_VIV | 0.997727 | 0.999355 | 0.964181 | 1.09349 |
| multiway_Hellinger | group_log2_depth | GC_vs_VIV | 1.00204 | 0.999249 | 0.958567 | 1.09364 |
| multiway_Hellinger | group_RIN_log2_depth | FLT_vs_GC | 1.00291 | 0.998591 | 0.976247 | 1.10476 |
| multiway_Hellinger | group_RIN_log2_depth | FLT_vs_VIV | 1.00707 | 0.991299 | 0.951293 | 1.19086 |
| multiway_Hellinger | group_RIN_log2_depth | GC_vs_VIV | 1.0024 | 0.997386 | 0.940242 | 1.12882 |
| binary_RMS_pp | group_only | FLT_vs_GC | 1 | 1 | 1 | 1 |
| binary_RMS_pp | group_only | FLT_vs_VIV | 1 | 1 | 1 | 1 |
| binary_RMS_pp | group_only | GC_vs_VIV | 1 | 1 | 1 | 1 |
| binary_RMS_pp | group_RIN | FLT_vs_GC | 0.997557 | 0.999106 | 0.972652 | 1.04314 |
| binary_RMS_pp | group_RIN | FLT_vs_VIV | 0.999548 | 0.990775 | 0.951686 | 1.11957 |
| binary_RMS_pp | group_RIN | GC_vs_VIV | 1.00523 | 0.997607 | 0.949747 | 1.07243 |
| binary_RMS_pp | group_log2_depth | FLT_vs_GC | 1.00034 | 0.999999 | 0.976256 | 1.07613 |
| binary_RMS_pp | group_log2_depth | FLT_vs_VIV | 0.994515 | 0.999452 | 0.957099 | 1.09296 |
| binary_RMS_pp | group_log2_depth | GC_vs_VIV | 1.0044 | 0.99945 | 0.958322 | 1.10284 |
| binary_RMS_pp | group_RIN_log2_depth | FLT_vs_GC | 1.00009 | 0.99894 | 0.967801 | 1.10672 |
| binary_RMS_pp | group_RIN_log2_depth | FLT_vs_VIV | 0.999468 | 0.990784 | 0.944827 | 1.19151 |
| binary_RMS_pp | group_RIN_log2_depth | GC_vs_VIV | 1.00703 | 0.997649 | 0.936906 | 1.12825 |

### Noise-corrected linear moments

| frame | metric | estimate | interval_low | interval_high |
| --- | --- | --- | --- | --- |
| multiway_Hellinger | FG_squared | 0.0042919 | 0.00137593 | 0.00634404 |
| multiway_Hellinger | FV_squared | 0.00128855 | -0.00162488 | 0.003138 |
| multiway_Hellinger | GV_squared | 0.00121507 | -0.00215201 | 0.00333857 |
| multiway_Hellinger | FV_dot_GV | -0.000894143 | -0.00332283 | 0.000545376 |
| binary_RMS_pp | FG_squared | 22.6278 | 11.0121 | 31.1255 |
| binary_RMS_pp | FV_squared | 4.75575 | -5.79446 | 11.2126 |
| binary_RMS_pp | GV_squared | 6.78115 | -6.21638 | 15.0689 |
| binary_RMS_pp | FV_dot_GV | -5.54545 | -13.8167 | -0.467262 |

## Manuscript-ready compact result

The deposited design was rank deficient when flight and collection date were modeled separately (rank 3 for 4 columns), confirming that date adjustment cannot identify an isolated flight effect. The aggregate junction profile was nevertheless insensitive to measured RIN and sequencing depth: simultaneous adjustment changed the multiway FLT–GC and FLT–VIV centroid norms by +0.29% and +0.71%, with adjusted/unadjusted vector cosines of 0.999 and 0.991. Under an explicit additive bias model, complete erasure of the observed differences required a bias norm equal to 1.269 and 1.010 times the same-date GC–VIV contrast, respectively. GC–VIV was used solely as a hypothetical scale benchmark, not to estimate collection-date effects. Joint animal-bootstrap uncertainty substantially weakened these sensitivity bounds. Distinct-animal centroid moments also showed that shared-control sampling noise affected apparent profile geometry; corrected nonlinear directions were unstable and were not interpreted biologically. These results support technical-covariate robustness of the observed association while quantifying the assumptions and sampling uncertainty that prevent causal attribution.

## Limitations and appropriate use

This is an additional retrospective, cohort-conditional sensitivity family. There is no prospective error guarantee across all analyses, no unmeasured-confounding estimate, no causal identification, no independently replicated junction mechanism and no claim that the bias bound is biologically plausible. Animals may share unrecorded cage-level dependencies; no usable cage randomization factor was available, so animal-bootstrap intervals assume independence at the deposited animal unit. The rank deficiency is exact and assumption-free; bootstrap intervals and noise corrections are conditional on sampling assumptions. These distinctions should be retained when integrating the figure or text into the manuscript.
