# Frozen retrospective comparator and identification analysis

15 October 2026. Saved before numerical analysis in this module. This is an analysis lock, not prospective registration. Author: Isarar Siddique. Affiliation and email blank.

## Questions and scope

Quantify what the three-arm design identifies, how observed profiles relate to the same-date GC–VIV direction, and what explicit bias assumptions would be needed for an interpreted residual association. No statistical adjustment can separate date from flight. Use the frozen 923-cluster square-root composition frame as primary; the frozen 5,249 binary fraction frame is a sensitivity measurement scale. No new feature selection, local discovery testing, or optimized endpoint search.

## Design

Extract the 26 unique animal identities and euthanasia/dissection dates from deposited metadata, joining the frozen sample order. Build [intercept, FLT, GC, early-date] and verify the exact null vector [0,1,0,-1] and rank. Save animal design, rank, singular values and two distinct coefficient vectors with identical fitted values.

## Geometry and finite-sample noise

For multiway features set X_i=concat(sqrt(p_ie))/sqrt(2E), giving squared centroid norm identical to mean squared Hellinger centroid distance. For binary features use X_i=100*f_ie/sqrt(E), yielding RMS percentage-point differences. Compare all three centroid vectors. Place VIV at origin and GC on the first axis: t=<mu_F-mu_V,mu_G-mu_V>/||mu_G-mu_V||²; orthogonal distance and fraction are descriptive.

Shared group sampling noise contributes to naive inner products. Form a Gram matrix of group means with diagonal corrected by tr(S_g)/n_g. This is the unbiased U-statistic Gram estimator under independent animals and a fixed feature frame. Derive corrected contrast squared norms and inner products from that matrix; corrected squared norms may be negative and will not be forced positive for inference. Use corrected geometric ratios only when all required squared norms are positive, report bootstrap invalid fraction. These do not remove batch bias.

Generate 10,000 joint within-arm animal bootstrap replicates (PCG64 seed 2026101531) for both frames using the same count weights. Save all replicate diagnostics. Percentile intervals for smooth geometry are approximate and exploratory. For contrast norms use joint error balls: the 95th percentile of the maximum over the three contrasts of ||delta*−delta|| supplies one simultaneous bootstrap radius per frame. Do not use a strictly positive percentile interval of a norm to claim significance.

## Explicit bias-bound sensitivity

Write delta_observed=theta+b for either FLT contrast. Under the assumed bound ||b||≤B, norm(theta) is within [max(0,||delta||−B), ||delta||+B]. B is not identified. Show B=lambda*||mu_G−mu_V||, lambda from 0 to 2 in steps 0.025, explicitly a hypothetical benchmark, not an estimate of date effects. The plug-in tipping lambda is ||delta||/||mu_G−mu_V||. Also show a conservative lower bound allowing the joint bootstrap radius r: max(0,||delta||−r−lambda*(||control contrast||+r)). This is approximate sampling sensitivity conditional on the bound, not a causal confidence interval. At unrestricted B no nonzero causal effect is identified. Bias ball is ambient Euclidean, so any additional compositional constraints could tighten it; none are assumed.

## Measured quality sensitivity

Reconstruct animal QC from 29 library records: sum read depth, mean RIN after technical-replicate sample collapse. Verify against the existing animal QC table. Fit multivariate OLS with VIV reference under four frozen models: group only; group+RIN; group+log2 summed read depth; group+both. Standardize continuous covariates once in the full 26 animals. Report three adjusted contrast norms, cosine alignment with unadjusted vectors, relative changes, rank, residual degrees of freedom and maximum leverage. Joint within-arm bootstrap refits for approximate uncertainty, excluding rank-deficient draws with their count reported. No adjusted p values or claim of causal rescue: RIN/read depth can be consequences of treatment or technical processing, and collection date remains aliased.

## Outputs and checks

Complete input hashes, reproducible script, tables, exact mathematical Methods, manuscript-ready Results, independent algebraic verifier, and a four-panel native 183-mm figure if the quantitative findings add interpretability. Vector PDF/SVG and RGB 600 dpi PNG; font at least 8 pt. Verification includes raw source design, distinct aliased fits, kernel-feature equality, direct U-statistic pairwise calculations, bias-ball extrema, QC aggregation, and bootstrap/index integrity. No guarantee of publication or forced positive finding.
