# Frozen retrospective count-model sensitivity plan

Date: 2026-10-15. Written before new benchmark results; prior manuscript outcomes are already known. This is retrospective sensitivity analysis, not preregistered confirmatory research.

## Data and immutable scope

Use all 923 previously selected RR23 multialternative clusters and all 26 biological animals, with technical libraries already summed. Use the exact frozen 3,376 alternative occurrences and their real junction counts. Retain the three FLT–GC, FLT–VIV and GC–VIV contrasts as one 2,769-test local family. Neither features, candidate choice, model choice nor thresholds will be tuned on new significance. Tcaf1 is the previously discussed coordinate-defined candidate. Animal exchangeability, and its collection-date limitation, remain exactly as in the existing analysis.

## Model and exact calibration

For counts x_i with total N_i, use the standard Dirichlet–multinomial likelihood DM(N_i, alpha), fitted to all 26 animals without reading group labels. Fit log(alpha) by bounded L-BFGS-B, from two deterministic starting concentration totals (10, 1,000) and the pooled count proportions. Bounds on each alpha are [10^-6, 10^5]. Keep the converged higher-likelihood solution; record both fit diagnostics, boundary hits, concentration A=sum(alpha), and intraclass overdispersion rho=1/(A+1). A cluster with no converged finite fit is retained with p=1; it is never dropped from correction.

At the pooled fit calculate the standard derivative of log likelihood with respect to mean-composition logits:
S_ij = alpha_j { [digamma(x_ij+alpha_j)-digamma(alpha_j)] - sum_l pi_l[digamma(x_il+alpha_l)-digamma(alpha_l)] }, pi=alpha/A.
All K score coordinates are retained and their one redundant dimension is removed numerically. For each contrast, center all pooled animal score vectors, eigendecompose their empirical covariance, and whiten eigenvalues greater than 10^-10 times the largest eigenvalue. Report every rank and condition number. The local statistic is the squared Euclidean distance between group-average whitened score vectors. Because fitting, centering and whitening do not depend on the tested labels, enumerate all fixed-size animal assignments (48,620 or 24,310) for an exact reference p value. Equality tolerance is 10^-10 relative to max(1, observed). These are permutation-calibrated DM-score sensitivity tests; they are not likelihood-ratio tests or a reimplementation of LeafCutter/rMATS. The empirically whitened score deliberately protects calibration against relying on the small-sample chi-square approximation; validity still requires animal exchangeability.

Compare the exact-calibrated DM score with (a) the frozen original square-root-composition test and (b) the multinomial-limit mean-logit score x-N*pi, passed through identical pooled whitening and exact permutation calibration. Use BY and BH across all 2,769 local tests independently for each named model. The primary benchmark is model agreement and informative failure, not an enlarged union of discoveries. An equal-cluster average of local statistics divided by rank supplies a descriptive global score; test its exact null and apply Holm across the three contrasts separately for each model. No asymptotic p values will support discoveries.

Report complete fit diagnostics, all exact p and adjusted values, native animal-average composition effect magnitudes (unchanged across models), correlations of ranks/evidence, discoveries and overlaps, rho versus depth, and all three Tcaf1 contrasts. Fits pooled across groups estimate total heterogeneity, not a biological within-group dispersion parameter or causal flight effect.

## Frozen-depth and annotation sensitivity

Repeat the DM-score benchmark for each of the 10 already saved independent half-depth full-count thinnings, keeping the original selected alternatives. Use the original gates (all-animal selected count >=20, dominance >=0.90). Failed gates or fits receive p=1 in the entire 2,769-test family. These are technical downsampling sensitivity analyses, not independent biological replications. Repeat original full-depth DM local p values under the stricter existing 0.95 dominance gate, with p=1 for failed clusters. Report all outcomes, including losses. Fit convergence is audited in all replicates. Do not alter the independent-mission conclusion or causal interpretation.

## Verification and presentation

Implement an independent checker using scalar log-likelihood and finite-difference score derivatives, direct covariance pseudoinverse statistics, direct selected permutation enumeration, count-to-composition reconstruction, and independent BY adjustment. Record SHA-256 input and plan hashes before analysis; retain all results and failures. Verify bibliographic sources from publisher or bibliographic metadata. Produce a 183-mm four-panel vector figure with >=8-pt text; numerical results are primary artifacts. No claim that stronger statistics guarantee npj acceptance.
