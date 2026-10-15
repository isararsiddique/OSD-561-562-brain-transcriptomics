# Independent scientific audit and exploratory robustness extension — 15 October 2026

The source counts, animal collapse and central saved results passed this audit. The strongest new finding is an exact, fully refitted label-enumeration validation of the ground collection-score association. Sensor31 achieved crossfitted AUC 1.000 in each core region; only the observed assignment and its label complement reached that statistic among 400 age-blocked assignments. Exact p is 0.005, BH q is 0.01125 and BY q is 0.03932 in the complete 18-test audit family. These probabilities concern observational within-age label exchangeability. They do not turn collection cohorts into randomized groups, establish independent mission replication, or validate correction of flight biology.

The full dataset remains 129 libraries, 126 animal-region observations and 83 Source Name animals. There are no newly generated animals or raw-read alignments in this extension.

## What was audited

- All bundled provenance hashes, both archived official core STAR matrices, nonnegative integer STAR values in all seven deposits, and technical-replicate sums were checked independently.
- Complete BH adjustments for the 54,989 cross-cohort replication tests and 72,023 evaluable flight-by-collection gene tests were rebuilt from their component probabilities.
- The full and metadata-excluded panel interaction coefficients were reproduced from direct count normalization and full design regression.
- Every existing ground holdout was reconstructed from source processed counts, using training-only references and direct regression. The largest difference from the supplied original-score MSE ratios was 6.66 × 10⁻¹⁶.
- The main audit passed 412 checks. A separate verifier passed ten further checks, including 486 explicit regressions covering all 81 holdouts for three selected null assignments in each region. It checked the efficient matrix implementation against full intercept/age/group regression, all-label-complement invariance and exact null-mean AUC 0.5. The largest discrepancy was 1.11 × 10⁻¹⁶.

No consequential source-count, animal-unit, coefficient or arithmetic error was found in these targets. The existing caution about metadata disagreement, overlapping splits, conditional normalization and unavailable probabilities is warranted. This audit does not refit thousands of negative-binomial models or re-run the completed simulation studies.

## Exact refitted-label validation

The analysis was specified as a post hoc robustness audit in `protocol.json` before running. For each core region, ground animals comprise six animals in each of two age blocks, with three animals per ground collection cohort in each age block. All C(6,3)² = 400 cohort-label assignments were enumerated, including the observed labels.

For every assignment, all 3⁴ = 81 balanced holdouts were rebuilt: one animal from each assigned age-by-cohort cell was held out, leaving eight training and four test animals. Gene filtering (training mean count at least ten) and geometric-mean normalization references used training counts only. Their results were cached by training animal set because these operations do not use cohort labels. The supervised, age-adjusted contrast weights were refitted for every assignment and training set. The test statistic was one mean of the 162 within-age held-out pair orderings across the 81 overlapping splits.

The efficient implementation uses the expression Gram matrix for each training/test configuration. With balanced training labels within both age blocks, the age-adjusted contrast is Y_train z / 4, where z contains the assigned −1/+1 training labels. Comparing weighted test scores therefore gives the same ordering as explicit full design regression. The separate verifier recomputed both supervised weights and score orderings directly.

All eight inherited panels plus a whole-genome contrast comparator were retained for each of the two regions, forming 18 exact tests in one audit family. A p value is the fraction of all 400 assignments with a mean AUC at least as large as observed. BH and BY adjustments use that same complete family. The label-complement symmetry explains the attainable minimum p = 2/400. The repeated holdouts are aggregated inside each permutation statistic; they are never treated as independent experimental units.

Sensor31 and rPRG were perfect in both core regions and had exact p = 0.005, BH q = 0.01125 and BY q = 0.03932. In hippocampus, dPRG, SRG, Circadian and the whole-genome comparator were also perfect and attained those probabilities. Cerebellar HSP AUC was 0.593 with exact p = 0.360; its whole-genome comparator AUC was 0.796 with p = 0.080. These complete comparisons support a reproducible collection-cohort association while limiting claims of a uniquely specific Sensor31 signature.

## Prediction-error anatomy

The original score's mean split MSE ratio was reproduced as 1.09090 in cerebellum, 0.96189 in hippocampus and 1.03488 in RR-23. The rank-one prediction risk satisfies the exact identity

MSE = mean(R²) + mean(P²) − 2 mean(RP),

where R is held-out expression residualized with the training-fitted intercept/age model and P is the training loading times the available held-out sensor score. All decomposition outputs are retained. No prediction strength was selected using outer held-out errors.

After forming expression tertiles from training mean counts only, mean split MSE ratios were:

| Deposit | Lower | Middle | Higher |
|---|---:|---:|---:|
| Cerebellum, RRRM-2 | 1.0961 | 1.0897 | 1.0650 |
| Hippocampus, RRRM-2 | 1.0030 | 0.8987 | 0.8148 |
| Cerebellum, RR-23 | 1.0448 | 1.0078 | 1.0254 |

Thus the adverse mean prediction result in cerebellum and RR-23 persists across training expression levels. This is a descriptive decomposition using the existing overlapping holdouts, without a new independent-sample confidence claim. The average fraction of outside-sensor genes with worse prediction was 56.2%, 40.3% and 54.9%, respectively.

Animal-level repeated-holdout averages are heterogeneous: score prediction worsened in six of 12 cerebellar ground animals, five of 12 hippocampal ground animals and 13 of 17 RR-23 ground animals. These are 41 region-specific animal entries involving 30 unique ground Source Name animals. They explain why a mean adverse result should not be presented as uniformly adverse for every animal. Ratio of animal-averaged errors and mean split ratios are different summaries and are explicitly distinguished in the saved tables.

## Age dependence of the Circadian candidate

Full-data direct younger and older hippocampal difference-in-differences were −0.5100 and −0.3070 log2. The equal-age-standardized estimate was −0.4085, close to the additive-age estimate −0.4053. Its HC3 normal interval was [−0.6830, −0.1341], with approximate p = 0.00353. However, its BH q was 0.11297 in the complete 64-estimand saturated-age audit family, and no saturated-age result passed BH q ≤ 0.05 (family minimum q = 0.09338).

Full-data cerebellar effects were −0.4489 and −0.1309, with equal-age mean −0.2899. These diagnostics retain effect direction while showing the uncertainty cost of permitting age-dependent group effects. Absence of a significant age interaction is insufficient evidence that the additive-age restriction is correct.

After conflict-animal exclusion, the hippocampal age-specific point estimates stayed negative (−0.4418 and −0.5992), but the older FLT ISS-T cell has one animal. The cerebellar older FLT ISS-T cell is empty, so its direct older and equal-age-standardized interactions remain unavailable. No pseudo-estimate or favorable p value was substituted for the empty cell. All eight panels and all cell means are provided.

## Suggested manuscript addition (116 words)

An exploratory exact refitted-label audit further tested ground score association without treating repeated splits as independent. All 400 age-blocked label assignments in each core region rebuilt the supervised weights across 81 balanced holdouts, using training-only filtering and normalization. Sensor31 had mean age-stratified AUC 1.000 in both regions; two assignments reached that value (exact p = 0.005; BH q = 0.01125; BY q = 0.03932 across all 18 panel/comparator-by-region tests). The rPRG panel also achieved perfect discrimination. Original loading prediction remained worse on average in every training-expression tertile in cerebellum and RR-23. These conditional observational association results strengthen the distinction between score discrimination and genome-wide loading transfer, without establishing biological recovery or independent mission replication.

A compact age qualifier can state: “Allowing all age-by-flight-by-collection terms retained negative younger and older hippocampal Circadian estimates, with equal-age mean −0.409 log2. The equal-age contrast did not pass the complete 64-estimand audit correction (q = 0.113).”

## Figure and legend

`independent_audit_panel.png` is 4500 × 3960 pixels at 600 dpi (7.5 × 6.6 inches). `independent_audit_panel.svg` retains editable Arial text and vector geometry. Both use lowercase bold panel labels and line widths of at least one point. The final artwork was visually inspected.

**Suggested Supplementary Figure 5 legend.** Independent refitted-label and age-model robustness audit. (a,b) Sensor31 mean held-out age-stratified AUC under all 400 age-blocked ground-cohort label assignments in cerebellum and hippocampus; each assignment refits supervised weights across 81 balanced animal holdouts. Red lines show observed AUC 1.000. Filtering and normalization references use training counts only and are cached independently of labels. Exact p values and BH/BY probabilities use all 18 panel/comparator-by-region tests; complete null tables for all nine scores per region are supplied. Label exchangeability is an observational assumption. (c) Original score-loading outside-sensor prediction MSE relative to the training-fitted intercept/age baseline, within tertiles defined by training mean counts. Points show split means; lines show overlapping-split 2.5th–97.5th percentile robustness ranges, not independent-sample confidence intervals. (d) Full-data Circadian younger, older and equal-age interactions from an unrestricted age-by-flight-by-collection model. Lines show approximate HC3 normal 95% intervals. None of the 64 saturated-age estimands passes BH q ≤ 0.05. Complete metadata-excluded cell means and unavailable empty-cell estimates are retained in source tables.

## Files and reproduction

Use the original scientific interpreter `research/brain_integrative_2026/.venv/bin/python`; NumPy 2.5.3, pandas 3.0.5 and SciPy 1.18.1 match the revision. Set OPENBLAS_NUM_THREADS, OMP_NUM_THREADS and VECLIB_MAXIMUM_THREADS to one for consistent resource use.

1. `run_independent_audit.py` reconstructs source accounting, exact association, age diagnostics, risk anatomy and 412 checks.
2. `verify_exact_ground_null.py` independently verifies the null-score implementation using explicit full regressions.
3. `build_audit_figure.py` builds the panel from saved complete results and checks size, Arial and vector text.

`protocol.json`, `audit_summary.json`, `exact_ground_verification_summary.json`, `figure_specification.json`, `structured_audit_results.json` and `source_file_hashes.json` provide configuration, structured results and provenance. Complete CSVs preserve all tested panels, all 400 null assignments, all 362 risk reconstructions, expression tertiles, animal summaries, age cell means and all 64 saturated estimands. Existing shared scripts, original result tables and manuscript files were not intentionally modified by this extension.
