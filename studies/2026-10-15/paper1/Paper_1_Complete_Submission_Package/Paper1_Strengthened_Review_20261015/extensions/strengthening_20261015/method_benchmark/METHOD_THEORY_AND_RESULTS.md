# Sensor-based prediction, gain constraints and matched-grid results

This strengthening analysis adds a concrete positive result: a gain-constrained multivariate sensor regression can improve held-out outside-sensor expression prediction while providing a numerical amplification guarantee. Its mean prediction error was below the training-fitted nuisance baseline in all three original ground settings and all four preselected external datasets. The corresponding ordinary ridge regression frequently predicted more accurately at substantially larger operator gains. No result establishes that an exposure-associated expression pattern is an artifact or that removing it recovers causal microgravity biology.

The analysis proceeds from a nine-method benchmark to a certified convex fit, followed by a matched-grid comparison that removes a hyperparameter-search difference. All complete result tables, failures and configurations remain available. Results below are descriptive held-out prediction outcomes, with no split-independent confidence or superiority claim.

## 1. Data, estimand and fair prediction design

The original benchmark uses all 362 existing ground holdouts: 81 each in RRRM-2 cerebellum and hippocampus, plus 200 prespecified RR-23 splits. These overlapping splits use 30 unique ground source animals. The two core regions share many animals. They are not 362 independent experiments.

Every method uses the same available inherited Sensor31 genes as predictors and all remaining training-retained genes as multivariate responses. Canonical Ensembl genes must have training mean count at least ten. Geometric-mean normalization references use positive training genes only. Each held-out sample's median-ratio factor is then computed against that fixed reference. Training-fitted intercept/age effects are removed from both predictors and responses; predictor RMS scaling is estimated within training. Outcome genes are not standardized individually.

A normalized held-out sensor value is available to the predictor. Held-out outside-sensor responses never train a score, loading, component, regression coefficient or hyperparameter. As in standard whole-transcriptome normalization, the held-out sample's library factor uses its counts across the training reference genes. This is a full-RNA-seq prediction workflow, rather than a sensor-only diagnostic assay that can be normalized without the rest of the transcriptome.

Every inner leave-one-training-unit-out fold rebuilds filtering, normalization, nuisance coefficients, predictor scaling, supervised sensor weights and regression estimates. Tuning minimizes summed outside-sensor prediction SSE over the same inner targets. Varying inner gene universes arise from the training-only filter and are shared by all compared methods within that fold. No outer response selects a hyperparameter.

The training cache is keyed by dataset and sorted training units. It stores training-derived references and estimates; per-target factors are evaluated separately. Direct preprocessing comparison, held-out count perturbation and sorted-key algebra checks passed. The score, cap-two score and earlier grid-CV scalar results reproduced all original outer ratios to at most 8.88 × 10⁻¹⁶.

## 2. Original nine methods

The fixed protocol (`protocol.json`) retained all of the following methods:

1. **Intercept/age baseline:** predict only training-fitted nuisance means.
2. **Original score:** regress genome-wide residual expression on a supervised weighted Sensor31 score.
3. **Gain-two score:** bound the perpendicular loading component while retaining exact weighted-score removal.
4. **Existing grid-CV scalar:** select subtraction strength from 0, 0.25, 0.5, 0.75 and 1.
5. **Risk-CV anchored loading:** use a continuous inner-CV risk estimate while keeping exact weighted-score removal.
6. **Ridge CV:** regress outside-sensor expression on all available scaled Sensor31 predictors.
7. **PCR CV:** learn principal components from training sensor predictors and regress outcomes on those scores.
8. **PLS2 CV:** use standard multivariate covariance-maximizing predictor/response deflation learned from training data.
9. **Gain-two clipped ridge CV:** clip the fitted raw-unit coefficient singular values inside every candidate inner fit before evaluating its prediction.

Ridge candidates were baseline, 1000, 100, 10, 1, 0.1, 0.01 and 0.001. PCR/PLS candidates were zero through four components; rank limits and zero-variance predictors are retained in output diagnostics. PCR is a sensor-predictor PCA regression; it is not labelled SVA or RUV. PLS2 uses training outcome information as intended by multivariate PLS, while outer outcomes remain unavailable to its components.

The complete benchmark has 3,258 outer method rows and 4,706 passing numerical checks. Separate independent tests passed 369 checks: dense primal/dual ridge, explicit PCR score regression, PLS2 against scikit-learn NIPALS, spectral clipping, dense operator norms, exact sensor and signal identities, reference perturbation and cache algebra. The largest PLS2 prediction difference from the iterative NIPALS reference was 8.70 × 10⁻⁶ in the independent small-matrix tests. The production implementation solves each covariance direction by eigendecomposition; the reference uses iterative convergence.

### Original mean split MSE / nuisance baseline

| Method | Cerebellum RRRM-2 | Hippocampus RRRM-2 | Cerebellum RR-23 |
|---|---:|---:|---:|
| Intercept/age baseline | 1.00000 | 1.00000 | 1.00000 |
| Original score | 1.09090 | 0.96189 | 1.03488 |
| Gain-two score | 0.99681 | 0.96405 | 0.98862 |
| Grid-CV scalar | 1.00493 | 0.93137 | 0.99185 |
| Risk-CV anchored loading | 1.00490 | 0.92481 | 0.99111 |
| Ridge CV | 0.94049 | 0.86639 | 0.90414 |
| PCR CV | 1.06440 | 0.95394 | 1.02238 |
| PLS2 CV | 1.04905 | 0.96852 | 0.97747 |
| Gain-two clipped ridge CV | 0.98974 | 0.96284 | 0.96748 |

Ordinary ridge reduced mean prediction MSE by approximately 5.95%, 13.36% and 9.59%, respectively, while its mean operator gains were 28.78, 18.56 and 11.92. Clipped ridge gave smaller mean improvements of 1.03%, 3.72% and 3.25% while enforcing gain two. PCR and PLS2 failures in cerebellum and the adverse cerebellar result for scalar risk adaptation are retained.

Repeated-animal averages clarify heterogeneity. Original clipped ridge improved in all 12 cerebellar and all 12 hippocampal animals, and 16 of 17 RR-23 animals. Ordinary ridge improved in eight of 12, 12 of 12 and 16 of 17, respectively. Those are region-specific animal entries; the shared core animals must not be counted twice as independent source animals.

Floating-point ratios within 10⁻¹² of one are represented as ties when counting improvement. This prevents a baseline from being spuriously counted as improving itself.

## 3. Scalar risk adaptation: guarantee and limitation

Let the fixed unit sensor weight be w, the original regression loading be γ, and wᵀγ = 1. Write γ = w + h with wᵀh = 0. For any nonnegative λ,

\[
\gamma_\lambda=w+\lambda h,\qquad A_\lambda=I-\gamma_\lambda w^T.
\]

Then Aλ² = Aλ, wᵀAλ = 0, and for at least two retained genes,

\[
\|A_\lambda\|_2=\sqrt{1+\lambda^2\|h\|^2}.
\]

Every vector orthogonal to w is preserved exactly. A cap L ≥ 1 is obtained by restricting λ ≤ √(L²−1)/‖h‖. These identities apply to the fixed normalized-expression operator; they do not identify the scientifically correct loading.

For inner held-out responses Rᵢ and original score predictions Pᵢ, the CV objective is

\[
E(\lambda)=\sum_i\|R_i\|_F^2-2\lambda\sum_i\langle R_i,P_i\rangle
                 +\lambda^2\sum_i\|P_i\|_F^2.
\]

The continuous choice is

\[
\widehat\lambda=\operatorname{clip}\left(
\frac{\sum_i\langle R_i,P_i\rangle}{\sum_i\|P_i\|_F^2},0,1\right),
\]

with λ = 0 if the denominator vanishes. This minimizes the declared scalar inner-CV objective. Fold-omission λ ranges are supplied as sensitivity ranges, without a claim of calibrated parameter uncertainty.

Because w has no outside-sensor coordinates, γλ,outside = λ γoutside. Its outside-sensor prediction therefore equals ordinary scalar shrinkage at the same λ. Exact weighted-score removal changes the geometry inside the sensor set; it does not create an additional outside-prediction advantage. The mean cerebellar ratio 1.00490 illustrates this limit. This is a constrained interpretation of ordinary shrinkage, not an asserted novel statistical estimator.

## 4. Multi-sensor projection and the gain certificate

Order the normalized expression vector as sensor coordinates followed by outside coordinates. Let C be the sensor-by-outside prediction coefficient matrix. Define

\[
A_C=\begin{pmatrix}0&0\\-C^T&I\end{pmatrix}.
\]

The outside prediction is Cᵀ times the available sensor residual. As a correction operator, AC removes all individual sensor coordinates, is idempotent, and preserves every signal supported entirely outside those coordinates. The kernel is the span of the columns of (I, Cᵀ)ᵀ; its range is the sensor-null subspace.

For a nonempty outside space,

\[
A_CA_C^T=\begin{pmatrix}0&0\\0&I+C^TC\end{pmatrix},\qquad
\|A_C\|_2=\sqrt{1+\|C\|_2^2}.
\]

Thus gain L is equivalent to the raw-unit coefficient constraint ‖C‖₂ ≤ √(L²−1). Clipping every singular value of an unconstrained C to that limit is the unique closest coefficient under Frobenius distance, by unitary invariance and Euclidean projection of its singular values. It need not minimize prediction loss.

For a biological vector b with sensor component bS,

\[
A_Cb-b=\begin{pmatrix}-b_S\\-C^Tb_S\end{pmatrix},\qquad
\|A_Cb-b\|\le L\|b_S\|.
\]

Every sensor-disjoint vector is preserved, but genuine biology in sensor coordinates is removed and can alter outside coordinates. The multi-sensor condition is stronger than the scalar condition wᵀb = 0.

For any positive-semidefinite covariance Σ,

\[
\lambda_{\max}(A_C\Sigma A_C^T)\le L^2\lambda_{\max}(\Sigma).
\]

At L = 2 this permits up to fourfold directional variance. It is an amplification ceiling, not a theorem of variance contraction, nuisance elimination or biological FDR control. Units, retained gene universe and normalization metric determine the numerical norm.

## 5. Convex constrained regression differs from Euclidean clipping

The strengthened protocol (`strengthened_protocol.json`) directly minimizes

\[
\min_{\|C\|_2\le\sqrt3}
\frac{\|Y-SC\|_F^2}{2n}
+\frac\alpha2\|DC\|_F^2,
\]

where S contains raw-unit training sensor residuals, Y contains training outside residuals, and D is the diagonal training predictor RMS scaling matrix. This is the same scaled-sensor ridge objective expressed in raw coefficient units, with the gain constraint inside the fit.

Set H = SᵀS/n + αD² and F = SᵀY/n. With α > 0 and positive scale entries, H is positive definite and the objective is strongly convex. Its unique constrained minimizer exists on the compact spectral ball. The unconstrained solution is C₀ = H⁻¹F. Completing the square shows that the constrained fit minimizes

\[
\tfrac12\operatorname{tr}\{(C-C_0)^TH(C-C_0)\},
\]

rather than the unweighted Frobenius distance minimized by post-fit singular-value clipping. A feasible C₀, or H proportional to the identity, is sufficient for the two constructions to agree. Correlated sensors generally create a different training-loss geometry.

At a fixed α and the same training data, the constrained minimum has objective no larger than any feasible clipped coefficient. That is a training-objective guarantee. It is not a guarantee of smaller held-out error or better biology.

### Response row-space reduction

For Y = UΣVᵀ, coefficient components orthogonal to V add loss and penalty without improving prediction; projecting them away cannot increase the spectral norm. Therefore one may optimize C = K Vᵀ over the response row space, using E = UΣ as the response. The full coefficient is reconstructed through Q = KΣ⁻¹Uᵀ and C = QY. The algorithm's numerical row rank drops eigenvalues of YYᵀ below 10⁻¹² times max(λmax, 1). The exact mathematical equivalence concerns the full nonzero row space; the production threshold is a numerical-rank approximation.

### Optimization and certificate

The production solver uses projected accelerated gradient with monotone restart. Its gradient in reduced coordinates is HK−F, with step reciprocal to the largest eigenvalue of H. Projection onto the spectral ball clips coefficient singular values.

For a feasible K, the first-order convex gap is

\[
g(K)=\langle\nabla f(K),K\rangle+
\sqrt3\,\|\nabla f(K)\|_*.
\]

This follows from the duality between the spectral and nuclear norms. Convexity implies 0 ≤ f(K)−f(K*) ≤ g(K). The declared stopping rule is g(K) ≤ 10⁻⁷ max(1, training baseline objective). All 16,084 saved fits converged, with maximum relative gap 9.99995 × 10⁻⁸. These are numerical optimality certificates conditional on the fitted training data, not statistical confidence intervals.

The independent full-coefficient CVXPY/Clarabel conic formulation passed 96 checks. Maximum objective disagreement with the production row-space solver was 2.46 × 10⁻⁹. Dense feasibility, operator norm, idempotence, sensor annihilation, sensor-disjoint signal preservation and independent gap reconstruction were checked.

The earlier SLSQP diagnostic remains available as `convex_slsqp_diagnostic.json`: four cases reached its iteration limit and five numerical checks failed. It is retained as an unsuccessful independent numerical attempt and is not used as the optimization reference. The full coefficient conic formulation is distinct from the production projected-gradient algorithm.

## 6. External datasets and strengthened outcomes

Four cohorts were selected from their design and count availability, before examining these predictive outcomes. OSD-202 contributes 40 Earth-analogue brain animals; OSD-352 contributes six RR-3 bulk-brain animals. OSD-374 contributes eight two-animal hippocampal pools and OSD-611 ten retained two-animal pools. OSD-611 has 12 count columns, but two lack the required NASA SampleTable support and remain excluded by the recorded design rule. Pool members are not treated as separate replicates.

The additional 64 outer holdouts are therefore 46 individual animals plus 18 pool-level units. They are not 64 newly generated animals. Each external unit is held out once. All methods are refitted and tuned within each external cohort; this tests algorithm behavior in independently acquired datasets, rather than transport of frozen original loadings. OSD-202 adjusts for collection time and unloading; the other added cohorts use an intercept. Their flight/radiation labels orient exposure-associated scalar scores and must not be reclassified as known nuisance labels.

The completed strengthened run has 1,002 outer rows: one constrained method for the 362 original splits plus ten methods for the 64 external holdouts. Its complete tables contain all original comparators and failures.

## 7. Matched-grid check of constraint versus clipping

The initial constrained fit searched baseline, 10, 1, 0.1 and 0.01; the original clipped comparator searched a larger eight-candidate grid. Those outer differences compare estimators with both algorithm and search-grid differences. To isolate this concern, a separately frozen amendment (`matched_grid_protocol.json`) rebuilt both ordinary ridge and clipped ridge using exactly the constrained method's five-candidate grid, with identical inner folds, preprocessing, targets, SSE weighting and per-candidate clipping.

This new run preserves 852 outer rows over all 426 original/external holdouts, plus full inner candidate tables. Earlier strengthened rows are unchanged. Twenty-eight checks compared all four nonzero α values on the first outer training set in each of the seven deposits. The constrained objective was strictly below the feasible clipped objective in every check; the smallest observed objective improvement was 36.72. Direct full-coefficient objective reconstruction agreed with the solver to at most 6.37 × 10⁻¹². Both coefficients obeyed the spectral limit within floating-point tolerance.

### Mean held-out MSE ratios with the matched grid

| Setting | Units/splits | Matched ordinary ridge | Matched gain-two clipping | Constrained gain-two fit | Constrained minus clipping |
|---|---:|---:|---:|---:|---:|
| RRRM-2 cerebellum | 81 overlapping splits | 0.94004 | 0.98958 | 0.98776 | −0.00183 |
| RRRM-2 hippocampus | 81 overlapping splits | 0.86639 | 0.96284 | 0.94835 | −0.01448 |
| RR-23 cerebellum | 200 overlapping splits | 0.90372 | 0.96748 | 0.96284 | −0.00463 |
| OSD-202 | 40 animal holdouts | 0.80280 | 0.95103 | 0.92845 | −0.02259 |
| OSD-352 | 6 animal holdouts | 1.00880 | 0.99356 | 0.98565 | −0.00791 |
| OSD-374 | 8 pool holdouts | 0.70438 | 0.96062 | 0.95013 | −0.01048 |
| OSD-611 | 10 pool holdouts | 0.97091 | 0.97382 | 0.97186 | −0.00196 |

All constrained mean ratios are below baseline one and below their matched clipped comparator. Mean baseline reductions are approximately 1.22%, 5.16%, 3.72%, 7.16%, 1.43%, 4.99% and 2.81%, respectively. The differences are descriptive; the very small RR-3 cohort and overlapping original splits limit population-level inference.

The constrained fit had lower error than matched clipping in 56/81 cerebellar splits, all 81 hippocampal splits, 185/200 RR-23 splits, 34/40 OSD-202 animals, all six OSD-352 animals, seven of eight OSD-374 pools, and three of ten OSD-611 pools. It improves a mean without uniformly winning each held-out case.

Ordinary ridge's mean gains in the matched run were 28.99, 18.56, 11.96, 30.65, 5.73, 29.93 and 27.14. Most settings illustrate the cost of the fixed gain cap: ordinary ridge achieves a larger prediction improvement at greater possible directional amplification. The matched OSD-352 ridge result also illustrates sensitivity to the search grid; its original wider-grid ratio was 0.99013, while the matched narrower-grid ratio was 1.00880. All outcomes are retained rather than choosing the favorable grid.

### Error weighting and heterogeneous outcomes

A mean of per-unit MSE ratios weights each unit's relative change equally. The pooled SSE ratio weights units by their baseline error and retained response counts. Both are supplied.

For OSD-374, the constrained mean ratio 0.95013 corresponds to pooled SSE ratio 0.97438. For OSD-611, the constrained equal-unit mean 0.97186 is lower than clipped 0.97382, but pooled ratios favor clipping: 0.98358 versus 0.98091. Its paired median constrained-minus-clipping difference is +0.00307, and constrained prediction wins only three of ten pools. Thus a universal superiority statement is contradicted by the full outputs.

Repeated-unit constrained averages in the original datasets improved in ten of 12 cerebellar animals, all 12 hippocampal animals and 16 of 17 RR-23 animals. The original clipped method was more uniformly beneficial in cerebellum (12/12), despite the constrained method's smaller overall mean error. This is a useful stability–prediction tradeoff, not a reason to hide the adverse animals.

## 8. Defensible interpretation and manuscript language

The contribution is an audited comparison of inductive sensor regression and a directly optimized amplification constraint on real transcriptomes. It strengthens the operational result beyond demonstrating failures of the original score: a feasible, numerically certified fit improves mean held-out prediction in all tested settings while exposing the accuracy cost and heterogeneity of gain control.

The optimization problem, spectral/nuclear duality and singular-value projection are standard convex mathematics. No novelty claim is made for ridge, PCR, PLS2, ordinary scalar shrinkage or these projection identities. Any novelty claim should concern the specific biological-resource benchmark, constraints, verification and observed stability tradeoffs, after comparison with prior literature.

A possible results paragraph is:

> In a fully nested Sensor31 predictor benchmark, ordinary ridge reduced mean outside-sensor prediction error by 5.9–13.4% in the three original ground settings, with mean gains of 11.9–28.8. A directly optimized gain-two constrained fit gave smaller but positive mean reductions of 1.2%, 5.2% and 3.7%. With identical hyperparameter grids and folds, it improved mean prediction relative to post-fit clipping in all seven original and added settings. All 16,084 convex fits met the declared dual-gap tolerance, and full-coefficient conic verification passed. Individual and error-weighted comparisons retained exceptions, including the ten-pool OSD-611 cohort. These results establish a practical prediction–stability tradeoff and require independent interpretation of exposure biology.

The retained counterexamples are essential: training-objective optimality is not test-error dominance; a norm cap does not prove correct nuisance attribution; a preserved sensor-disjoint vector does not imply calibrated differential-expression inference. The external radiation/flight predictors can carry genuine biology. Their success in predicting expression validates an estimator's predictive behavior, not removal of artifacts.

## 9. Complete files and reproducibility

- **Original design and outputs:** `protocol.json`, `run_method_benchmark.py`, `outer_split_results.csv`, `nested_CV_candidates.csv`, `nested_CV_risk_fold_components.csv`, `outer_animal_prediction_errors.csv`, `animal_average_prediction_summary.csv`, `method_summary.csv`, `verification_checks.csv`, `run_summary.json`.
- **Independent original-method verification:** `verify_methods.py`, `independent_method_verification.csv/.json`; representative factorization snapshots for the three original deposits.
- **Convex method and full added benchmark:** `constrained_regression.py`, `run_strengthened_methods.py`, `strengthened_protocol.json`, `strengthened_outer_results.csv`, `strengthened_unit_errors.csv`, `strengthened_inner_CV.csv`, `strengthened_split_roles.csv`, `strengthened_summary.csv`, `convex_optimality_certificates.csv`, `strengthened_run_summary.json`.
- **Independent convex verification and retained failure:** `verify_constrained_regression.py`, `convex_independent_verification.json` (96 checks PASS), `convex_slsqp_diagnostic.json` (earlier failed numerical diagnostic retained).
- **Matched search-grid amendment:** `matched_grid_protocol.json`, `run_matched_grid_comparator.py`, `matched_grid_outer_results.csv`, `matched_grid_inner_CV.csv`, `matched_grid_unit_errors.csv`, `matched_grid_fixed_alpha_objectives.csv`, `matched_grid_summary.csv`, `matched_grid_run_summary.json` (852 rows, 426 holdouts, 28 fixed-α checks PASS).
- **External-source accounting:** `../external_data/cohort_inventory.csv`, `harmonized_metadata.csv`, `sample_to_source_membership.csv`, `acquisition_provenance.json`, each dataset's metadata/count snapshot and `output_manifest.json`. The selection and available-count inventories remain distinct from the number of included outer units.

The benchmark and matched-grid code require NumPy, pandas and SciPy; the scientific Python environment used here matches NumPy 2.5.3, pandas 3.0.5 and SciPy 1.18.1. Independent PLS verification additionally uses scikit-learn; convex verification uses CVXPY/Clarabel. BLAS thread counts were restricted to one. None of the reported split quantiles should be presented as confidence intervals from independent split observations.
