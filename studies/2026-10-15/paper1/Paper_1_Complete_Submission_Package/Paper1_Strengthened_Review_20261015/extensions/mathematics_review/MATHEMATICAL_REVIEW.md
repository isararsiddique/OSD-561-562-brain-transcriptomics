# Independent mathematical review and amended stress analysis

This review concerns the revision's linear expression-space transformations. The bounds below apply to normalized expression and coefficient vectors, not directly to raw counts, biological truth, adjusted t statistics or causal effects. Results were verified from saved outputs.

## 1. Oblique score projection

Let \(Y\in\mathbb R^{G\times N}\), let \(w\in\mathbb R^G\) be a unit sensor-weight vector, and assume \(w^T\gamma=1\). The score operator is

\[
A=I-\gamma w^T,\qquad Y^*=AY.
\]

Then \(A^2=A\), \(w^TA=0\), \(\ker(A)=\operatorname{span}(\gamma)\), and \(\operatorname{range}(A)=w^\perp\). Therefore this operator is an oblique projection onto the sensor-score-null hyperplane along \(\gamma\). It is orthogonal only when \(\gamma=w\). The identity \(w^T\gamma=1\) follows algebraically from the ground regression \(\widehat\gamma=Y_0s/(s^Ts)\), with \(s=w^TY_0\), provided \(s^Ts>0\).

Write \(\gamma=w+h\), with \(w^Th=0\). In the orthonormal basis \(w,h/\|h\|\), the operator's nontrivial block is

\[
\begin{pmatrix}0&0\\-\|h\|&1\end{pmatrix}.
\]

For \(G\ge2\), its largest singular value is consequently

\[
\|A\|_2=\sqrt{1+\|h\|_2^2}=\|\gamma\|_2.
\]

The remaining singular values are one (multiplicity \(G-2\)) and zero. The one-dimensional exception is \(A=0\), whose norm is zero; the equality to \(\|\gamma\|\) requires the stated dimension qualification. Score cancellation does not mean cancellation of every sensor gene or reduced sensor-panel energy. The saved counterexample \(w=(1,0,0)^T\), \(\gamma=(1,2,3)^T\), and \(b=(1,0,0)^T\) has \(Ab=(0,-2,-3)^T\): the first-coordinate sensor score disappears while energy on sensor coordinates 1 and 2 increases from 1 to 4.

## 2. Bounded score operator and optimization interpretation

For an a priori gain cap \(L\ge1\), define

\[
c_L=\min\{1,\sqrt{L^2-1}/\|h\|\},\quad
\gamma_L=w+c_Lh,\quad A_L=I-\gamma_Lw^T.
\]

Use \(c_L=1\) when \(h=0\). Then \(w^T\gamma_L=1\), \(A_L^2=A_L\), \(w^TA_L=0\), and

\[
\|A_L\|_2=\|\gamma_L\|_2\le L.
\]

Moreover \(\gamma_L\) is the unique Euclidean closest admissible loading to \(\gamma\):

\[
\gamma_L=\operatorname*{argmin}_{g:\ w^Tg=1,\ \|g\|\le L}\|g-\gamma\|^2.
\]

This follows by writing \(g=w+k\), \(k\perp w\), and projecting \(h\) onto the ball \(\|k\|\le\sqrt{L^2-1}\). This is a geometric safeguard. The selected value \(L=2\) is not an estimate of the biologically optimal correction strength. It must be described as selected in advance for the sensitivity analysis. The numerical scale depends on the gene universe and expression metric.

For any positive-semidefinite covariance \(\Sigma\),

\[
\lambda_{\max}(A_L\Sigma A_L^T)\le L^2\lambda_{\max}(\Sigma).
\]

This bounds the largest possible directional variance, not gene-wise variance contraction. It allows a fourfold directional variance increase at \(L=2\). The same bound holds for squared Euclidean energy ratios for arbitrary coefficient vectors. If sensors are supported entirely on a set \(S\), and the outside set \(O\) is nonempty, \(\|P_OA\|_2=\sqrt{1+\|\gamma_O\|^2}\). Thus amplification outside the sensor is possible even though \(w^TAY=0\).

## 3. What biology is preserved and what remains unidentified

For any fixed biological effect vector \(b\),

\[
A_Lb=b-\gamma_L(w^Tb),\quad
\|A_Lb-b\|=\|\gamma_L\|\,|w^Tb|\le L|w^Tb|.
\]

Every vector in \(w^\perp\) is preserved exactly. This includes any effect supported entirely outside the sensor. Biological effects with a nonzero weighted sensor score are altered. Energy retention alone is insufficient to show preservation: energy can increase while direction or individual gene effects are distorted. Report relative distortion \(\|A_Lb-b\|/\|b\|\), and, where useful, angle, in simulations with known \(b\).

Suppose the target group-mean difference is \(\beta=b+\delta u\), where \(u\) is a nuisance direction. Then

\[
A_L\beta-b=-\gamma_L(w^Tb)+\delta\{u-\gamma_L(w^Tu)\}.
\]

The first term is loss/distortion of biology and the second is residual nuisance. Ground association alone does not distinguish these components or guarantee that target coupling equals ground coupling. If the original loading is exactly the target nuisance direction, \(u=(w^Tu)\gamma\), bounding it deliberately leaves residual nuisance:

\[
A_Lu=(w^Tu)(1-c_L)h.
\]

Thus a strong norm cap can sacrifice nuisance removal even under perfect transfer. The limit \(L=1\) becomes the orthogonal projection \(I-ww^T\) and leaves all outside-sensor coordinates unchanged.

## 4. High-dimensional loading uncertainty: an explanatory special case

Assume a fixed unit \(w\), selected independently of the data, isotropic Gaussian residual expression, and residual sample rank \(r=N-\operatorname{rank}(X_0)\). Population \(\gamma=w\) has norm one. After rotation to \(w=e_1\), the score norm is \(s^Ts\sim\chi_r^2\), while each of the \(G-1\) estimated perpendicular loading coordinates has conditional variance \(1/(s^Ts)\). Therefore, for \(r>2\),

\[
\mathbb E\|\widehat\gamma\|^2=1+\frac{G-1}{r-2}.
\]

The perpendicular squared norm has law \(\chi^2_{G-1}/\chi^2_r\), with independent numerator and denominator. For \(r>4\), its variance is

\[
\operatorname{Var}(\|\widehat\gamma\|^2)=
\frac{2(G-1)(G+r-3)}{(r-2)^2(r-4)}.
\]

This demonstrates a possible amplification mechanism from high dimensionality and small ground sample size without transferable outside-sensor nuisance. It is not the exact sampling law of the present estimator, whose \(w\) is estimated from the same ground data and whose expression is correlated and heteroskedastic. In an independent 1,000-repetition test with \(G=1200\), \(N=24\), and intercept residual rank \(r=23\), the theoretical mean squared gain was 58.0952 and the observed mean was 59.1409, within 1.77 theoretical Monte Carlo standard errors.

## 5. Numerical audit

`independent_operator_tests.py` tested 100 independently constructed dense geometries, dimensions 2, 3, 8, 31 and 80, perpendicular loading norms 0, 0.3, 1, 8 and 20, and caps 1, 1.1, 2 and 10. Dense SVD, idempotence, score cancellation, noise covariance bounds, sensor-null biology preservation, Euclidean closest-loading construction, outside-row norm and partial-correction singular values all passed. Maximum absolute discrepancy was approximately 1.1 × 10⁻¹⁴. Tests compare dense matrices with production transformations; they do not simply repeat the production expression.

## 6. Audit of the existing validation extension

1. **Correct scope of conditional OLS.** With a fixed operator trained independently of the test samples, independent identically distributed Gaussian sample errors and a correct sample design, each transformed gene can be tested by ordinary OLS. Cross-gene covariance does not invalidate the marginal t test. Biological false discoveries instead compare to original \(b\); when residual technical group differences remain, these are not null hypotheses for the transformed group means. Label the simulated metric biological FDP or recovery error, not demonstrated biological FDR control.
2. **BH dependence assumptions.** Signed latent gene loadings and the oblique transform create general dependence. Classical BH validity is not guaranteed for arbitrary dependence. The simulated mean FDP is empirical for the tested model and 200 repetitions, with Monte Carlo uncertainty; it is not a universal error-control theorem.
3. **Oracle nomenclature.** The root oracle uses the true ground nuisance direction, not the true target direction when coupling changes. It is an oracle-ground reference under domain shift, not a uniformly best biological estimator.
4. **Gain mismatch.** Root Gaussian learned gains average approximately 2.14–2.15, versus real score gains 8.13–15.31. The amended high-gain experiment below explicitly addresses this mismatch and is labelled an exploratory amendment.
5. **Joint bootstrap.** Resampling collapsed within-region animal columns by flight/cohort/age is appropriate for the stated per-region bootstrap. Relearning normalization and anchor and resampling shared OSD-525 HGC once address those estimated steps and reuse. Intervals remain percentile intervals conditional on the initially filtered gene universe and observational sampling assumptions. Joint normalization across ground and flight means the ground anchor is not fully independent of target expression in the real-data pipeline, even with disjoint animals. The joint bootstrap accounts for this shared computation, but independent-training claims should be restricted to the simulations or an explicitly separately normalized training pipeline.
6. **Small sample calibration.** Pooled CR1 approximate gene p values, fixed-anchor p values and label permutation assumptions do not become causal or exact biological inference because numerical code checks pass. Unknown collection confounding and metadata conflicts remain material.
7. **Degenerate residuals.** `Contrast.fit` presently sets t=0 when standard error is zero. For a nonzero exactly fitted coefficient this should be signed infinity or explicitly flagged as an exact fit. Real inputs are unlikely to trigger it, but deterministic signal benchmarks can. A non-finite-input check in `Contrast.__init__` would also make failures more transparent.

## 7. Exploratory amended high-gain stress benchmark

The amended script saves its exact configuration before executing. All 16 combinations are retained: true score gain 8 or 12; training size 24 or 96; target size 24 or 48; unchanged or flipped outside-sensor coupling. There are 100 independent repetitions per combination, 1,600 total repetitions and 6,400 method rows. Four methods are compared without target tuning: raw, genome-wide orthogonal projection, original score, and bounded score with fixed \(L=2\). The 1,200-gene Gaussian model has 31 sensors, 1,120 diffuse nuisance loadings, 60 biological effects of magnitude 0.8 entirely outside the sensor, gene noise SD 0.65–1.15, and five signed correlated-noise factors. True biology is exactly sensor-orthogonal; it need not be orthogonal to the genome-wide nuisance direction. Flipping reverses every outside-sensor nuisance loading while preserving sensor loadings.

Learned unbounded gains ranged from 6.712 to 12.573. Bounded gain never exceeded 2 beyond floating-point roundoff (maximum 2.0000000000000004). Both score operators preserved the known sensor-orthogonal biological vector exactly in every repetition. This operator preservation did **not** imply accurate recovery of the full observed biological contrast in the presence of nuisance.

Across the eight unchanged-coupling settings, equally weighted descriptive mean effect RMSE was raw 1.648, projection 0.388, score 0.401, and bounded score 1.377. Corresponding mean biological FDP was 0.956, 0.219, 0.215 and 0.953. Bounded correction left considerable diffuse nuisance because the true required oblique gain was far above two. Across the eight flipped-coupling settings, descriptive mean RMSE was raw 1.633, projection 0.502, score 3.093, and bounded score 1.889. Mean biological FDP was 0.957, 0.680, 0.948 and 0.954. Bounding substantially reduced the original score's worst domain-shift amplification but did not produce reliable biological recovery. Every individual configuration and Monte Carlo standard error is in the saved summary, and all repetitions are available.

The four-panel PNG at 600 dpi and editable SVG show both coupling regimes and all sample sizes. Error bars are ±1.96 Monte Carlo standard errors of the mean over 100 independent repetitions, not confidence intervals for real biological effects. Figure labels and legend were visually inspected. The figure's conclusion is a tradeoff between a guaranteed amplification limit and incomplete nuisance removal, not superiority of bounded score.

## Actionable manuscript conclusions

The defensible stronger contribution is a transparent diagnosis of when ground-derived expression adjustments fail, a proven geometric amplification safeguard, and an explicit distinction between sensor-null signal preservation and biological recovery. The available simulations and real data do not justify universal recovery, corrected biological FDR, or restored strong causal claims about spaceflight-induced inflammation. The high-gain failures are a useful scientific result and should be included, not filtered out.

A ground-only operational rule can predeclare a gain cap, report the learned gain and amount of shrinkage, assess held-out ground stability, and retain raw effects with bounded/projection sensitivity results. It cannot choose the scientifically correct biological adjustment from ground associations alone. External matched collection controls or valid negative controls would be needed to validate that choice.

## Files and reproduction

- `independent_operator_tests.py` → `independent_operator_test_results.json`.
- `fixed_sensor_noise_inflation.py` → `fixed_sensor_noise_inflation_results.json`.
- `high_gain_benchmark.py` → exact config, checks, repetition CSV, summary CSV, panel PNG and SVG.

Run each with the repository's `research/brain_integrative_2026/.venv/bin/python` and `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1`. This extension does not modify original root benchmark tables.
