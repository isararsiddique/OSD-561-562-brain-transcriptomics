# Novelty assessment for the strengthened correction study

15 October 2026. The defensible contribution is an experimentally grounded **application of gain-constrained sensor regression with verifiable optimization accuracy**, coupled to training-only validation and explicit biological failure tests. Ridge regression, spectral norm constraints, accelerated projected optimization, control-based unwanted-variation estimation and cross-validation are established ideas. No comprehensive first-priority claim or acceptance prediction follows from this literature check.

## What is established

Five newly retrieved Crossref records verify exact DOI/title/authorship metadata:

| Source | Relevant established idea |
| --- | --- |
| Hoerl and Kennard, [Ridge Regression: Biased Estimation for Nonorthogonal Problems](https://doi.org/10.1080/00401706.1970.10488634), 1970 | Quadratic regularization of unstable regression is established; do not name ordinary ridge a new algorithm. |
| Wold, Sjöström and Eriksson, [PLS-regression: a basic tool of chemometrics](https://doi.org/10.1016/S0169-7439(01)00155-1), 2001 | PLS prediction from high-dimensional correlated measurements is established; PCR/PLS are meaningful prediction comparators. |
| Johnson, Li and Rabinovic, [Adjusting batch effects in microarray expression data using empirical Bayes methods](https://doi.org/10.1093/biostatistics/kxj037), 2007 | ComBat uses known batch information with empirical-Bayes borrowing; it is not the same information setting as a control-only transferred loading. |
| Boyd and Vandenberghe, [Convex Optimization](https://doi.org/10.1017/CBO9780511804441), 2004 | Convex norm-constrained quadratic optimization, dual bounds and spectral constraints are standard. |
| Beck and Teboulle, [A Fast Iterative Shrinkage-Thresholding Algorithm for Linear Inverse Problems](https://doi.org/10.1137/080716542), 2009 | Accelerated proximal/projected methods are established. Cite this origin for the computational approach without claiming its original convergence theorem automatically applies to every modified restart implementation. |

The existing manuscript also correctly situates [SVA](https://doi.org/10.1371/journal.pgen.0030161), [control-gene/control-sample factor analysis](https://doi.org/10.1038/nbt.2931), [RUV negative controls](https://doi.org/10.1093/biostatistics/kxr034) and [batch-adjustment confidence concerns](https://doi.org/10.1093/biostatistics/kxv027). These papers predate the proposed application. Prior brain-specific work by [Gerstner et al. 2016](https://doi.org/10.1186/s12864-016-3065-8) already used RUV across mouse cortex/hippocampus/whole-brain sleep studies with time-of-day-matched controls. Collection vulnerability was independently established by Dachet/Marsh and explicitly discussed by Holley; the earlier journal audit preserves those primary-source comparisons.

`reference_1.json` through `reference_5.json` and `retrieval_provenance.json` are DOI metadata, not full-text access to all five works. The two saved priority searches return broad, partly unrelated records and do not establish exhaustive absence of similar previous methods. The earlier cached primary full texts support the unwanted-variation/collection context. Keep this distinction in any claim about the depth or exhaustiveness of the literature search.

## What the present construction adds

For training sensor matrix S, outside-sensor response Y and training sensor RMS scaling D, the strengthened method solves

`min_C ||Y − SC||_F²/(2n) + α||DC||_F²/2, subject to ||C||₂ ≤ sqrt(L² − 1)`.

The expression-removal operator is the block matrix `A = [[0, 0], [−Cᵀ, I]]`, with sensor coordinates listed first. It is idempotent, annihilates **individual sensor coordinates**, and has spectral gain `sqrt(1 + ||C||₂²) ≤ L`. The original rank-one ground operator removed a weighted sensor score; that is a different constraint and must not be conflated with removing every sensor gene.

This is a genuine constrained optimization problem, rather than clipping unconstrained ridge coefficients and calling that the constrained optimum. Spectral clipping is a valid projection onto the feasible set, but generally does not minimize a quadratic loss with an anisotropic Hessian. The applied advance is to combine the biological-resource audit, raw-unit gain constraint, fully nested preprocessing/model choice, optimization certificates and measurable generalization behavior in this spaceflight setting. Neither the spectral-ball projection nor FISTA is a newly invented mathematical primitive.

The computed Frank–Wolfe duality gap bounds suboptimality for the convex feasible problem. The run summary reports 16,084 converged certified fits with maximum relative gap below 1e−7. Independently implemented CVXPY/Clarabel checks agree in the retained verification cases. These are numerical/optimization guarantees; they do not identify the correct biological nuisance direction or certify validity of downstream biological discoveries.

The guarantee depends on the stated expression metric, retained gene universe and raw-unit sensor coefficients. It is not invariant to arbitrary gene scaling. With L=2 the covariance bound permits up to fourfold worst-direction variance, rather than guaranteeing contraction. Sensor-disjoint biology is preserved algebraically; biology with sensor support is changed or removed. Estimated parameter uncertainty and domain shift remain material even when the optimization problem is solved accurately.

## Evidence and appropriate strength of wording

The saved `method_benchmark/strengthened_summary.csv` reports mean held-out MSE ratios below the nuisance/intercept baseline in all seven settings for the constrained method. Original-ground mean ratio improvements are approximately 1.2–5.2%; the four external-context improvements are approximately 1.4–7.2%. Report these as **mean ratios over the defined holdouts**, not reductions in unknown biological confounding, independent-split confidence intervals, or universally better performance than all comparators. Pooled SSE ratios weight observations differently and must not silently replace the mean-ratio estimand. The matched-grid comparator amendment appropriately tests ordinary/clipped ridge under identical grids and preprocessing; preserve those results even when another method predicts better.

The independent datasets contribute 64 primary experimental units, comprising individual animals and two-animal pools. They include one new flight mission (RR-3) and three Earth experiments, not four new spaceflight missions. OSD-352 has only 3 + 3 animals and nonoverlapping flight/ground RIN ranges. The Earth contexts use radiation/exposure-associated scores; such scores can contain genuine biological responses and cannot be presumed nuisance. Leave-one-unit-out expression prediction tests algorithm behavior in new data, not removal of exposure biology. The repeated original 362 holdouts overlap and remain within-resource validation.

The count benchmark retains 1,400 independent simulated datasets across seven scenarios and two training sizes, with 14 methods. It tests count-moment-calibrated mechanisms, shifted loadings, sensor-supported biology, confounded batch structure and dispersion differences. The linear log-expression estimator is deliberately misspecified for negative-binomial observations. Keep all scenarios, failures and biological FDP/power outcomes, and distinguish biological-support FDP from a universal theorem about statistical FDR. It provides a much stronger failure analysis than a selected favorable simulation; it cannot establish real-data ground truth or universal dominance.

## Fair comparison and remaining boundaries

ComBat/known-batch OLS use target batch labels and protect a specified biological design. Ground methods use training controls without target truth, while the RUVg benchmark receives genuinely null synthetic control genes. Report the information available to each method. The benchmark separately retains naive post-ComBat and batch-aware downstream inference, which avoids attributing every inferential difference solely to the correction algorithm.

The RUVg comparator is a documented Python algebra implementation of a particular RUVSeq log-matrix procedure at k=1, with matrix correction and covariate-adjustment versions. It is not the entire family of RUV methods. The target PCA factor comparator is **not SVA**. SVA may be discussed as established literature, but the manuscript must not say it was experimentally benchmarked unless an actual documented SVA implementation was run. These scoped comparisons support transparent application-level evidence, rather than a general ranking of all unwanted-variation approaches.

Suggested framing: “We evaluate whether a reproducible collection signature supports transferable transcriptomic correction and implement a gain-constrained regression that makes the numerical tradeoff explicit.” Strong supported claims are reproducible collection association, measurable transfer failure, modest held-out prediction gains under an explicit amplification constraint, and retained failure under biological overlap/domain shift. Avoid “safe biological recovery,” “first-ever algorithm,” a validated universal neuroinflammatory signature, or a promise of Q1 acceptance.

The official npj Microgravity format evidence remains in the earlier `extensions/journal_audit_20261015/JOURNAL_AND_NOVELTY_AUDIT.md`: ≤150-word unstructured abstract, ≤15-word title, no strict main-text word limit, all Methods in the main file, funding inside Acknowledgments, appropriate reporting checklists, and legible RGB figures ≥300 dpi. The new optimization problem and numerical certificates belong in main Methods; supplementary files may hold full results and verification tables. Novelty and stronger figures do not replace author review, source-data access, a real archived code release or collection-controlled biological validation.
