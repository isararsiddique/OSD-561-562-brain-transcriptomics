# Effect stability is distinct from thresholded discovery stability

This retrospective diagnostic uses all923 frozen multialternative clusters and all2769 tests. The ten original read-thinning outcomes were known when this follow-up was planned. Every corresponding complementary half was analyzed without outcome selection; no new favorable thinning seed was chosen.

## Concrete result

The original ten A halves still have zero local BY discoveries. Among the complementary B halves, five have20,20,23,21 and20 calls and five have zero:104 contrast–cluster instances, not104 unique independent discoveries. All60 global comparisons across20halves retain their separate three-contrast Holm corrections. Thus the earlier claim about the original ten thinnings remains true, but a universal claim that no half-depth partition can retain local calls is unsupported.

The original candidate effects remain highly repeatable within these data. Among jointly measurable candidate split pairs, all354 retain positive alignment:325FLT–GC and29FLT–VIV instances. Median pairwise cosines are0.99643 and0.99330. All751 measurable original-to-half candidate instances and all771 single-animal-omission candidate instances retain positive projections. The whole923-cluster panel includes weak effects and lower concordance; all outcomes are preserved, not just selected candidates.

## Meaning and limits

This separates three questions: whether a measurement is retained, whether its observed effect direction/magnitude changes, and whether a multiplicity-corrected decision crosses0.05. None of those questions is independent biological replication. Disjoint read halves share the same animals; conditional on the original count their noises are negatively correlated. Repeated omissions also share animals. The original candidate subset was selected for strong effects and its greater technical repeatability is subject to that selection. No p value or effective biological sample size is assigned to these instance-level distributions. The result supports technical repeatability of candidate point estimates, not a stable externally validated biomarker panel.

## Mathematics

For cluster c, let d_c be the original native mean-composition contrast and h_c its perturbed counterpart. Cosine is d_c'h_c/(||d_c||||h_c||). Relative projection is d_c'h_c/||d_c||²;100% retains the original directed magnitude,0% loses it, and a negative value reverses direction. Native total-variation magnitude is50sum_j|h_cj| percentage points. All metrics use frozen alternatives. Direction and magnitude are evaluated separately from the original square-root exact test.

For ordered raw p values, BY rejects iff max_i [i alpha/(mH_m p_(i))]>=1. With m=2769,H_m=8.50363775 and the minimum exact p=1/24310, the floor-block adjustment is0.96859617/k. Twenty tied minimum-p tests giveq=.04842981; nineteen giveq=.05097875. Every positive B-half decision here is supported by20–23floor-level tests, while A halves have13–19. The complete rank envelopes are supplied. This is not permission to discard hypotheses or relax correction. It describes a discrete decision boundary, not a proof that all nonsignificant candidates are true or a quantitative estimate of power against an unknown truth.

## Full summaries

| analysis | contrast | subset | instances_total | instances_measured | positive_alignment_fraction | median_cosine | cosine_q05 | cosine_q95 | median_relative_projection | projection_q05 | projection_q95 | median_TV_pp |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| original_to_half | FLT_vs_GC | all923 | 18460 | 11816 | 0.973172 | 0.98265 | 0.40323 | 0.999848 | 0.9994 | 0.292266 | 1.69174 | 5.593068 |
| original_to_half | FLT_vs_GC | original_BY_tests | 800 | 692 | 1.0 | 0.999116 | 0.988637 | 0.999987 | 0.999975 | 0.825773 | 1.161069 | 17.133047 |
| original_to_half | FLT_vs_VIV | all923 | 18460 | 11816 | 0.960308 | 0.967657 | 0.1584 | 0.999678 | 0.996536 | 0.116767 | 1.870122 | 4.326505 |
| original_to_half | FLT_vs_VIV | original_BY_tests | 60 | 59 | 1.0 | 0.998263 | 0.989694 | 0.999949 | 0.993576 | 0.739977 | 1.256575 | 15.055911 |
| original_to_half | GC_vs_VIV | all923 | 18460 | 11816 | 0.96327 | 0.971929 | 0.218214 | 0.999655 | 0.998242 | 0.152012 | 1.848745 | 4.627012 |
| half_to_half | FLT_vs_GC | all923 | 9230 | 5337 | 0.888514 | 0.921859 | -0.586982 | 0.99934 | 0.826614 | -0.465285 | 1.878625 | 5.502282 |
| half_to_half | FLT_vs_GC | original_BY_tests | 400 | 325 | 1.0 | 0.996434 | 0.95789 | 0.999952 | 0.990544 | 0.697972 | 1.33473 | 17.280512 |
| half_to_half | FLT_vs_VIV | all923 | 9230 | 5337 | 0.838673 | 0.835286 | -0.732606 | 0.998477 | 0.69327 | -0.787505 | 2.043877 | 4.249821 |
| half_to_half | FLT_vs_VIV | original_BY_tests | 30 | 29 | 1.0 | 0.993302 | 0.964656 | 0.999582 | 0.949193 | 0.542155 | 1.270336 | 14.536017 |
| half_to_half | GC_vs_VIV | all923 | 9230 | 5337 | 0.852352 | 0.862664 | -0.714649 | 0.998296 | 0.724243 | -0.730823 | 2.041093 | 4.530291 |
| leave_one_animal_out | FLT_vs_GC | all923 | 16614 | 16614 | 0.994703 | 0.996762 | 0.881588 | 0.999974 | 1.000272 | 0.689644 | 1.303557 | 5.750245 |
| leave_one_animal_out | FLT_vs_GC | original_BY_tests | 720 | 720 | 1.0 | 0.999848 | 0.996777 | 0.999999 | 1.003341 | 0.92912 | 1.064339 | 17.026751 |
| leave_one_animal_out | FLT_vs_VIV | all923 | 15691 | 15691 | 0.990632 | 0.993039 | 0.778542 | 0.999929 | 1.000987 | 0.567836 | 1.433929 | 4.527502 |
| leave_one_animal_out | FLT_vs_VIV | original_BY_tests | 51 | 51 | 1.0 | 0.999868 | 0.998855 | 0.999997 | 1.010788 | 0.903683 | 1.078795 | 14.620385 |
| leave_one_animal_out | GC_vs_VIV | all923 | 15691 | 15691 | 0.989548 | 0.99389 | 0.793108 | 0.999942 | 1.000237 | 0.569895 | 1.418826 | 4.794982 |

## Verification and reproduction

Run scripts/analyse.py, scripts/verify.py and scripts/build_report_figure.py. The verifier independently checks all count partitions, selected coordinate reconstructions,405metric/correction/exact-permutation checks, complete dimensions and retained candidate identities. Original-half exact statistics were separately reproduced against the pre-existing tables. All results, ordered p curves, gates, Tcaf1 cases and negative values are retained. The analysis does not infer causal microgravity effects or resolve the exposure/date alias.
