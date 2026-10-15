# Retrospective effect stability and discrete-resolution audit

Frozen on 15 October 2026 before this module's new calculations. Previous results are known: the primary multialternative family has 43 significant contrast-cluster tests /41 clusters; all ten existing half-depth replicates have zero BY calls. This is a diagnostic follow-up, not preregistration or independent replication.

## Objective and frozen data

Distinguish stability of measured effect vectors from thresholded discoveries. Use all923 original clusters, all2769 contrast-cluster hypotheses, all26 RR23 animals, the original selected alternatives and the ten already frozen full-matrix binomial half-depth draws. No new candidate selection, altered correction family or favorable replicate substitution.

## Calculations

1. Recover native count proportions for every animal and frozen alternative. Compare original contrast vectors with each of the ten thinned vectors and with each complementary count half (original minus thinned matrix). Report signed projection onto original direction, relative projection magnitude, cosine, total-variation effect magnitude and measurement eligibility. Describe all clusters and the original43-test subset separately; the latter is post hoc and subject to original selection.
2. For each paired complementary split, directly compare the two independent read partitions' contrast directions and relative difference magnitude. They are disjoint subsets of the same observed reads, not independent animals; conditional on fixed observed counts their read noises are negatively correlated. No biological replication claim or inferential p value is assigned to read-split concordance.
3. Run the existing exact whole-animal family calculation for all ten complementary halves with original2769 slots (p=1 for lost measurements), retaining every positive and negative result. Existing half-draw family values are unchanged. Global comparisons use their valid frozen cluster set and Holm3 per replicate. Native vector diagnostics do not replace the primary square-root statistic.
4. Decompose the BY decision into measurement gates, raw exact p-value resolution and the full-family ordered-p envelope. BY rejects iff max_i i*0.05/(m*H_m*p_(i))>=1. Report that ratio, family-size/harmonic penalty, support rank, actual q minimum and exact allocation floor. Illustrate the theoretical number of minimum-p tests required to cross the boundary, but never reduce the family, change alpha, relabel failures or use a different correction to rescue a desired discovery.
5. Leave each contributing biological animal out once for each of three contrasts, recomputing native effect vectors for all923 clusters. These are descriptive influence diagnostics, without new p values or independent-validation claims.
6. Report Tcaf1 for both originally significant flight contrasts as a fixed candidate, with all replicates and omitted animals. No replacement of this candidate if results disappoint.

No new causal claim. Mathematical definitions, all metrics, gates and unfavorable results will be exported. Independent verification must reconstruct selected counts and several contrast metrics directly, check all complement count identities, recalculate correction arithmetic without importing analysis functions, and check dimensions/full families.
