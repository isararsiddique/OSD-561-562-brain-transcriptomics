# Retrospective independent-launch validation plan

Recorded 2026-10-15 before reading new-panel RRRM2 outcomes. This is a dated retrospective analysis plan, not a prospective preregistration. Historical results on a different RR23 panel are known; they will not be substituted for these results.

Author: Isarar Siddique. Affiliation: [blank]. Email: [blank].

## Question and fixed inputs

Does the FLT–GC junction-fraction profile learned from RR23 cerebellum transfer to independently launched RRRM2/OSD-561 cerebellum? RR10 shares SpaceX-21 with RR23 and is excluded from independent-launch validation. RR23 flight status and collection day are aliased, so transfer cannot identify a causal flight effect or functional injury.

Use the frozen RR23 all26 label-invariant definitions (5,724), animal fractions and measurement gates in `research/paper2_beta_20261014`. No significant-event, large-effect, gene-annotation, or direction-positive subset will be selected. All 26 RR23 animals establish feature definitions; the FLT–GC contrast uses the 9 FLT and 9 GC animals. VIV animals remain part of the frozen definition and measurement gate.

Reconstruct the same junction pairs from animal-aggregated OSD-561 junction counts in `research/priority_osdr_study_2026/study_1/source_data`. Preserve all source files. Check metadata and the historical inclusion rule excluding ISS_T age-29 records with unresolved conflicting assignment. Retain distinct biological animals only, with technical libraries already summed. RRRM2 contrast is weighted within collection-route × age strata; each stratum receives its fraction of the 19 included animals. Verify each retained stratum contains both exposures.

## Measurement and estimand

Common primary panel requires the frozen RR23 `measured_primary_all26` gate and, in every included RRRM2 animal, at least 20 reads across the pair and pair dominance at least 0.90 among all retained annotated canonical alternatives sharing the boundary. These gates do not use exposure labels. Absence or inadequate coverage is a measurement failure, not zero biological effect. Use fraction junction_a / (junction_a + junction_b). The RR23 effect is mean(FLT) − mean(GC). The RRRM2 effect is the animal-count-weighted average of within-stratum FLT − GC means.

Report measurement transportability for every frozen definition and for the RR23 measured panel; show coverage and dominance failures separately. Describe the resulting panel and do not imply that results apply to unmeasured features.

## Primary held-out inferential test

Freeze one equal-weight learned-direction score from RR23 on the entire common measured panel: for animal i, S_i = mean_e[sign(delta_RR23,e) × (psi_i,e − mean_RR23_GC,e)]. Events with exactly zero training contrast get weight zero. Direction and centering use RR23 only; all eligible events remain in the score. There is no training on RRRM2 exposure outcomes. Test the weighted stratified FLT − GC score contrast by enumerating every label assignment preserving FLT counts within RRRM2 route × age strata. Report exact one-sided p for the prespecified positive transfer direction and exact two-sided p, including the observed assignment. This is one primary test, alpha 0.05. No gene-independence or read-level inference.

## Profile and uncertainty

Report Spearman correlation, Pearson correlation and same-direction fraction between complete common-panel event effects descriptively. Do not use correlation p values that treat events or genes as independent replicates. Use 2,000 independent animal bootstraps within RR23 exposure groups and RRRM2 route × age × exposure cells, preserving fixed stratum weights and panel, to obtain percentile 95% intervals for Spearman correlation, directional concordance, and score contrast. RR23 directions stay fixed for the primary held-out score uncertainty; bootstrap RRRM2 animal scores only for its interval. Save all bootstrap draws and seeds.

Report route × age score contrasts with bootstrap intervals descriptively, without selecting favorable strata. Exact event-level directional p values can be reported as a secondary map of the full common panel, with Holm correction over every common-panel event and directions frozen from RR23. This is not a significant RR23 subset. Null replication is a valid result and must be stated plainly.

## Prespecified sensitivity and checks

Repeat complete-panel profile and learned-direction score with pair dominance ≥0.95 in both studies and ≥50 pair reads in both studies. These are secondary measurement sensitivities; retain the primary test regardless of their results. Record exact primary p and Holm-adjust the two sensitivity one-sided p values for descriptive safeguards; do not use a favorable sensitivity to replace a failed primary test.

Independently recalculate pair counts, fractions and dominance on sampled/all matched rows; verify RR23 fractions directly from the all-animal junction counts; verify one row per animal; enumerate assignment counts combinatorially; independently recompute exact score p with a separate label-combination implementation; verify score projection from the event-effect table; recompute SHA256 checksums of inputs before and after processing. Save scripts, input manifest, complete mapping/gates, effects, animal scores, exact score null, bootstrap draws, summary JSON, independent check results and a four-panel validation figure (measurement transport, common-panel effect scatter, animal-bootstrap correlation uncertainty, held-out scores by route × age).

## Interpretation

Validation, if observed, supports transport of a measured junction-fraction association across these sampled cerebellar cohorts. Lack of transfer does not prove no biological response or invalidate every individual event; it bounds the generality of the aggregate learned profile. Neither result supports causal flight or functional injury claims. No outreach or publication is authorized by this analysis.
