# Retrospective plan: frozen multiway composition transport

Recorded 2026-10-15 before inspecting RRRM2 multiway transport outcomes. This is a dated retrospective plan, not prospective preregistration. Previous binary-junction transport results and RR23 multiway results are already known; neither will be substituted for this new analysis.

Author: Isarar Siddique. Affiliation and email: blank.

## Fixed family and independent cohort

Test all 923 frozen RR23 multiway clusters, defined without group labels using all 26 animals in the verified `extensions/multiway_splicing` analysis. Do not select clusters by RR23 significance, effect size, favorable direction, annotation, or prior candidate status. Alternative order, boundary identity, and cluster membership remain frozen. Reconstruct the same alternative counts in RRRM2/OSD-561 cerebellum from its actual verified animal-count matrix (668,237 junction rows × 19 included animals) and complete inclusion metadata, frozen under `extensions/independent_validation/inputs`.

RRRM2 has 10 FLT and 9 GC animals across route × age strata ISS_T_12 (5/3), LAR_12 (2/3), and LAR_29 (3/3). Preserve the prior exclusion of the entire disputed ISS_T_29 stratum. Each retained column must represent one distinct animal. Technical libraries are already collapsed. RR10 shares the SpaceX-21 launch context with RR23 and supplies no independent-launch claim here.

## Composition, measurement and direction normalization

For frozen cluster c with J alternatives, native composition p_cji is junction-j unique reads divided by the sum of reads for its J fixed alternatives in animal i. Its denominator is the frozen selected-alternative sum, while dominance compares that sum with ALL retained annotated canonical junction reads sharing the boundary. The cluster is measurable in a cohort only if every included animal has ≥20 selected-alternative reads and ≥90% dominance. RR23 measurement is evaluated in all 26 animals. The common primary panel requires both cohort gates, independent of exposure labels. An absent alternative may legitimately have zero fraction when the whole cluster passes; an absent or inadequate denominator is measurement failure.

Train one native-composition contrast vector d_c solely from the RR23 9 FLT minus 9 GC group means. Before any held-out outcome, fix its normalization as u_c = d_c / ||d_c||_2. This unit Euclidean direction removes the training contrast magnitude from cluster weights. No sign or magnitude is learned from RRRM2. If ||d_c||_2 ≤ 1e−12, use the zero vector, retain the cluster in the measurement family and score denominator, and assign its directional p=1. All nonzero directions should sum to approximately zero because their underlying compositions sum to one.

For held-out animal i, the fixed cluster score is z_ci = u_c · (p_ci − mean_RR23_GC,c). Centering uses RR23 only and does not change a within-stratum contrast. The cluster transfer effect is the animal-count-weighted average of RRRM2 within-stratum FLT−GC score means. Positive transfer is the prespecified direction. Native directional projections are linear fraction contrasts, not Hellinger distances or percentages of genes altered. Scaling them by 100 expresses projected fraction units in percentage-point units.

## Primary global test and complete local family

The primary global equal-cluster score S_i averages z_ci over EVERY common measurable cluster with equal cluster weight, preserving dependence within an animal. Test its weighted within-stratum FLT−GC contrast using all 11,200 exact allocations preserving FLT counts separately in each route×age stratum. Include the observed allocation; ties use absolute tolerance 1e−12. Report exact one-sided p for positive transfer and exact absolute two-sided p. This is one primary global test at alpha 0.05, and will not be replaced by a favorable local result or sensitivity.

Each common measurable cluster also receives an exact one-sided directional p using its RR23-frozen unit vector and the same 11,200 allocations. Keep all 923 frozen clusters in the local family, assigning p=1 to failed measurement/zero-direction clusters. Report BH, BY and Holm adjustment across all 923. BH is an exploratory dependence-assumption sensitivity; BY and Holm provide the conservative evidence map under valid conditional animal-level p values. No cluster, gene, junction read, technical library, bootstrap draw or allocation is treated as an independent biological replicate.

## Uncertainty and descriptive characterization

Use 2,000 whole-animal bootstrap draws (seed 2026101503) within RRRM2 exposure×route×age cells, keeping panel, RR23 directions and stratum weights fixed. Save all draws; report percentile 95% intervals for the primary global effect and every common-cluster directional effect. These conditional intervals do not include training-direction uncertainty. Also report a secondary two-cohort animal bootstrap that re-estimates RR23 native contrast directions within its FLT/GC groups, using the same fixed common panel and fixed stratum weights; label this training-uncertainty analysis explicitly and do not replace the primary exact test.

Describe the entire common panel using native effect-vector norms, cosine alignment with the fixed RR23 direction, and the distribution of directional effects; do not compute event-independence correlation p values. Report the fraction of positive directional effects descriptively, with animal-bootstrap uncertainty. Strongly correlated/overlapping clusters remain intact through all animal allocations and resampling.

Tcaf1 is a historically selected RR23 illustration. Report its frozen cluster identity, selected alternative counts/dominance and compositions ONLY if measurable in RRRM2; otherwise state the failed gate explicitly. Its held-out inference remains within the complete 923-cluster correction family. Do not force a positive interpretation, choose a favorable alternative, or substitute a renamed candidate if Tcaf1 is unobservable.

## Planned measurement sensitivity, provenance and verification

Repeat the equal-cluster score on two stricter label-invariant common panels: ≥95% dominance and ≥50 selected-alternative reads in both cohorts. Report two sensitivity one-sided p values and their joint Holm adjustment; retain the primary result regardless. Compare measurement transport against the prior binary panel descriptively, not as a significance contest between different observables.

Freeze copies of required RR23 cluster definitions/features, count/meta sources, exact annotation lookup and RRRM2 count/meta/provenance under this extension so the full analysis reproduces independently of the source repository. Use script-relative paths; preserve all originals and SHA256 hashes. Independently reconstruct fixed-alternative compositions and full-boundary dominance; independently recalculate learned unit vectors/projections; enumerate score allocations with a separate direct-combination implementation; check representative local exact p values selected by row order; independently reproduce all BH/BY/Holm values and all measurement decisions. Save complete local results (923 rows), animal/alternative compositions and scores, exact null, bootstrap draws, design/gates, software versions, numerical checks and report.

The publication figure uses 183 mm width, fonts ≥7 pt, lowercase a–d labels, no title/footer inside the canvas, vector PDF/SVG and 600 dpi PNG. Show common-panel measurement transport, complete directional-effect summary/corrected evidence, bootstrap global uncertainty and held-out animal scores or all measurable Tcaf1 alternatives if appropriate. Captions stay in a separate file.

Before execution, the parent verified journal guidance and tightened presentation requirements to Arial/Helvetica, fonts ≥8 pt at final size where feasible, lines ≥1 pt, bold lowercase `a)`, `b)` panel labels, and a white background. This changes presentation only; the statistical plan above is unchanged.

## Interpretation

A positive result would support transport of a measured composition association across these sampled cohorts under conditional label exchangeability; a null result would limit such transport. Sex, strain, age, collection/recovery route, tissue preservation and coverage differ between cohorts. RR23 flight and collection day are aliased. No isolated causal flight effect, functional injury, protein mechanism, cell-intrinsic mechanism or universal conserved signature is identified. The local and global families are retrospective analyses of existing data; no broad study-wide prospective error guarantee is claimed. No outreach or publication action is part of this analysis.
