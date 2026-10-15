# Frozen multiway composition transport: scientific report

Author: Isarar Siddique  
Affiliation:  
Email:  
Analysis date: 2026-10-15

## Main result

The measured multiway RR23 FLT–GC composition profile does not show aggregate positive-direction transfer to independently launched RRRM2 cerebellum. The equal-cluster held-out effect is **-0.277 projected percentage points**, with exact one-sided **p=0.73625** and absolute two-sided p=0.521875 across all 11,200 stratified assignments. The 2,000-draw fixed-direction animal-bootstrap 95% interval is [-0.984, 0.371] pp. Independently resampling both cohorts and re-estimating RR23 directions gives a wider interval [-1.061, 0.739] pp. Neither interval excludes zero; the negative point estimate does not establish a reversed response.

Only **281/923 clusters (30.44%)** are measurable in both cohorts. Of these, 46.98% have positive learned-direction effects (conditional animal-bootstrap 95% interval 41.99% to 51.96%). Mean native cosine alignment is -0.0436, reported descriptively without cluster-independence p values.

## Design and estimand

The retrospective plan was recorded before inspecting these RRRM2 outcomes; previous binary transport and RR23 multiway results were already known. All 923 all26-defined clusters and all 3,376 alternative occurrences were retained. The analysis did not use RR23 significance, favorable effect size, positive held-out direction, gene annotation or Tcaf1 status to select the testing family.

Native composition divides each frozen alternative count by the sum over that cluster's fixed alternatives. Coverage and dominance were evaluated in every animal: ≥20 selected-alternative reads and ≥90% dominance among all retained annotated canonical alternatives sharing the boundary. RR23 gates use all 26 animals; RRRM2 uses all 19 included animals. Coverage fails in 604 clusters and dominance in 312, with 274 failing both. The disjoint failures are 330 coverage only, 38 dominance only and 274 both. Inadequately measured clusters receive p=1 and no interpretable held-out effect; an absent single alternative can have a valid zero fraction when its complete cluster passes.

For each cluster, the RR23 9-FLT minus 9-GC native mean vector is normalized by its Euclidean norm. No training direction is estimated from RRRM2 labels. The projected held-out score subtracts the RR23 GC mean composition and takes a dot product with that frozen unit direction. The global score averages all common clusters equally. None of the 923 training directions met the predeclared near-zero threshold. Unit normalization avoids weighting a cluster by its RR23 effect magnitude but can leave uncertain directions for small training effects; the separate two-cohort bootstrap exposes this uncertainty rather than silently selecting stable or larger effects.

The prior RRRM2 inclusion rule was preserved: the entire ISS_T×29-week stratum, including its corresponding controls, remains excluded because its four FLT records have conflicting location metadata. Included strata are ISS_T_12 (5 FLT/3 GC), LAR_12 (2/3), and LAR_29 (3/3). Within-stratum contrasts receive fixed animal-count weights 8/19, 5/19 and 6/19. Whole-animal allocation preserves every within-animal alternative/cluster dependency. Its interpretation requires conditional exchangeability within these strata; retrospective labels were not physically randomized.

## Complete local correction and exact-test resolution

All **923** directional local p values were corrected together, with failed measurements retained at p=1. There are **zero BH and zero BY discoveries** at 0.05; the smallest raw directional p is 0.001429. BH is an exploratory dependence-assumption sensitivity, while BY is the conservative step-up map under valid p values. Zero local calls do not establish that every true effect is zero.

Zero Holm calls are structurally guaranteed at the declared family size and label space: the smallest attainable raw one-sided p is 1/11,200=0.000089286, whereas the first Holm step requires p≤0.05/923=0.000054171. Its minimum attainable first-step adjusted p is **923/11,200=0.0824107**. This first-step bound applies regardless of the observed effect and must not be described as biological evidence of absence. BY can still reject through its step-up rule if sufficiently many small p values occur; none did here.

The earlier full 1,512-common-feature binary directional family has the same issue: its minimum attainable first-step Holm-adjusted p is **1,512/11,200=0.135**, exceeding 0.05. This statement concerns the first sorted p value and the consequent inability to initiate any Holm rejection, not a universal minimum of every unadjusted or BH/BY p value. It does not affect the separate one-score global tests. No correction family was reduced to rescue significance.

## Measurement comparison and candidate illustration

The binary panel measured 1,512 of 5,249 RR23 measurable features (28.81%); this multiway panel measures 281 of 923 (30.44%). These are distinct observables with different selection and geometry, so the proportions are descriptive and do not show that one assay has statistically superior transport or that their effect magnitudes are directly comparable.

The historically selected Tcaf1 cluster `6:2:42663893:start` is measurable (three fixed alternatives; RRRM2 minimum selected reads 76, minimum dominance 99.48%). Its held-out directional effect is -5.389 pp (conditional 95% CI -14.237 to 4.691), one-sided p=0.859464, two-sided p=0.275714; BH/BY/Holm q values are all 1. Its native held-out contrast has cosine alignment -0.973 with RR23, a descriptive opposite alignment with wide uncertainty, not evidence of a reversed causal mechanism. All three alternatives and all 45 RR23/RRRM2 animal compositions are preserved; no alternative was chosen to make the comparison favorable.

## Sensitivities and strata

Two planned stricter gates retain the primary direction definition and use the complete stricter common panel. Neither rescues aggregate positive-direction transport; their two one-sided sensitivity p values have joint Holm values of 1.

| Measurement panel | Common clusters | Directional contrast (pp) | Exact one-sided p |
|---|---:|---:|---:|
| primary | 281 | -0.277 | 0.736250 |
| dominance95 | 180 | -0.255 | 0.703750 |
| depth50 | 91 | -0.358 | 0.782946 |

Stratum effects are descriptive, based on very small cells, and are not independent validation datasets:

| Stratum | FLT/GC | Directional contrast (pp) | Conditional bootstrap 95% interval (pp) |
|---|---:|---:|---:|
| ISS_T_12 | 5/3 | 0.058 | [-1.147, 1.124] |
| LAR_12 | 2/3 | -0.481 | [-1.794, 0.832] |
| LAR_29 | 3/3 | -0.554 | [-1.461, 0.353] |

## Interpretation and reproducibility

The additional real-data contribution is a frozen multialternative composition-direction test and an explicit audit of transportability, measurement failure and finite exact-test resolution. It provides an independent-launch boundary on the generality of RR23 observations, not a new conserved flight-injury mechanism. Mission, sex, strain, age, recovery/collection route, tissue sampling, preservation and sequencing coverage differ between RR23 and RRRM2 and cannot be separately identified as causes of failed transfer. RR23 flight and collection date are aliased. No protein, functional injury, cell-intrinsic or isolated causal flight claim follows. RR10 shares SpaceX-21 with RR23 and is not the independent-launch validation cohort.

All 10 independent numerical checks passed. A separate implementation read all 805,562 RR23 and 668,237 RRRM2 count rows in chunks and independently recomputed every frozen alternative count, composition, full-boundary total and dominance gate. It checked every RR23-only unit vector and held-out projection, directly enumerated all 11,200 global allocations, reproduced 17 row-order-selected local exact p values and all 923 BH/BY/Holm values, and verified both whole-animal bootstrap series and intervals. Maximum independent unit-direction arithmetic difference was 4.26e-14, below the declared 1e−12 tolerance. All frozen-source hashes and the dated plan remain unchanged.

The extension is self-contained: 15 frozen input files include actual count matrices, full animal metadata, all cluster memberships/features, annotation lookup and NASA publisher provenance. The inherited RRRM2 provenance documents 26 raw junction files matching publisher MD5 and complete reconstruction of the 19 retained columns. Reproduction requires no path to the source repository. Run `scripts/analyse.py`, `scripts/verify_independently.py`, then `scripts/build_outputs.py` from any directory. The figure uses 183 mm width, Arial/Helvetica at ≥8 pt, ≥1 pt lines, lowercase bold panel labels, editable vector PDF/SVG and 600 dpi PNG, with title/caption outside the canvas.

Complete outputs include `tables/all_923_cluster_directional_results.csv`, `all_923_cluster_measurement_gates.csv`, `all_3376_alternative_group_composition_summaries.csv.gz`, all animal/alternative counts and fractions, frozen training vectors, all held-out animal scores, the full 11,200-assignment null, saved 2,000-draw conditional and two-cohort bootstraps, Tcaf1's complete observable case, and independent check records. `figures/Figure_Multiway_Transport_caption.txt` is the separate manuscript caption.
