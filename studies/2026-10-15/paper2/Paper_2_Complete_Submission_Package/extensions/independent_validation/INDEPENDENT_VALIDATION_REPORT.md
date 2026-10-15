# Independent-launch validation of the RR23 all26 junction panel

Author: Isarar Siddique  
Affiliation:  
Email:  
Analysis date: 2026-10-15

## Conclusion

The complete RR23 FLT–GC junction-fraction profile did not show aggregate transfer to independently launched RRRM2/OSD-561 cerebellum under the dated retrospective plan. Among 1,512 events measurable in both cohorts, Spearman ρ was 0.0102 (2,000 animal-bootstrap percentile 95% CI -0.0689 to 0.0739); Pearson r was -0.0109. The prespecified held-out learned-direction score gave a stratified FLT–GC contrast of 0.128 percentage points (95% CI -0.137 to 0.422), exact one-sided p = 0.2525 and two-sided p = 0.5022. This supports a clear limit on cohort-level transport, rather than a generalized conserved flight signature.

This is a new calculation for the all26-defined panel. The historical different-panel RRRM2 result was not imported or substituted. No significant, large-effect, direction-positive, or annotated-gene subset of RR23 was selected. A null result does not prove absence of a biological response or absence of all event-specific replication. RR23 flight status is aliased with collection day; no causal flight or functional injury claim is identified by these data.

## Design and measurement

The RR23 feature definition used all 26 biological animals (9 FLT, 9 GC, 8 VIV). Its learned contrast used the 9 FLT and 9 GC animals. The frozen 5,724 event definitions are ordered pairs of annotated canonical splice junctions sharing a boundary. Fractions are junction_a / (junction_a + junction_b); they are local junction fractions, not complete isoform abundances. The source STAR junction counts have technical libraries summed before analysis, and every retained column maps to one distinct animal.

RRRM2 contains 26 animal metadata records. The entire ISS_T × 29-week stratum was excluded under the existing rule because its four FLT records have conflicting deposited versus source-inferred euthanasia location; the three corresponding GC animals were excluded with that unusable stratum. The included 19 animals comprise 10 FLT and 9 GC across ISS_T × 12-week (5/3), LAR × 12-week (2/3), and LAR × 29-week (3/3) strata. The effect averages within-stratum FLT–GC contrasts with fixed weights 8/19, 5/19, and 6/19. These labels retain their original cohort identifiers; they do not resolve all collection and processing differences.

Measurement required ≥20 reads across the fixed pair and ≥0.90 pair dominance among all retained annotated canonical junctions sharing its boundary, in every included animal. RR23 measurements include all 26 animals. RRRM2 boundary alternatives were independently reconstructed from its complete animal-level junction-count file.

- Frozen RR23 definitions: 5,724.
- Measured across all RR23 animals: 5,249.
- Measured in RRRM2 among the frozen definitions: 1,557.
- Common primary panel: 1,512, 28.81% of the RR23 measured panel.
- RR23 measured events failing RRRM2 coverage: 3,662; failing dominance: 1,301; failing both: 1,226. These failure counts overlap. The disjoint failure counts are 2,436 coverage only, 75 dominance only, and 1,226 both.

The shared panel is determined by outcome-independent measurement thresholds, not inferred significance. Poorly measured features were excluded from transfer estimation rather than interpreted as biologically unchanged. Requiring all animals to pass restricts conclusions to the measured panel and can disproportionately retain abundant or stable features.

## Held-out score and inferential unit

For each common event, the sign of the RR23 FLT–GC mean difference was learned once. Each held-out animal score is the mean signed difference between its fraction and the RR23 GC mean across every common event. Directions and centering use RR23 only; RRRM2 exposure labels do not train the score. Events are equally weighted and their correlation is preserved within each animal. Subtracting a common reference does not change a FLT–GC contrast.

Every one of the 11,200 assignments preserving FLT counts within route–age strata was enumerated. The observed assignment is included. 2,828 assignments were as or more positive than the observed contrast; 5,625 were as or more extreme in absolute value. The primary expected-direction test was one-sided at alpha 0.05. This is an animal-level permutation inference under conditional exchangeability within strata; labels were not physically randomized by this retrospective analysis, and residual collection/processing confounding remains possible.

Profile intervals resampled whole animals within RR23 exposure groups and RRRM2 exposure × route–age cells, keeping panel and stratum weights fixed (2,000 draws, seed 20261015). RR23 training effects were re-estimated in the profile bootstrap. The held-out score interval keeps the primary RR23 directions fixed and resamples RRRM2 animals only. No event-independence correlation p values, gene-level independent-replicate p values, read-level binomial replication, or technical-replicate inference were used. Same-direction concordance was 51.98% (animal-bootstrap 95% CI 47.55% to 53.97%).

## Secondary results

The entire common panel also received held-out exact directional event p values with directions frozen from RR23 and Holm correction across all 1,512 common events. Zero events passed Holm-adjusted p < 0.05. The smallest unadjusted directional event p was 0.000893; it is not a corrected replication claim.

All route–age score contrasts are shown descriptively; none was selected to rescue the aggregate result. Small cell sizes make these percentile intervals unstable and they should not be interpreted as independently powered stratum validation.

| Route–age stratum | FLT/GC | Score difference (pp) | Animal-bootstrap 95% CI (pp) |
|---|---:|---:|---:|
| ISS_T_12 | 5/3 | -0.037 | [-0.476, 0.454] |
| LAR_12 | 2/3 | 0.192 | [-0.365, 0.749] |
| LAR_29 | 3/3 | 0.295 | [-0.091, 0.679] |

The two planned label-invariant measurement sensitivities also do not replace the primary result. Their one-sided p values are additionally Holm-adjusted together in the saved sensitivity table.

| Panel | Common events | Spearman ρ | Score difference (pp) | Exact one-sided p |
|---|---:|---:|---:|---:|
| primary20_dominance90 | 1,512 | 0.0102 | 0.128 | 0.2525 |
| dominance95 | 1,222 | 0.0100 | 0.134 | 0.2402 |
| depth50 | 472 | 0.0509 | 0.261 | 0.1646 |

## Interpretation and limitations

Independent launch helps test generality across missions but is not biological matching. RR23 metadata describe male C57BL/6J mice aged 16–17 weeks at launch, while RRRM2 metadata describe female C57BL/6NTac mice at 12 or 29 weeks, with distinct recovery routes, tissue sampling and preservation. Any loss of transfer can combine mission, sex, strain, age, collection, processing, and measurement effects; this analysis cannot assign the cause. Conversely, a positive score could not have established a flight-specific causal effect.

Only 28.8% of the RR23 measured panel can be evaluated under the same all-animal thresholds in RRRM2. The near-zero profile correlation has an interval spanning modest positive and negative associations. The primary null therefore establishes neither exact equivalence to zero nor a universal absence of response. The scientifically supportable contribution is a rigorously bounded, retrospectively planned test of junction-profile transportability, including its measurement limitations. Any paper should present this null validation prominently and should not claim a conserved cerebellar injury mechanism.

## Verification and outputs

Independent verification passed all checks. A separate script did not import the analysis code: it read raw gzip CSV counts, summed boundary alternatives with dictionaries, reproduced all 5,724 RR23 fractions and all 5,724 RRRM2 fractions/pair counts/dominance rows, recalculated every held-out animal score and common event contrast, independently enumerated the 11,200 score assignments, verified 13 deterministic event nulls and all Holm values, and checked saved bootstrap interval quantiles. SHA256 hashes of all source inputs and the dated plan are unchanged.

A standalone copied `beta/extensions/independent_validation` layout reproduced the complete summary exactly with `--base ../..`, using only the five required beta files and frozen extension inputs. All independent checks passed in that copied layout; the result is saved in `tables/portable_layout_check.json`. The temporary test copy was removed.

The figure has four panels: measurement transport; all common-panel event effects; correlation bootstrap uncertainty; and individual held-out scores grouped by route–age. Panel D subtracts a common mean within each stratum for display only; all inferential calculations use the original fixed score. The mean bars are descriptive. The complete feature and animal tables preserve every eligible feature and animal.

Main artifacts are `figures/independent_launch_validation_four_panel.png` (300 dpi) and `.svg`; `tables/summary.json`; `tables/complete_feature_mapping_and_gates.csv.gz`; `tables/complete_common_panel_effects.csv.gz`; `tables/RRRM2_heldout_animal_scores.csv`; `tables/RRRM2_exact_score_null.csv.gz`; `tables/animal_bootstrap_draws.csv.gz`; `tables/measurement_sensitivities.csv`; `tables/input_manifest_sha256.csv`; and `tables/independent_checks.json`.

## Reproduce

Run the three scripts in this order with Python 3.11+ and numpy, pandas, scipy, and matplotlib, and `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1`. Paths below are relative to this extension. For the first two scripts, supply `--base /path/to/beta` pointing to the beta analysis/package root containing `inputs/` and `results/`:

1. `scripts/analyze_validation.py --base /path/to/beta`
2. `scripts/verify_independently.py --base /path/to/beta`
3. `scripts/build_report_and_figure.py`

Input paths and byte-level SHA256 checksums are recorded with portable beta/extension scopes in `tables/input_manifest_sha256.csv`. Frozen RRRM2 count and metadata inputs, publisher inventory and MD5 checksum file, original file URLs and hash verification, and independent retained-STAR-to-collapsed-matrix verification are saved inside `inputs/`. The raw STAR files are not needed for reproduction from the checked collapsed matrix. The analysis reads frozen source files without modifying them. Its frozen plan is `RETROSPECTIVE_ANALYSIS_PLAN.md`; prospective preregistration is not claimed.
