# Paper 2: methods, command order and interpretation

This guide summarizes the current manuscript Methods and links the analysis families to their source code. All analyses are retrospective. Dated plans and amendments document what was decided when; they are not preregistrations.

## 1. Resolve source files and animals

Use RR23/OSD-525 (9 FLT, 9 GC, 8 VIV), RRRM2/OSD-561 (10 FLT, 9 GC) and RR10/OSD-563 (5 FLT, 5 GC). Sum technical junction libraries within source animal before calculating fractions. Complete animal vectors move together in every allocation or bootstrap.

RR23's 29 STAR junction-file MD5 values matched the processed-file sidecar. RR10's ten junction files and processing-information ZIP passed official MD5 checks; its malformed software-version sidecar remains an unresolved exception. Do not describe every metadata file as checksum-verified. RRRM2's source-location conflict was handled by the recorded full-stratum exclusion.

RR23 and RR10 are different animals from one launch; RRRM2 supplies the independent launch. No usable cage-level blocking factor or genotype fingerprinting was available. Animal exchangeability and metadata completeness are assumptions, not demonstrated experimental randomization.

## 2. Construct label-invariant junction choices

Retain canonical annotated STAR junctions on chromosomes 1–19, X, Y and MT, with strand 1/2, motif 1–6, maximum overhang at least eight and positive unique reads. A junction qualifies with at least five reads in at least three of all 26 RR23 animals.

Group junctions by chromosome, strand and shared start/end boundary. Binary choices use exactly two qualifying alternatives, lexical orientation, at least 20 selected reads in every animal, pooled mean fraction strictly between 0.05 and 0.95, and at least 90% dominance against all retained boundary reads in every animal. A 95% gate is a sensitivity. Construct features from all animals without group labels.

The final article reports 5,249 primary binary choices. Complete declared test families retain unmeasured definition slots with p=1 rather than reducing multiplicity after a gate failure. Earlier FLT/GC-only feature construction is retained as an archived first pass and is not the primary analysis.

## 3. Match genes by exact coordinates

GRCm39 release 112 transcript exons define genomic introns. Intersect gene assignments across junction alternatives and assign a gene only when the intersection is unique. Leave absent/ambiguous mappings unresolved. A coordinate match does not establish a complete transcript, reading-frame change, protein effect or pathway mechanism.

## 4. Test all three group comparisons

For FLT–GC, FLT–VIV and GC–VIV, use mean squared differences of group-average fractions over the common panel. Enumerate 48,620 allocations for 9-versus-9 and 24,310 for 9-versus-8 animals. Include the observed assignment and ties, with the recorded floating-point tolerances. Holm controls the three-test global family.

Event tests use absolute mean differences. The complete joint event family has 17,172 slots; separate contrast families have 5,724 each. Report BH and arbitrary-dependence BY. Correction depends on valid individual probabilities; label-invariant filtering does not establish animal exchangeability or a study-wide confirmatory error guarantee.

Exploratory control equivalence requires a predeclared fraction margin and approximate Welch two-one-sided tests. The conjunction takes the maximum of difference/equivalence probabilities and assigns p=1 to incompatible directions. A nonsignificant difference is not evidence of equivalence.

## 5. Localize at gene level and multialternative clusters

Gene statistics use all measured binary choices within each uniquely annotated gene. A gene Gram matrix allows exact squared-distance inference while moving entire animal vectors. Joint BH/BY covers all 5,949 gene tests, without selecting genes by event significance. Fixed-family 95% sensitivity retains failed measurements with p=1.

Multialternative clusters require at least three qualifying alternatives, 20 selected reads per animal, 90% boundary dominance, minimum pooled composition constraints and no dominant pooled alternative above 0.90. Transform within-animal compositions by square roots and build an equal-cluster kernel over 923 clusters. Global allocation tests use the complete vectors; all 2,769 local tests enter joint BH/BY. Stricter dominance, mapped-subset and equal-gene weighting analyses are separately reported.

The displayed scaled centroid distance is not the Hellinger distance between pooled mean compositions. Kernel identities, alternative sums, labels, feature selection and selected exact probabilities have independent verification scripts.

## 6. Compare relative abundance without claiming mechanistic independence

Verify deposited STAR unnormalized gene counts, collapse the same technical libraries, calculate label-invariant median-ratio size factors, and analyze `log2(normalized_count + 0.5)` over 1,977 pre-existing junction genes with at least 20 raw counts in every animal.

Exact abundance tests use 5,931 joint slots and a separate three-test global family. Approximate Welch equivalence margins are `log2(1.2)` and `log2(1.5)`, both reported. Junction/equivalence conjunctions use the maximum component probability and complete-family correction. Separate discoveries or an approximate equivalence result do not establish a junction-specific mechanism.

## 7. Examine technical robustness and finite discovery resolution

Use whole-animal bootstrap resampling and descriptive single-animal omissions. Complementary read halves use one binomial split of each junction/animal count, with the remainder assigned to the other half. Reused junctions share the same draw. Freeze cluster membership and all 2,769 test slots; reevaluate coverage/dominance and retain p=1 for failures.

Read halves are measurements of the same animals, not independent biological validation. Compare effect alignment, projection and total variation without treating features or halves as independent replicates. Retain all twenty halves and all failures.

With `m=2769`, BY harmonic factor `H_m` and smallest two-sided exact probability `p_min=2/48620`, the first attainable tied minimum-p block needs `ceil(m * H_m * p_min / 0.05) = 20` tests. Evaluate the full ordered probability curve as well as this diagnostic. Do not shrink the family to rescue a short discovery list.

Use the centered bootstrap error radius for nonnegative global-distance intervals. Signed local intervals, dispersion and studentized centroid diagnostics have separate roles; they do not provide simultaneous discovery guarantees under arbitrary heteroscedasticity.

## 8. Freeze weights before external transport

Match validation junctions and alternatives by exact coordinate/strand to RR23's frozen definitions. Apply measurement gates without validation exposure labels. Learn directions and reference means only from RR23, separately for FLT–GC and FLT–VIV; do not select a favorable comparator blend or overlap using external outcomes.

RRRM2's three route/age strata yield 11,200 restricted allocations over 19 animals. RR10's 5-versus-5 block yields 252. Positive one-sided transfer is primary, and two-sided values are also reported. The expanded comparator/cohort analysis has a four-test Holm global family and a complete 3,692-slot local family, with p=1 for unmeasured clusters. A common-intersection sensitivity remains separate.

The original 923-cluster RRRM2 local Holm family cannot attain its first 0.05 rejection because `1/11200 > 0.05/923`; RR10 is coarser. Zero corrected local calls therefore add no biological absence evidence. Expanded intervals use 4,000 whole-animal draws; fixed-weight and resampled-discovery-weight sensitivities are both retained.

## 9. Compare external transport with an internal benchmark

Within RR23, compute directions from five training animals per FLT/GC group and project onto the held-out four-plus-four contrast for every split. Repeat for externally measurable subsets and the full cluster panel; use a leave-one-pair-out variant and a label-permuted reference.

A dated amendment adds the same-date GC–VIV housing direction with all 8,820 corresponding splits and a matched reference. Overlapping splits reuse animals, and held-out RR23 animals share cohort properties including collection date. The benchmark is descriptive and cannot distinguish a causal flight effect from a cohort-associated signal.

## 10. Retain supplementary analyses and audit history

Count-model benchmarks, quality/confounding sensitivity, transcript architecture, graded transport, enrichment and design power examine specific alternative explanations. The bounded outcome-independent cohort search records 73 accession decisions, detailed inventories for 22 plausible cohorts and corrected date-bounded literature queries. It is not an exhaustive systematic review or proof that no other relevant dataset exists.

Enrichment uses complete declared tested backgrounds and retains the frozen service outputs. Marker extraction needs externally obtained Kozareva supplementary tables. The design calculator reports attainability and scenario-dependent power; post hoc power is not used to declare absence of an effect.

## Command order

From the repository root, restore the full archive and create its environment:

```sh
cd .local-reproduction/paper2/Paper_2_Complete_Submission_Package
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements.lock.txt
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
.venv/bin/python scripts/reproduce_analyses.py
```

The primary driver runs core analysis, verification and summaries, then gene omnibus, multiway splicing, abundance, independent validation, robustness, multiway transport, transcript context, RR10 frozen-count verification, effect stability, count benchmarks, confounding sensitivity and comparator transport. It uses frozen source material rather than rerunning the original source search/download.

Added revision analyses run from the restored `Paper_2_Revision_Build/research/` directory with that same Python environment:

| Stage | Commands |
|---|---|
| Discrete BY floor | `python D_simulation/sim_by_floor.py` |
| Graded transport | `python G_graded_transport/graded_transport.py`, then `disjoint_checks.py` and `five_views_figure.py` |
| Internal benchmark | `python K_internal_transport/internal_transport.py`, then `figure_internal_transport.py` |
| Enrichment | Run `I_enrichment/scripts/01_build_gene_sets.py` through `08_figure_S6.py` in order; steps 03/07 query online services and step 02 needs external marker tables |
| Design calculator | `python H_design_calculator/design_calculator.py`; `--figure` redraws only |

Use the explicit `.venv/bin/python` path from the primary environment or activate it before those commands. Seeds, plans, logs, independent verifiers and complete outputs remain in each stage's folder. Arial/Helvetica availability affects figure typography.

## Results and limits

The main article is supported by the primary package's `results/`, `extensions/` and `supplementary_data/` (Supplementary Data 1–15), plus revision `research/` tables and figures. Large complete event tables, count matrices and bootstrap arrays are in the full ZIP; every archived path and hash is listed in the GitHub release manifest.

Comparator dependence, technical decision instability and limited external transport are the supported findings. Collection date, housing, protocol, sex/strain and measurement compatibility remain constraints. The work does not establish an isolated spaceflight-specific splicing mechanism, independently launched RR10 replication or a validated treatment target.
