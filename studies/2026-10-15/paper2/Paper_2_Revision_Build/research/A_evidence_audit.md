# Paper 2 — evidence audit

**Target:** npj Microgravity. **Manuscript audited:** `/Users/zain/Downloads/OSD-561-562-brain-transcriptomics-main/DOWNLOAD_PACKAGES_2026/Paper_2_Complete_Submission_Package/manuscript/PAPER_2_MANUSCRIPT.tex` (320 lines, 15 Oct 2026 build).
**Audit date:** 15 October 2026. **Scripts:** `A_audit.py` (entry point) + `A_claims.py` (claim table). **Claim-by-claim table:** `A_claims_check.csv`.

## Headline

| | |
|---|---|
| Numeric claims extracted and checked | **350** (Abstract, Results incl. Table 1 and all eight figure captions, Discussion) |
| MATCH | **349** |
| MISMATCH | **1** (a last-place rounding slip, below) |
| NOT-FOUND | **0** |

Every numeric claim in the Abstract, Results and Discussion resolves to a value in the package's own result files. Nine headline quantities were additionally **recomputed from the raw/complete tables** rather than read from a summary; all nine reproduce. The manuscript's numbers are sound. The weaknesses are in wording, not arithmetic — section headings and Abstract sentences that assert more than the verified numbers support.

Nothing in the package was modified. No pipeline step was rerun; only existing outputs were read (stdlib `csv`/`json`/`gzip`, no pandas — it is not installed, and host memory was tight).

## The one mismatch

| Section | Sentence | In text | Source value | Source |
|---|---|---|---|---|
| Results/robustness | "versus median rho=−0.4133 across disjoint four/five GC partitions" | −0.4133 | **−0.41324882** | `results/profile_dependence.json` → `disjoint_rho_median` |

−0.41324882 rounds to −0.4132 at four decimals, not −0.4133. Substantively irrelevant, but it is a visible arithmetic error in a paper whose contribution is precision about thresholds. **Fix:** write `rho=-0.4132`.

## Independent recomputations (part 2 of `A_audit.py`)

All reproduced from the complete families, not from the summary JSONs that the manuscript cites.

| Recomputed quantity | Manuscript | Recomputed | Source read |
|---|---|---|---|
| BY harmonic factor, m=2769 | 8.50364 | **8.50363775** | `sum 1/k, k=1..2769` |
| Exact p floor | 2/48,620 = 1/24,310 | **4.1135335253e-05**, both forms identical | arithmetic |
| q for a 20-test floor block | 0.04843 | **0.04842981** | `m·H_m·p_min/20` |
| q for a 19-test floor block | 0.05098 | **0.05097875** | `m·H_m·p_min/19` |
| Smallest floor block crossing 0.05 | 20 | **20** | solved for k |
| Table 1 Holm p, all 6 rows | see Table 1 | **all 6 reproduce to 1e-12** | Holm step-down from the stored exact p |
| Multialternative magnitude = sqrt(T/2) | 0.08313 / 0.06615 / 0.06551 | **identical** | `global_exact_results.csv` statistic column |
| Gene BY positives | 47 FLT–GC, 0 / 0 vivarium | **47 / 0 / 0**, min q 0.048256988, all 47 on the p floor | BY recomputed over all 5,949 tests in `gene_omnibus_complete.csv.gz` |
| 47-gene RMS range | 3.42–22.77 pp | **3.42–22.77 pp** | `arbitrary_dependence_corrected_genes.csv` |
| Cluster BY positives | 40 + 3 = 43, 41 unique | **40 + 3 = 43, 41 unique clusters**, 30 tests on the p floor | BY recomputed over all 2,769 tests in `local_exact_results_complete.csv.gz` |
| Split-pair positive alignments | 354 of 430 possible | **354 measured (325 + 29), 354 positive**; medians 0.996434 / 0.993302 | `all_complementary_half_concordance.csv.gz` |
| B-half BY calls | 20, 20, 23, 21, 20; total 104; 5 of 10 | **[0,20,0,20,0,0,23,0,21,20]; total 104; 5 of 10** | `complement_B_local_exact_full_families.csv.gz` |
| A-half floor ties all < 20 | yes | **[18,13,15,13,17,19,15,17,15,18], max 19** | `BY_resolution_envelopes.csv` |
| RR10 measurable clusters | 702/923 (76.06%) | **702/923, 76.0563%** | `RR10_all923_measurement_gates.csv` |
| RRRM2 measurable clusters | 281/923 (30.44%) | **281/923, 30.4442%** | `RRRM2_all923_measurement_gates.csv` |
| Tcaf1 long-intron fractions | 78.13 / 61.72 / 65.87% | **78.1322 / 61.7234 / 65.8661%** (per-animal means, n=9/9/8) | `Tcaf1_all26_alternative_read_support.csv`, 23,000-bp alternative |
| Tcaf1 contrasts | 16.41 and 12.27 pp | **16.4088 and 12.2661 pp** | same |
| DM overdispersion | median 0.01724, IQR 0.01024–0.02753 | **0.01724208, 0.01024055–0.02752504** | `all_model_fit_diagnostics.csv.gz`, 923 full-depth fits |
| Alpha-bound hits | 1 full-depth, 44 A-half | **1 and 44** | same |
| Holm floors, external cohorts | 923/11,200 = 0.08241; 1,512/11,200 = 0.135 | **0.08241071 and 0.13500000** | arithmetic |
| Animals / libraries / MD5 | 26 / 29 / 29 matched | **26 / 29 / 29 `official_md5_match=True`** | `inputs/junction_library_manifest.csv` |
| Depth and RIN ranges | 16.25–45.09 M, 8.8–9.5 | **16.25–45.09 M, 8.8–9.5** | `results/animal_qc.csv` |
| Three-cohort animal total | 55 | **26 + 19 + 10 = 55** | `animal_qc.csv` + `comparator_transport/summary.json` |

One recomputation did **not** land on the stated integer: see over-claim (h) below (`three A-half fits had discrepant starting-solution likelihoods` — 3 at a >0.01 log-likelihood-gap tolerance, 4 at >0.001; the manuscript states no tolerance).

## The 7 strongest defensible findings, ranked for a Q1 reviewer

### 1. A discrete multiplicity boundary that makes "the discovery vanished" and "the effect vanished" provably different things
The decisive arithmetic: in the 2,769-test family the attainable minimum two-sided p is 1/24,310, the BY harmonic penalty is **8.50363775**, so a tied block of **20** floor-level tests gives **q = 0.04842981** and **19** gives **q = 0.05097875**. All ten original A halves carry **≤19** floor ties (max 19) and yield **zero** BY calls; the five B halves that carry 20–23 ties yield **20, 20, 23, 21, 20** calls. Meanwhile **all 354** jointly measurable candidate split-pair instances keep a positive effect dot product, median cosine **0.9964** (FLT–GC) and **0.9933** (FLT–VIV).
*Files:* `extensions/effect_stability/results/BY_resolution_envelopes.csv`, `effect_stability_summary.csv`, `complement_B_local_exact_full_families.csv.gz`.
*Why it is the strongest:* it is exact, not asymptotic; it is fully recomputable from two integers; and it is the only result here that generalises beyond this dataset.

### 2. The aggregate three-arm association survives four analytically distinct tests, not four re-runs of one test
Holm-adjusted p across separate families: exact composition **1.2341e-04 / 1.3163e-03 / 2.1802e-03**; Dirichlet–multinomial score **1.2341e-04 / 3.7022e-03 / 3.7022e-03**; studentized centroid **1.2341e-04 / 1.6454e-03 / 2.9206e-03**; and **60 of 60** read-half global comparisons Holm-positive within their own three-test families.
*Files:* `extensions/multiway_splicing/tables/global_exact_results.csv`, `count_model_benchmark/results/all_global_results.csv`, `multiway_robustness/results/dispersion_and_studentization_diagnostics.csv`, `effect_stability/results/summary.json`.

### 3. Habitat and vivarium controls are not interchangeable — a positive result on the control–control contrast
GC–VIV is itself significant on the same label-invariant panel: **4.0746 pp**, exact p **0.00222131**, Holm **0.00279720** (binary); **0.06551**, exact p **0.00218017**, Holm **0.00218017** (multialternative).
*Files:* `results/global_contrasts.csv`, `extensions/multiway_splicing/tables/global_exact_results.csv`.
*Why it matters:* it is a directly actionable warning for every future spaceflight reanalysis that treats one ground arm as "the" control, and it is a positive finding rather than a caveat.

### 4. A pre-frozen, comparator-symmetric transport test that fails cleanly — and cannot be blamed on measurement overlap
Both RR23 training directions were frozen before testing, and neither produced corrected transfer in either cohort: effects **−0.2770 / +0.1824 pp** (RRRM2) and **+0.4358 / +0.0596 pp** (RR10), all four conditional 95% intervals spanning zero, all four Holm-adjusted p **= 1**; the common-270-cluster sensitivity minimum Holm p **0.90476**. The "too little overlap" explanation is closed off: RR10 measures **702/923 (76.06%)** clusters against RRRM2's **281/923 (30.44%)**, and still shows nothing.
*Files:* `extensions/comparator_transport/results/all_comparator_global_transport.csv`, `summary.json`, `RR10_all923_measurement_gates.csv`.

### 5. The local discovery list is model- and gate-contingent while the measured effect is not
The DM count score agrees with the composition statistic on evidence ranking (Spearman of −log p: **0.9435 / 0.9207 / 0.9208**) yet returns only **20** BY-positive tests across **19** clusters, **16** shared with the original 43; the multinomial-limit score returns **23**. Under the 95% dominance gate the floor ties fall from 20 to **19**, minimum q rises **0.04842981 → 0.05097875**, and the positive count goes to **0** — with **identical retained raw p values**.
*Files:* `extensions/count_model_benchmark/results/model_concordance_summary.csv`, `exact_resolution_and_BY_threshold_audit.csv`, `summary.json`.

### 6. The causal limit is stated as exact linear algebra, and the bias bound that destroys the paper's own margin is reported
Design rank **3** for **4** columns, smallest singular value **4.32e-16**, exact null vector **(0,1,0,−1)**. Plug-in tipping multipliers **1.269** (multiway FLT–GC) and **1.010** (FLT–VIV); once the joint bootstrap sampling allowance is included the positive lower bound survives only below **0.094** times the GC–VIV benchmark, and FLT–VIV has a zero lower bound at zero assumed bias. Simultaneous RIN + log2-depth adjustment moves the flight-comparison norms by **+0.291%** and **+0.707%** (cosines 0.99859, 0.99130).
*Files:* `extensions/confounding_sensitivity/results/summary.json`, `identification_nullspace.csv`, `bias_tipping_points.csv`, `QC_adjustment_sensitivity.csv`.

### 7. Tcaf1 resolved to exact transcript membership, with the annotation gap disclosed
Three measured minus-strand introns sharing acceptor-side start 42,663,893, spans **5,974 / 22,570 / 23,000 bp**; the three protein-coding matches (Tcaf1-201/-202/-203) share **2,766** CDS bases, identical phase structure and **CCDS20076**, with the junction-flanking exonic bases in the 5′ UTR; the 22,570-bp alternative matches a retained-intron transcript with no CDS. The companion intron is supported in all **26** animals at **20–169** unique reads. Counter-evidence is reported in the same breath: **1,199 of 3,248** multiway junctions have no exact GTF intron match.
*Files:* `extensions/transcript_context/tables/Tcaf1_measured_junction_transcripts.csv`, `Tcaf1_transcript_CDS_UTR_summary.csv`, `summary.json`.

## Weakest / over-claimed sentences, with proposed rewrites

Each quoted sentence is numerically correct; the problem is what the wording implies beyond it.

**(a) Abstract — "20 complementary read halves" reads as 20 datasets.**
> "Across 20 complementary read halves, local discoveries ranged from zero to 23 despite concordant candidate directions, exposing a discrete multiple-testing boundary."

The 20 halves are ten disjoint partitions of the *same* reads from the *same* 26 animals, and "concordant candidate directions" holds for the post-hoc 43-test subset; across all 923 clusters the paired positive-alignment fractions are only **88.85%** and **83.87%** (`effect_stability_summary.csv`).
**Rewrite:** "Splitting the same libraries into ten disjoint read-half pairs left local BY calls ranging from zero to 23 across the 20 halves while the selected candidates' effect directions stayed positive in every measurable pair, exposing a discrete multiple-testing boundary; across all 923 clusters paired alignment was lower (88.9% and 83.9%)."

**(b) Abstract — the 47 genes and 41 clusters are threshold-adjacent, which the Abstract does not say.**
> "Joint dependence-robust correction identified 47 gene associations and 41 multialternative clusters."

All 47 genes sit on the exact p floor and share the identical **q = 0.04825699** — 0.0017 below 0.05 — and the whole set is one floor tie away from vanishing (the DM family demonstrates exactly that: 20 ties → q 0.04843, 19 → 0.05098).
**Rewrite:** "Joint dependence-robust correction identified 47 gene associations and 41 multialternative clusters, all at the exact-allocation p floor and within 0.002 of the 0.05 threshold."

**(c) Abstract — "reproducible" is the wrong word for read partitions of one library set.**
> "These results establish reproducible within-cohort measurements and technical sensitivity boundaries..."

The package's own record says otherwise: `effect_stability/results/summary.json` → `new_biological_replication: false`, `counts_are_disjoint_halves_of_same_animals: true`.
**Rewrite:** "These results establish technically repeatable within-cohort measurements and explicit sensitivity boundaries..."

**(d) Results heading — an untested comparative claim.**
> "Aggregate robustness exceeds local discovery stability"

Nothing in the section compares the two on a common scale; they are a Holm decision count and a BY decision count in different families.
**Rewrite:** "Aggregate associations persist where local discovery decisions do not"

**(e) Abstract — Tcaf1 presented without its three failures.**
> "Tcaf1 junctions mapped to alternative first-exon and untranslated architecture, with coding structure shared across three transcripts."

Tcaf1 was selected after the corrected results were known, and its second-comparator result fails three separate checks: 95% dominance (q **0.05381**), the DM score (q **0.54639**) and the multinomial-limit score (q **0.91089**); both external projections are negative (**−5.389** pp RRRM2, **−0.337** pp RR10, Holm = 1).
**Rewrite:** "A post hoc coordinate-defined candidate, Tcaf1, mapped to alternative first-exon and untranslated architecture with coding structure shared across three transcripts; its flight–ground association was retained under count-model and dominance sensitivity, while its flight–vivarium decision and its external projections were not."

**(f) Results — "strong" is an editorial word doing work the numbers do not.**
> "These disclosed numerical limitations do not alter the strong aggregate agreement, but the model sensitivity further restricts local and two-comparator Tcaf1 claims."

**Rewrite:** "These numerical limitations do not change the aggregate agreement (evidence-rank Spearman 0.92–0.94), but the model sensitivity further restricts local and two-comparator Tcaf1 claims."

**(g) Results heading — overstates one linear model.**
> "Measured technical covariates do not explain aggregate profiles"

The point estimates move ≤0.707%, but the bootstrap ratio intervals are an order of magnitude wider (multiway FLT–VIV: **0.9513–1.1909**, `QC_adjustment_sensitivity.csv`), and only two covariates were measured.
**Rewrite:** "Adjusting for the two measured quality covariates leaves the aggregate profiles nearly unchanged"

**(h) Results — an integer that depends on an unstated tolerance.**
> "three A-half fits had discrepant starting-solution likelihoods"

Recomputed from `all_model_fit_diagnostics.csv.gz`: **3** fits exceed a 0.01 log-likelihood gap (2.749, 0.580, 0.304), but **4** exceed 0.001 and **23** exceed 1e-4. The claim is defensible only with a stated threshold.
**Rewrite:** "three A-half fits had starting-solution log-likelihood gaps above 0.01 (maximum 2.75)".

**(i) Results — rounding slip.** `rho=-0.4133` → **`rho=-0.4132`** (source −0.41324882).

**(j) Results heading — only one of the two external cohorts is an independent mission.**
> "Frozen profiles show limited independent-mission transport"

RR10 shares the RR23 SpaceX-21 launch (`comparator_transport/summary.json` → `RR10_shares_discovery_launch: true`); only RRRM2 is independently launched. The body says this; the heading does not.
**Rewrite:** "Frozen profiles show limited transport to two external cohorts, one independently launched"

**(k) Abstract — no study-wide error control, and the Abstract implies there is.**
"All three aggregate comparisons survived exact and count-score analyses" is true, but every Holm family in the paper is separate (per panel, per model, per half), as `PUBLICATION_READINESS.md` concedes ("family-specific corrections do not provide one study-wide confirmatory error guarantee"). One clause in the Abstract should say so, because a reviewer will otherwise find it themselves and treat it as concealed.

## Caveats on this audit

- Four claims are verified only against a *stated* value rather than a derivable one, and should be read as provenance checks, not independent confirmations: "spanning two launches" (a boolean flag), "a minimum of 20 raw gene reads per animal" (a gate recorded in the normalization string), "Across all 18 FLT/GC single-animal omissions" (a scope string), and the age-at-launch range "16–17" (a metadata cell). They are marked by their `source_field` text in `A_claims_check.csv`.
- Matching tolerance: a claim passes when the source value rounds to the claim's own stated precision, with a sub-half-last-place allowance for decimal values and no allowance at all for integer claims. That is what surfaced (i) and would have surfaced any off-by-one count.
- The audit checks the Abstract, Results (including Table 1 and all eight figure captions) and Discussion. Methods, Supplementary Information and the Supplementary Data archives were **not** audited.
- Claims were checked against the package's own result files. Agreement therefore establishes internal consistency and correct arithmetic — not that the upstream pipeline computed the right quantity from the deposited reads.
