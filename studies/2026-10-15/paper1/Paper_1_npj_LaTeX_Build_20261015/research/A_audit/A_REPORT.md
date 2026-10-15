# Worker A: evidence audit and independent recomputation, Paper 1 (15 Oct 2026)

Source audited (read-only): `/Users/zain/Downloads/OSD-561-562-brain-transcriptomics-main/DOWNLOAD_PACKAGES_2026/Paper_1_Complete_Submission_Package/Paper1_Strengthened_Review_20261015/` (`MANUSCRIPT_COMPLETE.md`, `SUPPLEMENTARY_INFORMATION.md`, `inputs/`, `results/`, `extensions/`).
All outputs are in `/Users/zain/Downloads/OSD-561-562-brain-transcriptomics-main/DOWNLOAD_PACKAGES_2026/Paper_1_npj_LaTeX_Build_20261015/research/A_audit/`.
Unavailable this session: MCP servers `fusion`, `cloudflare-api`, `github` (none were needed).

## Headline verdicts

1. **Numbers.** I extracted 591 numbers from the Abstract, Introduction, Results, Discussion (none there apart from citations), Tables 1–6, all 13 figure legends, and the Supplementary results text and legends. 47 of them are conventions or definitions (95% levels, 0.05 thresholds, quantile labels, algebraic constants), so **544 are quantitative claims. 543 match their source file at the written precision, 1 is a rounding inconsistency, and there are 0 mismatches and 0 unsourced numbers.** 293 were tied to a source row automatically. The other 250 were resolved by explicit hand-check rules that read the named file and column (`A_claims_hand.py`).
2. **Independent recomputation from raw inputs (my own code).** 143 checks: **109 reproduce exactly** at the written precision. 3 agree only within Monte Carlo error, as expected for bootstrap/wild intervals with a different seed. **0 fail.** 23 are new analyses (marked NEW) and 8 are informational. The headline numbers I rebuilt from the bundled counts and metadata, without using package code, are:
   - **Accounting:** 129 libraries, 126 animal–region columns, 83 animals, 24 + 10 region pairs, 6 conflict animals, 30 unique ground animals.
   - **Panel effects:** all 136 primary design-contrast panel rows (17 contrasts × 8 panels), effects and exact permutation p, to 1.8×10⁻¹⁵.
   - **Operator gains:** all 15 gains (8.129–15.309).
   - **Ground holdouts:** AUC, cosine, MSE ratio and energy for all 81 + 81 core holdouts, to 4×10⁻¹⁶ per split.
   - **Intersection–union screen:** 54,989 tests, 0 BH and 0 BY survivors.
   - **Erythroid genes:** Hbb-bt −3.495 and Hbb-bs −3.009 with my own NB GLM.
   - **Circadian interaction:** −0.405 and its omission ranges.
   - **External units:** 64 = 46 + 18.
   - **Gain-2 constrained regression:** reimplemented end to end (own nested CV and own convex solver) from raw STAR counts for OSD-202, -352, -374, -611 and from raw counts for the 81 + 81 core holdouts. Every mean MSE ratio matches to 4 decimals (0.9878, 0.9484, 0.9284, 0.9857, 0.9501, 0.9719), as do the improved-unit counts.
3. **The real problems are design and wording, not arithmetic.**
   - **(a) Unnamed protocol factors.** The manuscript never names euthanasia agent, freezing method or kill location (0 hits for ketamine, xylazine, anaesthesia, isoflurane, Cryochiller, liquid nitrogen or freezing). In the deposited metadata, every RRRM-2 contrast it interprets is aliased with at least one of these (table in §3).
   - **(b) False statement about RR-23 dates.** The text says the "HGC-versus-VGC contrast includes different handling and harvest dates". Deposited metadata give the same euthanasia and dissection date (17 Jan 2021) for HGC and VGC. The date alias is actually flight (14 Jan) versus both ground groups (17 Jan), and the text never states it.
   - **(c) The "metadata conflicts" are cross-deposit contradictions.** For 5 of the 6 conflict animals, the cerebellum and hippocampus records of the same animal disagree with each other (same RFID, same euthanasia date). The OSDR study description says ISS-T mice were "sacrificed on-orbit". Only RRRM-2_FL_ISS_18 reads "On Earth" in both records.
   - **(d) Day counts come from the study description.** The 55–58 d and 32 d durations are quoted from the OSDR study description. The deposited per-animal `exposure duration` field gives 53–56 d (FLT ISS-T) and 33 d (FLT LAR).
   - **(e) Gain value.** Results says "about 2.14"; Table 5 says "≈2.15". The source mean is 2.147.
4. **Correction to the brief's working hypothesis.** In the ground-only contrast (GC ISS-T minus GC LAR, rPRG +0.889 log2), the deposited bulk-RNA-seq preservation is liquid nitrogen in both cohorts. The LAR string reads "for spatial transcriptomics: OCT on dry ice; for bulkRNAseq: snap frozen in liquid nitrogen". So for ground animals the aliased deposited factor is the **euthanasia agent** (ketamine/xylazine overdose versus gas anaesthesia), not the freezing method. Freezing (Cryochiller versus LN2) is aliased only in the flight ISS-T contrasts (FLT ISS-T versus FLT LAR; FLT versus GC within ISS-T; OSD-613).
5. **The strongest material for a Q1 reframe is a protocol-versus-flight contrast** (rPRG, exact whole-animal permutation tests):
   - **Ground protocol effect:** +0.889 log2 (p = 0.0014; bootstrap 0.576–1.214).
   - **Matched-protocol flight effects:** −0.03 to +0.24 log2, none significant. Pooled RRRM-2 within LAR +0.235 (p = 0.115); pooled within ISS-T −0.027 (p = 0.886); OSD-613 +0.019; RR-10 −0.05 to +0.20.
   - **RR-23 counter-signals, both strong:**
     - HGC − VGC: same date, euthanasia agent and preservation, differing in housing hardware. +0.641 (p = 0.0003, q = 0.012).
     - FLT − HGC: matched agent and RNAlater, but flight killed one day after splashdown and HGC three days later. −0.399 (q = 0.017); HSP −0.993 (q = 0.011).

     So the defensible claim is "the IEG/HSP panel tracks protocol and handling", not "spaceflight has no IEG effect".
6. **The 8–15 operator gain is not specific to the IEG sensor.** Random 30-gene sensors on the same anchors give larger gains: medians 23.1, 22.4 and 20.8, with Sensor31 below 97.8%, 100% and 99.6% of random sets. The gain is a p ≫ n property (about 18,000 genes, residual rank 10). This supports "cannot be removed post hoc", but the paper must not present the gain as a feature of the collection signature.

## 1. Claims check (`A_claims_check.csv`)

| Part | Numbers | MATCH | ROUNDING | MISMATCH | NO_SOURCE | Not a claim |
|---|---|---|---|---|---|---|
| Abstract | 7 | 7 | 0 | 0 | 0 | 0 |
| Introduction | 2 | 2 | 0 | 0 | 0 | 0 |
| Results | 311 | 286 | 1 | 0 | 0 | 24 |
| Discussion | 0 | – | – | – | – | – |
| Tables 1–6 | 109 | 109 | 0 | 0 | 0 | 0 |
| Figure legends | 58 | 49 | 0 | 0 | 0 | 9 |
| Supplementary results + legends | 104 | 90 | 0 | 0 | 0 | 14 |
| **Total** | **591** | **543** | **1** | **0** | **0** | **47** |

How the check works:
- **Automatic pass.** `A_claims_check.py` indexes 6.1 M numeric cells from 182 result files. It excludes the package's own manuscript-text and look-up tables, and tables over 3 MB. A number is matched only within the files expected for its section.
- **Source-row choice.** Among candidates, the script prefers the source row that holds the most numbers of the same sentence. It prefers summary tables over per-repetition or gene-level tables.
- **Low-precision numbers.** Integers and numbers with 1–2 significant digits are never accepted from an isolated coincidence.
- **Hand checks.** `A_claims_hand.py` holds 148 context rules backed by 130 evidence functions. Each reads the exact file and column and overrides the automatic match.
- **Audit of the automatic matches.** I reviewed all 87 low-precision automatic matches and a random sample of 45 others. About 45 rows had the right number in a coincidental row or file (for example a raw-count cell, or a per-repetition row instead of the summary). All were replaced by hand rules; none turned out to be a wrong value. After the tie-breaking was tightened, a second random sample of 40 automatic rows showed 2 residual wrong-row matches (both then fixed) and about 3 rows pointing at a sibling row of the right table (for example a different method row of `holdout_summary.csv`, where the cosine is identical across methods). Treat `source_file` for automatic rows as "value present in the expected table", not as a curated pointer.
- **Scope of MATCH.** MATCH means internal consistency between text and package tables. Independent truth of the headline values comes from §2.

The one ROUNDING row and the non-numeric factual problems the extractor cannot catch:

| Location | Text | Evidence | Action |
|---|---|---|---|
| Results, Known-truth para 4 | "learned gain was about 2.14" | `results/known_truth_benchmark_summary.csv`: mean score `spectral_gain_mean` over 5 scenarios = 2.147; Table 5 says ≈2.15 | Write 2.15 in both places |
| Results, Real-data attenuation para 2 | "the HGC-versus-VGC contrast includes different handling and harvest dates" | `inputs/OSD-525_metadata.csv.gz`: HGC and VGC Euthanasia Date = BSP Dissection Date = 17-Jan-2021; FLT = 14-Jan-2021 | **False as written.** Say: HGC and VGC were collected on the same day but differ in housing hardware; FLT was collected 3 days earlier, 1 day after splashdown (13 Jan, OSDR description) |
| Results, Ground collection cohorts para 2 | "ISS-T animals spent 55–58 days on orbit, while LAR animals returned after 32 days" | Matches the OSDR study description in `inputs/osd_561_public_metadata.json`. The deposited per-animal `Parameter Value[exposure duration]` = 53/54/55/56 d (FLT ISS-T) and 33 d (FLT LAR) | Cite the description, or reconcile with the per-animal field |
| Results, Sample identity para 3 | "Seven location entries from six core animals disagreed…" | `A_location_conflicts.csv`: 5/6 animals have contradictory locations between their own OSD-561 and OSD-562 records (same RFID, same date); only FL_ISS_18 is "On Earth" in both; HGC_ISS_20 (a ground animal) is "On ISS" in OSD-562 | Report as curation contradictions. The exclusion sensitivity, which removes all older FLT ISS-T animals, is then a very conservative bound |
| Intro, para 4 | "686 mouse samples from 13 organs" | Abstract in `inputs/reference_metadata/10.1038_s41467-026-68737-1.json`: "686 small RNA samples of female mice from 13 solid organs" | MATCH; consider adding "female" and "solid" |

## 2. Independent recomputation (`A_recompute.py`, `A_recompute_results.json`)

Own implementations throughout. Package scripts were read only for definitions. Run `…/p1env/bin/python A_recompute.py a b c d e f g h i j k l m` (each section can also be run alone). Total runtime is about 5 min; memory is small.

| Claim | Expected (paper) | Recomputed | Status |
|---|---|---|---|
| (a) libraries / animal–region columns / animals | 129 / 126 / 83 | 129 / 126 / 83 (RRRM-2 47 + RR-10 10 + RR-23 26) | MATCH |
| (a) Cb+Hip pairs; RR-10 pairs; OSD-613 other animals | 24; 10; 18 | 24; 10; 18 | MATCH |
| (a) genes retained per deposit (Table 1) | 19308 … 20655 | identical for all 7 | MATCH |
| (a) conflict entries / animals; older FLT ISS-T left after exclusion (Cb, Hip) | 7 / 6; 0, 1 | 7 / 6; 0, 1 | MATCH |
| (b) all 17 design contrasts × 8 panels (136 rows) | `results/panel_results.csv` | max \|Δeffect\| 1.8e-15; max \|Δp\| 1.1e-16; identical assignment counts | MATCH |
| (b) pooled GC ISS-T − LAR rPRG | 0.889, p = 0.001429 (700 assignments) | 0.8885, p = 0.001429 | MATCH |
| (b) Cb / Hip GC ISS-T − LAR rPRG | 0.704 (p = 0.005) / 1.118 (p = 0.010) | 0.7044 / 1.1180, same p | MATCH |
| (b) FLT ISS-T − LAR rPRG pooled / Cb / Hip | 0.656 / 0.332 (p = 0.21) / 1.004 (p = 0.0005) | 0.6558 / 0.3322 (0.2095) / 1.0045 (0.00051) | MATCH |
| (b) FLT − GC within ISS-T, Cb / Hip | −0.206 (0.36) / 0.164 (0.49) | −0.2060 (0.3556) / 0.1645 (0.4934) | MATCH |
| (b) FLT − GC within LAR, Cb / Hip | 0.188 (0.38) / 0.281 (0.06) | 0.1885 (0.380) / 0.2814 (0.060) | MATCH |
| (b) pooled GC rPRG bootstrap interval | 0.595–1.183 | 0.571–1.198 (own seed, 1,000 draws) | DIRECTION (Monte Carlo) |
| (b) NEW pooled FLT − GC within LAR (not in package) | – | rPRG +0.235 (exact p = 0.115, 400 assignments; bootstrap +0.023 to +0.435); HSP +0.065; SRG +0.087 | NEW |
| (b) NEW pooled FLT − GC within ISS-T | – | rPRG −0.027 (p = 0.886, 7,056 assignments; −0.341 to +0.285); HSP +0.123; SRG +0.338 (p = 0.027) | NEW |
| (c) ‖γ‖ for all 15 anchor–target geometries | 8.129–15.309 | identical to 3 dp for all 15; ‖I − γwᵀ‖₂ computed directly equals ‖γ‖; w·γ = 1 | MATCH |
| (k) NEW random 30-gene-sensor null for ‖γ‖ | – | OSD-561: Sensor31 15.31 vs random median 23.06 (2.5–97.5%: 15.5–38.6); OSD-562: 11.30 vs 22.45; OSD-525: 8.13 vs 20.83. Isotropic fixed-w expectation (Eq. 11) 47.7 / 47.2 / 37.1 | NEW |
| (d) OSD-561 / OSD-562 sensor AUC (81 splits each) | 1.000 / 1.000 | 1.000 / 1.000 (minimum 1.000) | MATCH |
| (d) train–test whole-genome cosine | 0.228 / 0.521 | 0.2277 / 0.5214 (per split \|Δ\| ≤ 1.7e-16) | MATCH |
| (d) original-score MSE ratio; outside energy removal | 1.091, −23.26% / 0.962, 16.85% | 1.0909, −23.26 / 0.9619, 16.85 | MATCH |
| (d) whole-genome AUC; sensor cosine (Supplement) | 0.796, 0.683 / 1.000, 0.843 | same | MATCH |
| (e) intersection–union screen: tests; BH/BY survivors | 54,989; 0/0 | 54,989; 0/0 (minimum BH q = 0.652) | MATCH |
| (e) concordant / nominal by region | Cb 5,047/3; Hip 5,001/0; LH 9,884/554 | identical | MATCH |
| (e) NEW: same screen with own log-expression OLS, left hemisphere | – | 9,865 concordant, 399 nominal, 0 BH | NEW |
| (f) Hbb-bt / Hbb-bs cerebellar interaction (NB, log2) | −3.495 / −3.009 | −3.495 / −3.009 (own NB2 MLE, own Cox–Reid dispersion; my gene-wise dispersion gives larger SE, p 1.7e-7 / 3.6e-6 vs package 1.0e-8 / 4.2e-7); OLS −3.654 / −3.112 | MATCH (coefficient); q not recomputed |
| (g) hippocampal Circadian interaction | −0.405 | −0.4053 (12 genes, 27 animals) | MATCH |
| (g) bootstrap interval (4,000) | −0.591 to −0.215 | −0.595 to −0.216 (own seed) | DIRECTION |
| (g) HC3 restricted wild p (499,999 draws) | 0.000838 | 0.000762 (package MC 95%: 0.000758–0.000920) | DIRECTION |
| (g) exclusion; one-animal and one-gene omission ranges | −0.460; −0.459 to −0.365; −0.452 to −0.360 | identical | MATCH |
| (h) external count columns; primary units | 66; 64 = 46 + 18 | 66; 64 = 46 + 18 (OSD-611: 12 columns, 10 in GLDS-592 SampleTable) | MATCH |
| (i) constrained mean MSE ratio / pooled SSE / improved, 7 settings | Table 3 | identical, recomputed from per-holdout SSE rows; max gain 2.000000 | MATCH |
| (i) error reduction ranges | 1.2–5.2%; 1.4–7.2% | 1.22–5.16; 1.43–7.16 | MATCH |
| (j) NEW own nested constrained regression, raw STAR counts | OSD-202 0.9284, 352 0.9857, 374 0.9501, 611 0.9719 | 0.9284 (37/40), 0.9857 (5/6), 0.9501 (8/8), 0.9719 (9/10); max per-unit \|Δ\| ≤ 4.8e-5 | MATCH |
| (m) NEW own nested constrained regression, core 81 + 81 holdouts | 0.9878 (78/81), 0.9484 (81/81) | 0.9878 (78/81), 0.9484 (81/81); max \|Δ\| ≤ 3.0e-5 | MATCH |
| (l) NEW margin over matched clipped ridge (per holdout) | "lower in every setting" | Δ mean MSE ratio −0.0018, −0.0145, −0.0046, −0.0226, −0.0079, −0.0105, −0.0020; constrained better in 56/81, 81/81, 185/200, 34/40, 6/6, 7/8, 3/10 | NEW (supports "on average", not "uniformly") |

Not independently rerun (reported as is, with reasons):
- The PyDESeq2 fits (pydeseq2 is not in the environment). Component Wald p values of the screen and the 72,023-test family q values were not refitted.
- The 1,400 count simulations and the 2,600 Gaussian datasets. Only their row counts and summary values were checked.
- The 400-assignment exact label test (p = 0.005).
- The 200 RR-23 holdouts.
- The 1,488-row biological family and the 104 paired tests. Only their counts, effects and q values were checked against files.

## 3. Design audit (`A_design_matrix.csv`, `A_contrast_aliasing.csv`, `A_design_library_level.csv`, `A_location_conflicts.csv`)

`A_design_matrix.csv` cross-tabulates all 7 core and 4 external deposits from the deposited fields, one row per combination of these factors with library and unit counts: exposure, collection cohort, euthanasia method, deposited and design euthanasia location, preservation (raw string and bulk-RNA-seq class), carcass preservation, dissection, euthanasia and dissection date, recovery, habitat, sex, strain and age. `A_contrast_aliasing.csv` scores every contrast the paper uses, factor by factor, as matched / ALIASED / partially overlapping / not deposited.

| Contrast (paper's label) | Matched by deposited fields | Aliased with the contrast |
|---|---|---|
| GC ISS-T − GC LAR (ground "collection") | sex, strain, age (balanced), kill location (Earth), bulk preservation (LN2), habitat | **euthanasia agent** (ketamine/xylazine overdose vs gas anaesthesia); recovery field (LAR 24 d); dates partly offset (18–20 vs 18 Sep 2019) |
| FLT ISS-T − FLT LAR | sex, strain, age | **agent; kill location (ISS vs Earth); freezing (Cryochiller vs LN2)**; time on orbit (55–58 vs 32 d) + 24 d Earth recovery |
| FLT − GC within ISS-T (Cb, Hip; also OSD-613) | agent (ketamine/xylazine), sex, strain, age | **kill location (ISS vs Earth), freezing (Cryochiller vs LN2)**; dates offset (16–19 vs 18–20 Sep); time of day at euthanasia not deposited |
| FLT − GC within LAR | agent (gas), location (Earth), LN2, 24 d recovery, sex, strain, age | flight history; euthanasia date (19 vs 18 Sep 2019) |
| RR-10 FLT − GC (563, 564, 612) | agent (ketamine/xylazine), carcass Cryochiller, LN2 after BSP dissection, sex, strain | **kill location (ISS vs Earth)**; euthanasia date (2 vs 5 Jan 2021); BSP dissection date (12 vs 14 May 2021). Times of day are clock-matched (OSD-612: 09:34–13:26 GMT vs 04:34–09:10 EST) |
| RR-23 HGC − VGC (ground anchor) | agent (isoflurane), RNAlater, euthanasia and dissection date (17 Jan) | habitat (flight hardware vs vivarium cage) and transport history |
| RR-23 FLT − HGC / FLT − VGC | agent, RNAlater | **euthanasia and dissection date (14 vs 17 Jan 2021**; FLT 1 day after splashdown); habitat for FLT − VGC |
| RR-3 FLT − GC (OSD-352) | ketamine/xylazine/acepromazine, LN2, frozen-carcass BSP dissection | sacrifice dates (17/19 vs 22/23 May 2016); RIN fully separated (6.3–6.8 vs 7.5–8.1); location not deposited |
| OSD-202 / 374 / 611 (Earth radiation) | within-study protocol fields uniform | euthanasia not deposited for OSD-202; dates not deposited |

In every RRRM-2 contrast the paper interprets, "collection cohort" bundles concrete protocol factors that the text never names. No RRRM-2 contrast varies the euthanasia agent with flight held fixed: within ISS-T both groups got ketamine/xylazine, within LAR both got gas. So **agent and freezing effects are separable from flight only across cohorts, never within one.** The cleanest flight contrast is FLT − GC within LAR, where everything deposited is matched except a 1-day date offset.

Older audit issues (`Paper_1_Complete_Submission_Package/REVIEW_AND_CORRECTIONS.md`), Paper 1 items:
1. **Author/affiliation:** the old names are gone (0 hits for Zainulabedin, Khan, Mumbai). The current front matter lists 5 authors, "Convolity AI Lab, Lucknow" and an email. The earlier instruction was to leave affiliation and email blank, so Isarar should confirm this is now intended.
2. **Publication identifiers:** resolved in the .md (no s41526-026 DOI, no received/accepted dates). The LaTeX template must not reintroduce them.
3. **Vivarium date:** partially resolved. "20 January" is gone, but a new incorrect claim of different HGC/VGC harvest dates has appeared (§1), and the flight/date alias (14 vs 17 Jan) is still not stated.
4. **Mixed regression summaries (22 contrasts, R² 0.6528, F test):** resolved; absent from the text.
5. **Kinetic/perimortem window:** resolved. The Introduction now disclaims "time since euthanasia", and no 15–60 min claim remains.
6. **221-check coverage:** superseded. The text now cites 1,533, 96 and 1,337 verifier checks, which are still the package's own self-checks. §2 is the independent check of the headline numbers.

Items 2, 5, 6, 8, 9 and 10 of that file concern Paper 2.

## 4. Overclaim scan (wording that exceeds the evidence → accurate rewording)

1. **Title:** "Gain constrained regression improves transcriptomic prediction across spaceflight and Earth brain datasets" → "A gain-capped sensor regression gives small held-out prediction gains in seven mouse-brain transcriptome settings". The gains are 1.2–7.2% versus a nuisance-only baseline. Unbounded ridge is better by 5–25 points in 5/7 settings, and the margin over clipped ridge is 0.2–2.3 points. A reframed title around the protocol finding is preferable.
2. **Abstract:** "a convex problem with verifiable optimality and stability bounds" → "a convex problem whose solutions carry numerical optimality certificates and a hard cap of 2 on worst-direction amplification". The "bound" is an imposed constraint, not a derived stability property.
3. **Introduction:** "The contribution is a validated stability-constrained analysis for these data" → "an audited, gain-capped analysis of these data; no step validates biological correction".
4. **Keywords:** "independent datasets; nested validation" → "external datasets (refitted within cohort); nested cross-validation".
5. **Results heading:** "Ground collection cohorts share response programs with flight contrasts", and "collection-cohort-associated expression" throughout → name the factors: "In unflown animals, the ISS-T protocol (ketamine/xylazine overdose) versus the LAR protocol (gas anaesthesia) was associated with +0.889 log2 rapid primary-response expression; the deposited fields cannot separate the agent from other cohort procedures."
6. **Results heading:** "A direct flight by collection interaction identifies a qualified circadian candidate" → "A circadian-panel difference between ISS-killed and Earth-killed animals that is specific to the ISS-T cohort". Within ISS-T, FLT − GC is aliased with kill location and Cryochiller versus LN2, and time of day at euthanasia is not deposited for RRRM-2. A phase difference at euthanasia would produce exactly this panel signal, so "flight" is not identified.
7. **Results heading:** "Constrained multi-sensor regression improves held-out prediction" → "…reduces held-out error by 1.2–5.2% while capping amplification at 2 (ordinary ridge: 6.0–13.4% at mean gains 12–29)".
8. **Results heading:** "The constrained algorithm generalizes to four external prediction settings" → "Refitted within four further cohorts (6–40 units each), the constrained estimator reduced mean held-out error by 1.4–7.2%". "Generalizes" implies transport of a fitted model; each cohort was refitted.
9. "Constrained regression had lower equal-split/unit mean MSE than matched clipped ridge in every setting" → add the per-unit picture: "by 0.2–2.3 points on average; it was lower in 56/81, 81/81, 185/200, 34/40, 6/6, 7/8 and 3/10 holdouts".
10. **Discussion:** "an independent full-coefficient solver verified the optimization and operator identities" → "an independent conic solver reproduced 12 small instances (maximum objective difference 2.5×10⁻⁹)".
11. **Discussion:** "stable prediction can be made constructive and auditable even when causal nuisance separation remains unidentified" → "a hard gain cap can be imposed at a modest cost in prediction; it neither identifies nor removes nuisance". The cap bounds the operator, not the stability of predictions.
12. "Direction concordance alone was frequent: 5,047 cerebellar, 5,001 hippocampal and 9,884 left-hemisphere genes" → "Direction concordance was at chance level: 28.0% and 28.2% of genes for three contrasts (25% expected) and 51.5% for two (50% expected)" (computed from the screen: 5,047/18,048, 5,001/17,759, 9,884/19,182).
13. "The learned spectral gains ranged from 8.129 to 15.309…" → add: "random 30-gene sensors on the same anchors give larger gains (medians 21–23), so the amplification reflects dimension and residual rank, not a property of the response signature" (§2 k).
14. **Count benchmark:** "With 24 training samples and changed loading, outside-sensor RMSE was 0.426 for constrained regression, 0.473 … 0.513 … RUVg 0.403" → state that unadjusted analysis had 0.4255, the same as constrained regression (`benchmark_summary.csv`, correlated_changed_loading, n_train 24). The constraint matched no correction; it did not improve on it.
15. "HGC-versus-VGC contrast includes different handling and harvest dates" → false; see §1.
16. "Seven location entries … disagreed" → state that 5/6 animals contradict themselves across deposits and that the OSDR description says ISS-T mice were sacrificed on orbit.
17. **Table 2 / external text:** RR-3 is described with RIN separation only → add that flight and ground sacrifice dates also differ (17/19 vs 22/23 May 2016).
18. "Sensor31 retained AUC 1 in both regions (exact p = 0.005…)" → state that each split's AUC rests on two within-age pairs (four test animals), and that p = 0.005 is 2/400: the observed assignment plus one other reached AUC 1 (`n_null_at_least_observed = 2`).
19. "All external datasets were chosen … before their new analysis outcomes were examined" and "Hyperparameter selection … fixed before these outer results were inspected" → NOT VERIFIED: no time-stamped protocol predating the outputs was found. Keep only as a statement of intent, or cite the protocol file hash and time.

## 5. Strength ranking (evidence × novelty × reviewer interest)

| # | Finding | Exact numbers and test | Weakest point | Shown in |
|---|---|---|---|---|
| 1 | Ground-only protocol (euthanasia-agent-aliased) IEG effect | GC ISS-T − GC LAR rPRG +0.889 log2, bootstrap 0.576–1.214; exact whole-animal permutation within age blocks p = 0.00143 (all 700 assignments); BH q = 0.019 (232 tests). Cb +0.704 (p = 0.005), Hip +1.118 (p = 0.010); HSP +0.332 (p = 0.011); SRG +0.032 (p = 0.84) | 13 animals, one mission; the agent is aliased with any unrecorded ISS-T ground-schedule procedure | Fig. 2 / Fig. 1b; `results/panel_results.csv`; reproduced exactly |
| 2 | Matched-protocol flight effects are small and non-significant against a large protocol effect | pooled FLT − GC within LAR +0.235 (p = 0.115; bootstrap +0.023 to +0.435); within ISS-T −0.027 (p = 0.886); OSD-613 +0.019 (p = 0.93); RR-10 Cb +0.204, Hip −0.047, LH +0.046 (all p ≥ 0.65). The protocol-effect lower bound (0.576) exceeds the LAR flight upper bound (0.435) | Absence of evidence; small cells (3 per age × group in LAR); percentile bootstrap narrower than the exact test; the two intervals share GC LAR animals | New (pooled rows not in package); `A_recompute_panel_effects.csv` |
| 3 | Housing/return handling moves the IEG/HSP panel in RR-23 | HGC − VGC rPRG +0.641 (p = 3.3×10⁻⁴, q = 0.012), HSP +0.726; FLT − HGC rPRG −0.399 (q = 0.017), HSP −0.993 (q = 0.011), SRG −0.369 (q = 0.020) | FLT − HGC is aliased with collection date (1 vs 3 days after return); one cohort; male C57BL/6J | Fig. 2; `results/panel_results.csv` |
| 4 | Perfect within-resource discrimination with poor transfer of genome-wide loadings | sensor AUC = 1.000 in 81/81 and 81/81 splits (label-enumeration p = 0.005, BH q = 0.011 of 18); whole-genome train–test cosine 0.228 / 0.521 / 0.165; original-score prediction worse than baseline in Cb (MSE ratio 1.091) and RR-23 (1.035); outside energy −23.3% in Cb | Test contrasts from 4 animals (df = 1), so the cosine is attenuated by noise; no noise ceiling was computed | Fig. 6; `extensions/ground_holdout/holdout_summary.csv`; reproduced per split |
| 5 | Score removal is an oblique projector with gain 8–15, generic to p ≫ n | ‖γ‖ = 8.129–15.309 for all 15 geometries (exact identity ‖A‖₂ = ‖γ‖); random-sensor medians 21–23; isotropic expectation 37–48 | The norm is a worst-direction bound; random-sensor null is new and not in the paper | Fig. 5a; `results/operator_geometry.csv`; `A_recompute_operator_gains.csv` |
| 6 | Gain-capped regression: small, consistent, fully reproducible gains at a bounded gain | mean MSE ratio 0.9878, 0.9484, 0.9628, 0.9284, 0.9857, 0.9501, 0.9719; improved 78/81, 81/81, 200/200, 37/40, 5/6, 8/8, 9/10; gain ≤ 2; reproduced from raw counts with independent code. Margin over clipped ridge 0.18–2.26 points | Gains are small; unbounded ridge predicts better in 5/7 settings (e.g. OSD-374 0.704 vs 0.950); the target is outside-sensor residual expression, not a biological endpoint; core splits overlap (30 animals) | Fig. 7–8, Table 3 |
| 7 | No reproducible flight-response gene across cohorts | 0/54,989 BH and BY (minimum q = 0.652); concordance at chance (28% vs 25%; 51.5% vs 50%); reproduced from component p; holds with an OLS estimator (0 BH, left hemisphere) | Component contrasts are themselves protocol-aliased (both left-hemisphere contrasts are ISS-killed vs Earth-killed) | Suppl. Fig. 1; `results/regional_gene_replication_summary.csv` |
| 8 | Beta-globin interaction tracks kill site | Cerebellar Hbb-bt −3.495, Hbb-bs −3.009 log2 (global q 0.000197 / 0.0014); FLT − GC in ISS-T −1.73 / −1.39 vs in LAR +1.76 / +1.62. ISS-killed animals carry less globin mRNA than Earth-killed controls, which fits exsanguination/dissection differences | Blood content not measured; hippocampus not significant (q 0.32 / 0.54) | Suppl. Fig. 4; `extensions/flight_collection_interaction/` |
| 9 | Hippocampal circadian-panel interaction | −0.405 log2 (bootstrap −0.591 to −0.215; wild p = 0.000838; q = 0.027 of 32; new family q = 0.036, no BY survivor); all one-animal/one-gene omissions negative; count model −0.367 | Time of day at euthanasia not deposited for RRRM-2; aliased with kill site and freezing; exclusion q = 0.786; RR-3 opposite sign (+0.023) | Suppl. Fig. 3, 6 |
| 10 | Exact identity and accounting audit | 129 → 126 → 83 animals; 24 + 10 pairs; STAR counts match upstream; RFIDs agree across regions for all 24 pairs; 5/6 conflict animals contradict themselves across deposits | Mostly a methods asset, not a finding | Table 1, Table 6, `A_location_conflicts.csv` |

How to use this ranking: findings 1–4 plus 7 support the "protocol, not flight" paper. Findings 5–6 are its constructive half: gain-capped prediction, with the honest caveat that the cap costs accuracy and does not identify nuisance. Findings 8–9 are protocol-aliased and should be presented as such, not as flight biology.

## 6. Remaining work and limits

- **PyDESeq2 refits not performed** (the package is not installed; install with a pinned version and refit OSD-561/562 to check the gene-level q families).
- **Simulations not rerun:** count benchmark (1,400 datasets) and Gaussian benchmarks (2,600 datasets). They are only cross-checked at summary level.
- **Exact AUC label test not rerun:** the 400-assignment test and the RR-23 holdouts were not rerun.
- **No noise ceiling for the holdout cosine:** for example, a split-half cosine of ground contrasts.
- **Pre-specification claims not verified:** statements that choices were fixed before outcomes were seen have no verifiable timestamp in the package.
- **External raw-file timestamps:** the files in `extensions/strengthening_20261015/external_data/*/raw_sources/` have filesystem dates of 11 Sep, while the references say "accessed 15 October 2026". This is probably preserved server timestamps; it is not verified.
- **No literature claims:** none are made here. Novelty and literature support (for example, prior work on anaesthetic effects on immediate-early genes) belong to Worker B.

## 7. Files

| File | Content |
|---|---|
| `A_REPORT.md` | this report |
| `A_claims_check.csv` | 591 rows: id, section, sentence, number_as_written, source_file, source_value, status, note, auto_status |
| `A_claims_check.py` / `A_claims_hand.py` | extraction and automatic matching / hand-check rules (each rule reads the named file) |
| `A_recompute.py` / `A_recompute_results.json` | independent recomputation (sections a–m) / 143 rows: id, claim, expected, recomputed, status, note |
| `A_recompute_panel_effects.csv` | all recomputed panel effects with exact p and package values side by side |
| `A_recompute_operator_gains.csv`, `A_recompute_holdout_OSD-561.csv`, `A_recompute_holdout_OSD-562.csv` | per-geometry gains; per-split holdout metrics |
| `A_design_matrix.py` / `A_design_matrix.csv` / `A_design_library_level.csv` | protocol-factor cross-tabulation (cell-level and library-level) |
| `A_contrast_aliasing.csv` | per-contrast, per-factor matched/aliased status |
| `A_location_conflicts.csv` | the 6 location-conflict animals with both regional records, RFID and date agreement |

Python: `/Users/zain/.kiro/crew/scratch/runtime-735c0e1b/p1env/bin/python`, BLAS single-threaded; no new packages installed. `A_claims_check.py` caches its 6.1 M-value index at `$KIROCREW_SCRATCH/A_claims_value_index.pkl` on first run (about 2 min to build; the cache was deleted after this run).
