# Worker B — novelty and prior-art verdicts for Paper 1

15 October 2026. Searches run 02:07–02:30 UTC; 198 logged queries/fetches (`B_search_log.json`). 32 prior works recorded, 29 read at source with verbatim quotes machine-checked against the cached text (`B_closest_prior_art.json`; mutation test: altered quotes are rejected).

## Headline

1. E1 (empirical) survives, **narrowed**. Collection confounding in rodent spaceflight RNA-seq is not new (Lai Polo 2020), and neither is anaesthesia or handling changing brain early-response genes (Staib-Lasarzik 2014; Marsh 2022). What is new, as far as these searches show, is the size of the effect measured in brain in **flight-naive** animals under the two RRRM-2 protocols, compared with flight effects matched for protocol.
2. T1 (theoretical) is **killed as stated**. The oblique-projector gain is a classical identity (Szyld 2006). The cosine evidence for "non-transferable loadings" does not hold: an independent re-run shows that the observed train–test cosines (0.22 / 0.52 / 0.15) are what an *identical* loading gives at these sample sizes (expected 0.21 / 0.50 / 0.15). What survives is a narrower claim about precision and gain (T1′ below). **This changes a conclusion in the current manuscript.**
3. M1 and C1 are **narrowed** to claims about application and reporting. Every component is already established. The one item with no prior hit in our searches is reporting a correction operator's spectral gain as a standard diagnostic.
4. P1 is **supported**, with a caveat about search sensitivity. OSDR links no publication to OSD-561/562/563/564/612/613/525. Europe PMC and PubMed return no paper that cites these deposit IDs.

## (1) Verdict table

| Claim | Verdict | Exact sentence for the paper (≤45 words) | Bounding references | Prior work that narrowed or overturned it |
|---|---|---|---|---|
| E1 | **NARROWED** (supported as worded) | "In flight-naive RRRM-2 ground controls (12 mice per region), the ISS-terminal cohort exceeded the live-animal-return cohort on a rapid primary-response gene panel by 0.70 (cerebellum) and 1.12 (hippocampus) log2 units, more than protocol-matched flight effects on that panel (−0.21 to +0.28)." | Lai Polo 2020; Staib-Lasarzik 2014; Tyssowski 2018; Grandke 2026 | Lai Polo 2020 already showed that collection technique eclipses flight effects (liver, quadriceps), so "first collection confound" is gone. It also found that ketamine/xylazine vs CO2 did not alter liver or muscle, and Chakkingal 2023 found no euthanasia-agent DEGs in pig hypothalamus/pituitary, so the agent cannot be named as the cause. Staib-Lasarzik 2014 and Marsh 2022 show that the class of artefact (early-response genes moving with anaesthesia or handling) is known. Grandke 2026 noted a TERM/LAR extraction effect in RRRM-1 miRNA. Dropped from the draft: "shifts" (causal), "any protocol-matched flight effect" (only this panel was compared), and folding "no gene replicates" into the novelty claim (see the note below the table). |
| T1 | **KILLED** as a theory claim and as a non-transfer claim. Replacement T1′ is **NARROWED/SUPPORTED** | T1′: "Near-perfect discrimination by a 30-gene nuisance sensor (mean AUC 0.96–1.00) did not make its genome-wide loading usable: with 8–11 training animals, train–test loading cosines (0.15–0.52) matched those expected for an identical loading, and score-removal operators had spectral gain 8–15." | Szyld 2006; Miao 2018; Wang 2017; Gagnon-Bartsch & Speed 2012 | Szyld 2006: ‖P‖ = ‖I−P‖ for any oblique projector, "proven repeatedly". The score-removal gain ‖g‖·‖w‖ follows directly from it, so there is no new theorem. Miao 2018: a single proxy does not identify the confounding adjustment. Wang 2017: genome-wide adjustment needs negative controls or sparsity. So "a reliable proxy does not license correction" is already established theory. The noise-ceiling check (below) removes the empirical "non-transfer" reading. |
| M1 | **NARROWED** | "We use a mission's own flight-naive animals, collected under both protocols, as a ground-only holdout with training-only preprocessing, nested tuning and a perfect-transfer noise ceiling, and fit sensor regressions under an explicit spectral-gain cap solved to a certified duality gap." | Varma & Simon 2006; Kapoor & Narayanan 2023; Gouk 2021; Ye 2023 | Nested CV (Varma 2006), leakage control (Kapoor 2023), projection under a Lipschitz/gain bound (Gouk 2021) and design-based batch remedies (Ye 2023) are all established. Lai Polo 2020 already used ground experiments to estimate collection effects. The "noise ceiling" clause is valid only if the parent adds that analysis (script provided). Without it, the holdout cannot tell non-transfer apart from imprecision. |
| C1 | **NARROWED** | "We propose reporting the spectral gain of every post hoc correction operator, the worst-case factor by which it can amplify residual noise, as a routine diagnostic; in a targeted search, we found omics correction studies reporting removal and calibration but not this quantity." | Gouk 2021; Szyld 2006; Gerard & Stephens 2021; Nygaard 2016 | A gene panel used as a proxy for unwanted variation is the RUV control-gene idea (Gagnon-Bartsch & Speed 2012), here with positive controls, so the "sensor" framing is not a novel concept. Lipschitz-bounded fitting is standard in ML (Gouk 2021). The "small-gain theorem" wording should be **dropped**: the correction is a feed-forward linear map with no feedback loop, so operator norm or Lipschitz constant is the correct term. Correction-induced false positives and overconfidence are established (Nygaard 2016; Zindler 2020). |
| P1 | **SUPPORTED** (absence claim, search-limited) | "To our knowledge, based on OSDR publication links and targeted Europe PMC, PubMed and Crossref searches on 15 October 2026, the RRRM-2, RR-10 and RR-23 brain bulk RNA-seq deposits analysed here (OSD-561, -562, -563, -564, -612, -613, -525) have no prior peer-reviewed analysis." | Overbey 2022; Grandke 2026; Masarapu 2024; Choi 2020 | Nothing overturned it. Closest: Overbey 2022 (RRRM-2 design; scRNA-seq of LAR immune tissues), Sakilam 2026 preprint (RRRM-2 *kidney*), Grandke 2026 (RRRM-1 small RNA incl. brain), Masarapu 2024 (RR-3, OSD-352), Rodriguez 2026 (RR-1 brain, protein). P1 is a resource fact, not a scientific advance, and should not lead. |

Note on "no gene replicates across cohorts after correction" (package `results/regional_gene_replication_summary.csv`, read, not re-run): 0 genes at BH or BY FDR < 0.05 in cerebellum (18,048 genes, 3 contrasts), hippocampus (17,759; 3) and left hemisphere (19,182; 2). With 3–5 animals per group this is a null with low power. Report it as "no gene reached FDR < 0.05 in every contrast", never as "no reproducible flight effect exists".

## New finding that changes the manuscript: the cosine "transfer failure" is explained by sampling noise

`analysis/noise_ceiling.py` and `analysis/noise_ceiling_OSD525.py` use ground animals only, with the package normalisation. For each deposit they compare the observed train–test cosine of the whole-genome nuisance contrast with the cosine expected if training and test animals shared **exactly the same** loading, given the per-gene residual variance: E[cos] ≈ S/√((S+N_tr)(S+N_te)). A Gaussian simulation confirms the formula.

| Deposit (design) | Observed cosine, this re-run (package value) | Expected under identical loading | Genome-wide SNR, test split / train split | Package outside-sensor MSE ratio (score) |
|---|---|---|---|---|
| OSD-561 cerebellum (train 8 / test 4) | 0.224 (0.228) | 0.211 | 0.19 / 0.38 | 1.09 (worse) |
| OSD-562 hippocampus (8 / 4) | 0.519 (0.521) | 0.502 | 0.73 / 1.47 | 0.96 (better) |
| OSD-525 RR-23 HGC vs VGC (11 / 6) | 0.152 (0.165) | 0.148 | 0.13 / 0.24 | 1.03 (worse) |

Sensor31 genes on their own: observed 0.69 / 0.84 vs expected 0.66 / 0.82 (consistent with identical loadings, SNR 1.4–3.3).

- **What this means:** the data cannot tell a non-transferable loading apart from a transferable but imprecise one. Perfect AUC reflects the 30-gene sensor's high SNR. The genome-wide loading has SNR < 1 in two of three deposits. Score removal applies that noisy loading with gain 8–15 (full-data, `results/operator_geometry.csv`; holdout means 9.3–18.2).
- **Constructive half (hypothesis from 3 points; needs formal check by Worker D/A):** in all three deposits, the score correction worsened held-out prediction where training SNR was < 1 (OSD-561, OSD-525) and improved it where it was > 1 (OSD-562). That suggests a pre-correction rule: estimate the loading SNR from ground controls, and correct only when it exceeds about 1. Otherwise, match collection by design.
- **Reproduced in passing:** rPRG ground protocol effect +0.696 (cerebellum; package +0.704) and +1.127 (hippocampus; package +1.118). Close in magnitude, same direction, small difference from the gene-filter scope.
- **Limits:** the noise-only Gaussian model ignores gene–gene correlation; S is estimated from the same 12 animals; the RR-23 splits are this re-run's own 200, not the package's.

## (2) Lead claims and novelty text

Lead with **E1** (the inversion: in brain, the IEG signal tracks collection protocol, not flight), then **T1′ + C1** (the constructive half: why post hoc removal fails at these sample sizes, and the gain and SNR to report). Use **P1** in Methods and the cover letter only.

**Introduction (117 words):**
Sample collection is known to confound rodent spaceflight transcriptomes, and anaesthesia or tissue handling can alter early-response genes in brain. What has not been measured is how large a collection-protocol effect is in brain relative to spaceflight itself, using animals that never flew. The RRRM-2 mission made this measurable: its flight-naive ground controls were collected under both the ISS-terminal and the live-animal-return protocols. Here we show that, in these ground controls, the protocol difference on a rapid primary-response gene panel exceeds the four protocol-matched flight effects on the same panel. We then ask whether such a signature can be removed post hoc, and find that the binding limit is estimation precision, quantified by operator gain, not detectability.

**Cover letter (75 words):**
We reanalyse brain RNA-seq deposits from RRRM-2, RR-10 and RR-23 that had no prior peer-reviewed analysis. In mice that never flew, the ISS-terminal versus live-animal-return collection protocol is associated with a 0.7–1.1 log2 difference in a rapid immediate-early gene programme, larger than protocol-matched flight effects. The signature is readily detected, but at available sample sizes it cannot be removed post hoc without amplifying noise eight- to fifteen-fold. Collection protocols must therefore be matched by design.

## (3) Positioning against the closest work

Lai Polo et al. (2020) established that collection technique can eclipse flight effects in rodent spaceflight RNA-seq. Their evidence came from RR-1 liver and dedicated ground experiments: carcass freezing and polyA selection drove the artefact, and ketamine/xylazine versus CO2 did not alter liver or quadriceps. They studied no brain and cautioned against extrapolating to other tissues. Paper 1 differs in four ways:

- tissue: cerebellum and hippocampus;
- gene class: a minutes-scale activity programme (Tyssowski 2018);
- variable: the ISS-terminal versus live-animal-return protocol, aliased with euthanasia agent; within the ground contrast both arms were LN2-frozen;
- comparator: protocol-matched flight effects in the same mission.

Holley et al. (2022) and Overbey et al. (2022) already flagged return and landing stress as a confound and proposed on-orbit collection. Grandke et al. (2026) saw a TERM/LAR extraction effect in RRRM-1 miRNA and read the LAR–TERM contrast as reentry stress. Kremsky et al. (2023) report Arc among the top flight DEGs in mice processed within 20 h of splashdown. Our ground-only measurement shows why such contrasts cannot be read as flight biology without protocol-matched controls.

Methodologically, the operator identity (Szyld 2006), proxy non-identification (Miao 2018), control-based adjustment theory (Wang 2017; Gagnon-Bartsch & Speed 2012) and correction-induced errors (Nygaard 2016; Zindler 2020) all precede this work. Paper 1's contribution there is to measure gain and loading precision in a real spaceflight resource.

What the paper must not claim:

- first to show collection or euthanasia effects;
- ketamine as the cause (the agent is aliased; Lai Polo and Chakkingal found null agent effects in other tissues);
- "non-transferable" loadings;
- the small-gain theorem.

## (4) Standing caveat

This was a targeted search, not a systematic review. Its sensitivity has three measured limits:

- Europe PMC full-text search sees only open-access full text it indexes. The positive control Masarapu 2024, which names OSD-352 in its Methods, is **not** retrievable that way.
- Rodriguez 2026 (published 4 Aug 2026) was not yet indexed in Europe PMC.
- Crossref bibliographic search is fuzzy.

Every novelty sentence above therefore carries "to our knowledge" or "in a targeted search". Before submission, each claim needs a full Europe PMC section search with positive controls, plus a Google Scholar or Dimensions pass, which these tools could not run.

What this pass overturned in the parent's own first drafts:

- "first collection confound" (by Lai Polo 2020);
- "new theory" (by Szyld 2006 / Miao 2018);
- "loadings do not transfer" (by this worker's own noise-ceiling computation);
- "small-gain" (by category error).

Other limits:

- Behrens & Scharf 1994 and Zames 1966 were verified as metadata only, so cite them for nothing more than their titles.
- Staib-Lasarzik 2014 was read as abstract only (the full-text API returned HTTP 500).
- Search logs and claim numbers taken from the package (permutation p, flight effects, gains, replication counts) were read at source. Only the rPRG ground effects and the cosines were re-run here.

## Files

- `/Users/zain/Downloads/OSD-561-562-brain-transcriptomics-main/DOWNLOAD_PACKAGES_2026/Paper_1_npj_LaTeX_Build_20261015/research/B_novelty/B_REPORT.md` (this report)
- `.../research/B_novelty/B_search_log.json`: 198 entries (endpoint, query, URL, UTC time, HTTP status, hit count, top hits, cache file). 11 early Europe PMC calls failed with HTTP 503, caused by the `sort=RELEVANCE` parameter. They were re-run without it; this probably also explains the 503s in the earlier `journal_audit_20261015` search.
- `.../research/B_novelty/B_closest_prior_art.json`: 32 entries (id, DOI, PMID, title, year, ≤60-word summary, verbatim quote, `quote_verified`, relation, claim ids).
- `.../research/B_novelty/analysis/noise_ceiling.py`, `noise_ceiling_results.json`, `noise_ceiling_OSD525.py`, `noise_ceiling_OSD525.json`
- `.../research/B_novelty/tools/bsearch.py` (logged query helper), `tools/build_prior_art.py` (builds and verifies the prior-art JSON)
- `.../research/B_novelty/cache/`: raw API responses, full texts and publisher pages (16 MB)
