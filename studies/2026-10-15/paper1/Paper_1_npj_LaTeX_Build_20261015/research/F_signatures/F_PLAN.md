# F — Published-signature test: analysis plan (written before any signature result was computed)

Status when written (15 Oct 2026): the gene lists below were built from the source files and their
overlap with the package panels was counted. No expression value has been read for any of these lists.
The package rPRG/HSP results (research/D_protocol) were already known. The SHA-256 of this file and of
`signatures.json` is recorded in `F_PLAN.sha256` before `F_analysis.py` is run. The analysis script refuses
to run if either hash differs. The timestamp is local and not externally certified.

## Question

The paper's lead finding uses the package's own rapid primary-response panel (rPRG). Do gene sets defined
by other laboratories, for other purposes, reproduce the same pattern? The pattern is a collection-protocol
effect in never-flown ground mice that exceeds the protocol-matched flight effect.

## Gene sets (signatures.json; source files and SHA-256 in sources/)

| ID | Source | Genes | Role |
|---|---|---|---|
| S1 | Hrvatin et al. 2018, Suppl. Table 3: early-response genes, increase (label a) in at least 1 of 30 cell types | 269 | confirmatory |
| S2 | Hrvatin et al. 2018, Suppl. Table 3: late-response genes, increase (label c) in at least 1 cell type and never early | 152 | time-scale control |
| S3 | Marsh et al. 2022, Suppl. Table 4: All CNS Cells Act. Score (ex vivo activation) | 19 | confirmatory |
| S4 | Marsh et al. 2022, Suppl. Table 4: Micro/Myeloid Shared Act. Score | 25 | confirmatory |
| S5 | van den Brink et al. 2017, SI R vector of dissociation-affected genes (identical to the copy in Marsh 2022 Suppl. Table 21) | 140 | confirmatory |
| S6 | Wu et al. 2017, Table S3: List of IEGs | 139 | confirmatory |

Two versions of each set:
- **full**: all mapped genes.
- **novel**: excluding every gene in the package's seven base panels (rPRG, dPRG, SRG, HSP,
  Glia_Complement, BBB_Endothelium, Circadian; 107 symbols). These versions share no gene with the package
  panels.

A version is tested only if at least 5 genes map in the analysed deposit(s).

## Mapping

- S1 and S2 are mapped by the source's Ensembl ID to count-matrix rows.
- S3 to S6 are mapped by exact symbol match through the deposit annotation, which is the package rule;
  there is no alias substitution.
- A gene must pass the package count filter (mean raw count of at least 10 in the deposit). Pooled analyses
  use genes present in both OSD-561 and OSD-562.
- Unmapped symbols are listed in the report.

## Estimator and inference (identical to research/D_protocol)

- **Estimator:** the signature effect is the mean, over mapped genes, of the OLS group coefficient on
  log2(count/size factor + 1) after adjusting for age and, when pooled, tissue (package equation 1).
- **Inference:** two-sided exact p from complete enumeration of whole-animal labels within age blocks, with
  all regions of an animal kept together. This uses package `animal_permutations`; a Monte Carlo plus-one p
  is used only if there are more than 100,000 assignments.

RRRM-2 contrasts:

| Contrast | Definition | Meaning |
|---|---|---|
| P_G | GC ISS-T − GC LAR | protocol effect, no flight |
| F_L | FLT LAR − GC LAR | protocol-matched flight effect |
| D | GC ISS-T − FLT LAR | dominance: positive means P_G > F_L relative to the shared GC LAR reference |
| P_F, F_I | | descriptive |

- **Analyses:** cerebellum (OSD-561), hippocampus (OSD-562), and pooled for P_G, P_F and D.
- **Cross-mission, descriptive only:** RR-10 FLT − GC (OSD-563, -564, -612); RR-23 FLT − HGC, FLT − VGC
  and HGC − VGC (OSD-525).

## Reproduction gate (must pass before any signature number is computed)

- The rPRG P_G, P_F and D effects and exact p, recomputed by `F_analysis.py`, must equal
  `research/D_protocol/D_tables` (|effect difference| < 1e-12, |p difference| < 1e-12) in all three
  analyses.
- The RR-10 and RR-23 rPRG effects must equal package `results/panel_results.csv` (design_cohort rows).

## Families and decision rules (fixed here)

- **Primary family:** pooled D for S1, S3, S4, S5 and S6, full and novel versions (up to 10 tests).
  BH-adjusted. Hypothesis: D > 0.
- **Secondary family:** pooled P_G for the same versions, BH-adjusted.
- **Per-region results** are descriptive; they are reported in full but enter no decision.

| Rule | Criterion |
|---|---|
| R1 | Supported if at least 3 of the 5 full confirmatory versions have pooled D > 0 with BH q ≤ 0.05 in the primary family. |
| R2 | The effect extends beyond the package panel genes if at least 3 of the testable novel versions have pooled D > 0 with q ≤ 0.05. Otherwise the paper must say that independent signatures confirm the pattern only through the canonical IEG/HSP genes they share with the package panels. |
| R3 | Time-scale specificity holds if S2 pooled P_G has exact p > 0.05 and S2 P_G < S1 P_G (full versions). |

If R1 fails, the paper states that independently defined signatures do not reproduce the dominance pattern,
and the lead claim is restricted to the package rPRG panel.
