# F — Published-signature test (15 Oct 2026)

Plan: `F_PLAN.md`, hashed at 2026-10-15T04:35:43Z before `F_analysis.py` was written (`F_PLAN.sha256`).
The script checks both hashes before running. Run log: `F_run.log`. Tables: `F_tables/`.

## Registered verdict

**R1 FAILED.** None of the five independently defined signatures reproduced the dominance pattern at the
registered threshold.

- **Reproduction gate passed:** 19 rPRG rows, maximum |effect difference| 1.4e-16, maximum |p difference|
  8.6e-17.
- **R2 failed:** none of the 4 testable non-overlapping versions passed.
- **R3 passed only as written and is uninformative:** S1, its comparator, is itself ≈ 0.

Per the plan, the paper must state that independently defined signatures do not reproduce the dominance
pattern, and the lead claim is restricted to the package rPRG panel.

Pooled RRRM-2, full versions (mapped genes; effect in log2; two-sided exact p; q in the 10-test primary
family):

| Signature | Genes | P_G (protocol, no flight) | D (dominance) | q (D) |
|---|---|---|---|---|
| S1 Hrvatin 2018, 1-h early-response genes | 255 | −0.001 (p 0.99) | −0.022 (p 0.88) | 0.88 |
| S3 Marsh 2022, ex vivo activation, all CNS cells | 14 | +0.577 (p 0.011) | +0.465 (p 0.031) | 0.14 |
| S4 Marsh 2022, ex vivo activation, myeloid | 18 | +0.366 (p 0.011) | +0.290 (p 0.021) | 0.14 |
| S5 van den Brink 2017, dissociation-affected | 131 | +0.049 (p 0.54) | +0.020 (p 0.82) | 0.88 |
| S6 Wu 2017, IEG list | 110 | +0.125 (p 0.15) | +0.073 (p 0.43) | 0.55 |
| S2 Hrvatin 2018, late-response only (control) | 146 | −0.063 (p 0.46) | −0.064 (p 0.50) | – |

**Non-overlapping versions** (the package's 107 panel genes removed): pooled D is −0.078, −0.161, −0.096
and −0.091 for S1, S4, S5 and S6, and none is significant. S3 is not testable: only 4 genes remain.

**Where a signal appears, it comes from shared genes.** The two Marsh signatures are the only ones with a
nominal signal, and 10 of their 14 and 12 of their 18 mapped genes are in the package rapid-response or
heat-shock panels. Per-region results (descriptive) are in `F_tables/F1_signature_results.csv`. The Marsh
signatures are strongest in hippocampus: D +0.627 (p 0.005) and +0.427 (p 0.010).

**Mapping notes:**
- Hrvatin lists several Ensembl IDs in one cell for 10 early-response genes, including Fos, Dusp1 and Dnajb1.
- Under the registered ID rule these genes did not map; the registered result stands.
- A labelled post-hoc sensitivity check mapped S1 by symbol instead (253 genes;
  `posthoc_S1_symbol_sensitivity.log`). It gives pooled P_G +0.011 (p 0.91) and D −0.008 (p 0.96), so
  the null does not depend on the mapping rule.
- Old histone and chemokine symbols in Marsh and Wu did not match the current annotation; the package rule
  allows no alias substitution.

## What the failure means

The protocol signal does not belong to broad activity-regulated, IEG or dissociation programmes. Three
lists of 110–255 genes defined in other tissues and settings show no mean protocol effect:
- visual cortex after light (Hrvatin 2018);
- a compiled list of 139 IEGs, many immune or vascular (Wu 2017);
- muscle stem cells after collagenase digestion (van den Brink 2017).

**Consequence for the paper:** "immediate-early genes" in the title over-generalises. The supported claim
concerns rapid primary-response genes.

## POST HOC (cannot change the registered verdict): the originating paper's own gene classes

Script: `posthoc_tyssowski.py`. Output: `posthoc_results.json`, `F_tables/PH1_tyssowski_classes.csv`.
The source is Tyssowski et al. 2018 Table S6, archived in `sources/` with its SHA-256 in
`posthoc_tyssowski_classes.json`. It gives 19 rapid PRGs, 116 delayed PRGs and 38 secondary-response genes.
15 of the 19 rapid PRGs are in the package rPRG panel.

| Class (pooled RRRM-2) | Genes | P_G | D | q (3 classes) |
|---|---|---|---|---|
| Tyssowski rapid PRGs | 19 | +0.705 (p 0.011) | +0.491 (p 0.026) | 0.077 |
| Tyssowski delayed PRGs | 116 | +0.005 (p 0.96) | −0.062 (p 0.76) | 0.76 |
| Tyssowski SRGs | 38 | −0.116 (p 0.39) | −0.168 (p 0.30) | 0.44 |
| package rPRG (reference) | 19–20 | +0.889 (p 0.0014) | +0.662 (p 0.0043) | – |

**Interpretation, not established by this analysis:**
- The originating paper's own rapid class gives the same direction as the package panel, and about 80% of
  its size. Its dominance contrast does not survive correction across the three classes (q 0.077).
- The pattern matches Tyssowski's finding that brief activity induces only the rapid wave. It is confined to
  rapid PRGs and absent from delayed PRGs, secondary-response genes and 1-h activity programmes. That is
  consistent with a brief stimulus in the minutes before freezing.
- The data cannot show what the stimulus was.

## Proposed changes to the paper (not applied; awaiting the authors' decision)

1. **Title (15 words, the npj limit):** "Collection protocol outweighs spaceflight effects on rapid brain
   primary-response genes and resists post hoc correction".
2. **Results:** one paragraph reports the registered test and its failure. The post-hoc Tyssowski classes
   follow, labelled as post hoc. A Supplementary Table carries all rows.
3. **Introduction:** the empirical-novelty bullet is narrowed to rapid primary-response genes.
4. **Discussion:** one sentence states that the signature has the temporal profile of brief activity, as an
   interpretation.
5. **Methods:** the registered plan is cited with its hash and time.
