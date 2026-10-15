# Parent diagnostic, run after Worker G finished (15 Oct 2026)

This is not part of Worker G's independent result. G's frozen numbers (`../G_results.json`) and report
(`../G_REPORT.md`) are unchanged. The parent read the package code to find why G's independent numbers
differed slightly. It then re-ran G's own code (copied to `code/`) with the two specification details found
there. Both details are missing from the paper's Methods text.

## The two missing details

1. **ERCC spike-ins are removed before filtering and normalisation.** The count matrices hold 92 ERCC rows.
   The package `normalize()` keeps only rows matching `^ENSMUSG\d+$`. G kept them, so a few ERCC rows that
   pass the count filter entered the median-ratio size factors.
   - Edit (in `code/g_common.py`, `load_deposit`): one line that drops non-ENSMUSG rows.
2. **The noise ceiling uses pooled within-cell residual variance.** The cells are cohort × age, with
   df = n − 4 in RRRM-2. G used the residual of the additive age + group model (df = n − 3). RR-23 has a
   single age, so the two definitions coincide there.
   - Edit (in `code/g_part2.py`, `coef`): within-cell variance when age varies.

Outputs were redirected to this folder.

## Result: G's independent code then reproduces the paper

**Part 1 (dominance test and rPRG contrasts), ERCC removed** (`G_part13_ercc.json`, `part13_ercc.log`):
- Across 34 rows (24 dominance tests and the rPRG contrasts), the maximum |effect difference| from the
  package is 1.5e-14 and the maximum |p difference| is 1.1e-16.
- The three 2-dp values G flagged are therefore correct as printed: CB F_I −0.21, HPC F_I 0.16 and
  HPC P_F 1.00.

**Part 2 (noise ceiling).** ERCC-only results are in `part2_ercc_only.log`; ERCC plus within-cell variance
in `G_part2_ercc_plus_cellvar.json`.

| Setting | Observed cosine | Expected cosine | S/N_tr, whole genome | S/N_tr, sensor |
|---|---|---|---|---|
| RRRM-2 cerebellum | 0.226 (paper 0.23) | 0.211 (0.21) | 0.383 (0.38) | 2.81 (2.8) |
| RRRM-2 hippocampus | 0.521 (0.52) | 0.502 (0.50) | 1.467 (1.47) | 6.60 (6.6) |
| RR-23 cerebellum | 0.165 (0.17) | 0.148 (0.15) | 0.236 (0.24) | 4.94 (not printed) |

- ERCC removal alone fixes all the observed cosines, including RR-23 (0.195 to 0.165), and the RR-23
  expected cosine and S/N.
- The within-cell variance fixes the RRRM-2 hippocampus S/N. With ERCC removal alone it is 1.539; with
  both details it is 1.467.

**Part 3** (published signatures) already matched under F's mapping rule.

## Consequence for the paper

Add both details to the Methods so the analysis can be reproduced from the text alone. Applied on
15 Oct 2026 in `src/sections/methods.tex`.
