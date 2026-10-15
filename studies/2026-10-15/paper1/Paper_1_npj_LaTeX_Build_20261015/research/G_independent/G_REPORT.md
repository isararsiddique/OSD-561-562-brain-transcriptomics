# G — Independent re-implementation: dominance test, noise ceiling, published signatures

> **Parent addendum (15 Oct 2026, after G finished).** G's results below are unchanged. Both unexplained
> differences trace to details the Methods text omitted:
> - ERCC spike-in rows are removed before normalisation;
> - the noise ceiling uses pooled within-cell (cohort × age) residual variance.
>
> With those two details, G's own code reproduces Part 1 to 1.5e-14 and every Part 2 value at the printed
> precision. The three flagged 2-dp values are correct as printed. See `PARENT_DIAG_ercc/README.md`. The
> Methods now state both details.

Worker G, 15 Oct 2026. Code: `G_code/`. Frozen primary numbers: `G_results.json`. They were written before
any paper value was read. Diagnostics run after the comparison are kept separately and change nothing:
`G_diag_part1.txt`, `G_diag_part2.json`.

## Verdicts

- **Part 1 (dominance test): MATCH at the paper's reporting precision, with small residual differences.**
  - Every Results-text value for the dominance test matches, and the same 9 of 24 tests pass BH.
  - Supplementary Table 2: 20/24 effects and 19/24 exact p values match at the printed precision.
  - Where they differ, the gap is at most 0.0007 log2 and at most 2 enumerated assignments.
  - The four rPRG contrasts differ by 0.0005–0.0025 log2. Three of these differences change the
    2-decimal value printed in the text (CB F_I, HPC F_I, HPC P_F).
- **Part 2 (noise ceiling): PARTIAL.**
  - The observed whole-genome cosines MATCH for both RRRM-2 regions. RR-23 DIFFERS (0.20 vs 0.17).
  - The expected cosines and every S/N_tr DIFFER. My S/N is systematically higher (by 4–25%), so my
    expected cosines are 0.007–0.028 higher.
  - The qualitative claim holds in my numbers: observed ≈ expected in all three settings.
- **Part 3 (published signatures): MATCH,** but only under F's mapping rule. F counts an Ensembl cell that
  lists several IDs as unmapped.
  - With that rule, all 24 numbers agree except S1 P_G (−0.002 vs −0.001).
  - My pre-registered primary rule kept every listed ID. It gives slightly different S1/S2 values; the
    conclusion is the same.

## Part 1 — dominance (GC ISS-T − FLT LAR), 24 tests, BH over 24

Assignments enumerated: 200 / 400 / 700 (MATCH with Methods and SI Table 2).

| Region | Panel | G effect | Paper | G p | Paper p | G q | Paper q | Class |
|---|---|---|---|---|---|---|---|---|
| CB | rPRG | 0.5036 | 0.503 | 0.0100 | 0.0100 | 0.0343 | 0.0343 | effect DIFFERS +0.0006; p, q MATCH |
| CB | dPRG | 0.136 | 0.136 | 0.4450 | 0.4400 | 0.4855 | 0.4855 | p DIFFERS (1/200) |
| CB | SRG | −0.443 | −0.443 | 0.1600 | 0.1600 | 0.2475 | 0.2475 | MATCH |
| CB | HSP | 0.2107 | 0.210 | 0.2000 | 0.1900 | 0.2824 | 0.2682 | DIFFERS (+0.0007; 2/200) |
| CB | Glia_Complement | −0.3132 | −0.314 | 0.1550 | 0.1550 | 0.2475 | 0.2475 | effect DIFFERS 0.0008 |
| CB | BBB_Endothelium | −0.048 | −0.048 | 0.8650 | 0.8650 | 0.9026 | 0.9026 | MATCH |
| CB | Circadian | 0.091 | 0.091 | 0.4250 | 0.4250 | 0.4855 | 0.4855 | MATCH |
| CB | Sensor31 | 0.396 | 0.396 | 0.0250 | 0.0250 | 0.0600 | 0.0600 | MATCH |
| HPC | rPRG | 0.837 | 0.837 | 0.0150 | 0.0150 | 0.0419 | 0.0419 | MATCH |
| HPC | dPRG | 0.133 | 0.133 | 0.4450 | 0.4450 | 0.4855 | 0.4855 | MATCH |
| HPC | SRG | 0.385 | 0.385 | 0.0100 | 0.0100 | 0.0343 | 0.0343 | MATCH |
| HPC | HSP | 0.351 | 0.351 | 0.0100 | 0.0100 | 0.0343 | 0.0343 | MATCH |
| HPC | Glia_Complement | −0.223 | −0.223 | 0.1650 | 0.1650 | 0.2475 | 0.2475 | MATCH |
| HPC | BBB_Endothelium | −0.359 | −0.359 | 0.0550 | 0.0550 | 0.1100 | 0.1100 | MATCH |
| HPC | Circadian | 0.247 | 0.247 | 0.0100 | 0.0100 | 0.0343 | 0.0343 | MATCH |
| HPC | Sensor31 | 0.6645 | 0.664 | 0.0100 | 0.0100 | 0.0343 | 0.0343 | MATCH within rounding (0.664502) |
| Pooled | rPRG | 0.662 | 0.662 | 0.0043 | 0.0043 | 0.0343 | 0.0343 | MATCH |
| Pooled | dPRG | 0.132 | 0.132 | 0.3857 | 0.3843 | 0.4855 | 0.4854 | p DIFFERS (1/700) |
| Pooled | SRG | −0.022 | −0.022 | 0.9157 | 0.9143 | 0.9157 | 0.9143 | p DIFFERS (1/700) |
| Pooled | HSP | 0.281 | 0.281 | 0.0286 | 0.0286 | 0.0623 | 0.0623 | MATCH |
| Pooled | Glia_Complement | −0.262 | −0.262 | 0.1129 | 0.1129 | 0.2084 | 0.2084 | MATCH |
| Pooled | BBB_Endothelium | −0.211 | −0.211 | 0.3143 | 0.3114 | 0.4190 | 0.4152 | p DIFFERS (2/700) |
| Pooled | Circadian | 0.174 | 0.174 | 0.0157 | 0.0157 | 0.0419 | 0.0419 | MATCH |
| Pooled | Sensor31 | 0.522 | 0.522 | 0.0071 | 0.0071 | 0.0343 | 0.0343 | MATCH |

**Results text:**
- Dominance 0.50 (p 0.010, q 0.034), 0.84 (p 0.015, q 0.042) and pooled 0.66 (p 0.0043, q 0.034): MATCH.
- "Nine of the 24 passed BH": MATCH. Same nine rows: CB rPRG; HPC rPRG, SRG, HSP, Circadian, Sensor31;
  pooled rPRG, Circadian, Sensor31.
- dPRG shows no dominance: MATCH.

### Four RRRM-2 rPRG contrasts

Paper values are from SI Table 3 (3 dp) and the Results text (2 dp).

| Region | Contrast | G effect | G p (assignments) | Paper (SI T3) | Paper text | Class |
|---|---|---|---|---|---|---|
| CB | P_G | 0.7035 | 0.0050 (400) | 0.704 | 0.70, p 0.005 | DIFFERS 0.0005 at 3 dp; text MATCH |
| CB | F_L | 0.1874 | 0.380 (200) | 0.188 | 0.19, p 0.38 | DIFFERS 0.0006; text MATCH |
| CB | F_I | −0.2037 | 0.359 (1960) | −0.206 | −0.21, p 0.36 | DIFFERS 0.0023; text effect DIFFERS (−0.20); p MATCH |
| CB | P_F | 0.3345 | 0.207 (735) | 0.332 | 0.33, p 0.21 | DIFFERS 0.0025; text MATCH |
| HPC | P_G | 1.1184 | 0.010 (400) | 1.118 | 1.12, p 0.010 | MATCH |
| HPC | F_L | 0.2817 | 0.060 (400) | 0.281 | 0.28, p 0.060 | DIFFERS 0.0007; text MATCH |
| HPC | F_I | 0.1653 | 0.491 (1960) | 0.164 | 0.16, p 0.49 | DIFFERS 0.0013; text effect DIFFERS (0.17) |
| HPC | P_F | 1.0054 | 0.00051 (1960) | 1.004 | 1.00, p 0.0005 | DIFFERS 0.0014; text effect DIFFERS (1.01) |
| Pooled | P_G | 0.888 | 0.0014 (700) | – | 0.89, p 0.0014 | MATCH |

**Diagnosis:**
- The p-value differences are near-tie flips caused by the sub-0.001 effect differences. They do not
  indicate a different enumeration: the assignment counts match.
- The largest gaps are in contrasts that include flight ISS-T animals (F_I, P_F).
- A post-comparison diagnostic tested four alternative preprocessing choices
  (`G_code/g_diag_part1.py`, `G_diag_part1.txt`). None closed the gap:
  - size factors rescaled to geometric mean 1;
  - size factors computed before the count filter;
  - both of the above together;
  - normalisation restricted to each contrast's animals (this one moves the effects further away).
- Panel-to-Ensembl mapping is one-to-one for every panel symbol. Arntl is absent from
  gene_annotation.csv.gz, Egr2 fails the filter in CB, and Fosl1 and Serpinb2 fail it in both regions.
- The residual difference is therefore unexplained from the Methods text. It is consistent with an
  undocumented detail in the authors' loader. My leading suspect is a difference in how flight ISS-T
  profiles are assembled or normalised, but I cannot test that without the authors' code.
- It changes no conclusion. It does change three 2-dp values in the text (CB F_I, HPC F_I, HPC P_F); the
  authors should confirm those against their tables.

## Part 2 — identical-loading noise ceiling (ground-only holdouts: 81 / 81 / 200)

| Quantity | CB: G | CB: paper | HPC: G | HPC: paper | RR-23: G | RR-23: paper | Class |
|---|---|---|---|---|---|---|---|
| Observed mean cosine, whole genome | 0.232 | 0.23 | 0.524 | 0.52 | 0.195 | 0.17 | CB, HPC MATCH; RR-23 DIFFERS +0.025 |
| Expected cosine, whole genome | 0.217 | 0.21 | 0.517 | 0.50 | 0.178 | 0.15 | DIFFERS +0.007 / +0.017 / +0.028 |
| S/N_tr, whole genome | 0.396 | 0.38 | 1.560 | 1.47 | 0.295 | 0.24 | DIFFERS (+4% / +6% / +23%) |
| S/N_tr, Sensor31 | 2.86 | 2.8 | 7.21 | 6.6 | 4.95 | – | CB MATCH at 1 dp; HPC DIFFERS +0.6 |
| Observed mean cosine, sensor | 0.688 | ≈0.69 (Fig. 4b, visual) | 0.843 | ≈0.85 (visual) | 0.813 | ≈0.81 (visual) | NOT RUN numerically (no printed value); visual agreement |
| Expected cosine, sensor | 0.661 | ≈0.66 (visual) | 0.829 | ≈0.82 (visual) | 0.780 | not available (caption) | NOT RUN numerically; visual agreement |

My 400 Gaussian simulations reproduce my own eq. 3 values within 0.001 for the whole genome and within
0.015 for the sensor (`G_part2.json`).

**Diagnosis:**
- The paper's own numbers are internally consistent with eq. 3. With r = S/N_tr and N_te/N_tr = 2
  (RRRM-2) or 1.82 (RR-23), r = 0.38, 1.47 and 0.24 give 0.21, 0.50 and 0.15. The disagreement therefore
  lies in the estimate of S relative to N, not in the formula.
- The Methods say only that S is "estimated without bias from the ground animals" and that N is "derived
  from per-gene residual variances and the group sizes of each split". The estimator, the residual degrees
  of freedom and the gene set are all unstated.
- A post-comparison variant estimated S and N inside each split's training animals
  (`G_diag_part2.json`). It gives a CB median S/N of 0.378, which matches 0.38. It does not match HPC
  (1.57) or RR-23 (0.27), so it is not the paper's estimator either.
- The RR-23 observed-cosine gap (0.195 vs 0.17) is not explained by fixing the gene filter on the whole
  deposit (0.199). It may reflect different RR-23 preprocessing. Candidates are the handling of technical
  libraries (GC1 and V5 each have two) or a different count source. That is untested.
- Under every variant I ran, the observed cosine stays within about 0.02 of the identical-loading ceiling.
  The paper's inference (low cosine does not imply non-transfer) is reproduced.

## Part 3 — published signatures, pooled RRRM-2, full versions

| Signature | G genes | G P_G (p) | G D (p) | F genes | F P_G (p) | F D (p) | Class |
|---|---|---|---|---|---|---|---|
| S1 Hrvatin ERG (single-ID rule) | 255 | −0.0016 (0.99) | −0.022 (0.88) | 255 | −0.001 (0.99) | −0.022 (0.88) | P_G DIFFERS ≤0.001; rest MATCH |
| S1 (G primary: all listed IDs) | 268 | +0.004 (0.96) | −0.015 (0.93) | – | – | – | mapping-rule difference; same null |
| S2 Hrvatin LRG only (single-ID rule) | 146 | −0.063 (0.46) | −0.064 (0.50) | 146 | −0.063 (0.46) | −0.064 (0.50) | MATCH |
| S2 (G primary: all listed IDs) | 150 | −0.067 (0.41) | −0.070 (0.46) | – | – | – | mapping-rule difference |
| S3 Marsh, all CNS | 14 | +0.577 (0.011) | +0.465 (0.031) | 14 | +0.577 (0.011) | +0.465 (0.031) | MATCH |
| S4 Marsh, myeloid | 18 | +0.366 (0.011) | +0.290 (0.021) | 18 | +0.366 (0.011) | +0.290 (0.021) | MATCH |
| S5 van den Brink | 131 | +0.049 (0.54) | +0.020 (0.82) | 131 | +0.049 (0.54) | +0.020 (0.82) | MATCH |
| S6 Wu IEG | 110 | +0.125 (0.15) | +0.073 (0.43) | 110 | +0.125 (0.15) | +0.073 (0.43) | MATCH |

F's q values (10-test family) are NOT RUN. The brief does not define the family, and BH over the ten
full-version P_G and D tests does not give 0.14 for S3/S4, so F's family evidently includes other tests.

## Ambiguities in the Methods text and the choices I made

1. **Normalisation scope for subgroup contrasts.** Filter and size factors use the whole deposit (all
   flight and ground animals); each contrast is then fitted on its own animals. The subset-normalised
   alternative fits worse.
2. **Size factors.** DESeq-style median of ratios over filtered genes that are positive in every sample.
   Not rescaled to geometric mean 1.
3. **Pooled gene set.** The intersection of genes passing the filter in each deposit (17,100). Pooled
   nuisance terms are intercept, age indicator and tissue indicator, with no age × tissue term.
4. **Exact p.** Two-sided |null| ≥ |obs| with a 1e-12 tolerance, over the full enumeration (which
   includes the observed assignment). Pooled relabelling is per animal, using the Source Name across
   regions.
5. **Panel mapping.** Symbol → all Ensembl IDs carrying that symbol in gene_annotation.csv.gz (all
   one-to-one here). The panel mean is over genes passing the filter.
6. **Holdout preprocessing.** Training animals set the gene filter (mean ≥ 10), the positive-in-all gene
   set and the geometric-mean reference. Test size factors are the median of test count / training
   reference in linear space, so a zero count gives ratio 0 instead of −∞.
7. **Cosine definition.** Uncentred cosine of age-adjusted per-gene OLS coefficients. The whole genome is
   all training-filtered genes. The sensor version uses Sensor31 genes present after the filter.
8. **S, N_tr, N_te.**
   - These come from one fit on all ground animals of the deposit, preprocessed from those animals only.
   - σ²_j = RSS/(n−3) for RRRM-2 and RSS/(n−2) for RR-23.
   - S = Σδ̂² − Σσ²/c_all, N = Σσ²/c, where c = g̃ᵀg̃ (3, 2 and 1 for RRRM-2; 4.24, 2.73 and 1.5 for
     RR-23).
9. **RR-23 setting.** Label 1 is HGC, 0 is VGC, taken from split_animal_roles.csv and checked against the
   metadata. There is no age covariate (single age). Technical libraries are summed per Source Name.
10. **Multi-ID Ensembl cells in S1/S2.** My primary rule kept every listed ID. A single-ID-only rule was
    run as a sensitivity before the comparison, and that is the variant that matches F.
11. **Exact symbol mapping (S3–S6).** No aliases. A symbol carried by several Ensembl IDs contributes all
    of them.

## Independence statement and files opened

I did not open, read, grep, import or execute any .py or .ipynb file under PKG or BUILD, other than my
own code in `G_code/`. I read no result table before my numbers were frozen. The only paper values I used
were those printed in the paper text, Supplementary Table 2 and 3, Fig. 4 and, after Part 3 was final,
F_REPORT.md.

**Files opened** (PKG = `.../Paper1_Strengthened_Review_20261015`, BUILD = `.../Paper_1_npj_LaTeX_Build_20261015`):
- BUILD/briefs/G_independent_brief.md
- BUILD/src/sections/methods.tex (before analysis)
- PKG/inputs/OSD-561_counts.csv.gz, OSD-561_metadata.csv.gz, OSD-562_counts.csv.gz,
  OSD-562_metadata.csv.gz, OSD-525_counts.csv.gz, OSD-525_metadata.csv.gz, gene_annotation.csv.gz
- PKG/results/panel_definitions.csv
- PKG/extensions/ground_holdout/split_animal_roles.csv
- BUILD/research/F_signatures/signatures.json (gene lists and mapping fields)
- After the numbers were frozen:
  - BUILD/src/sections/results.tex (lines 22–53 and 80–100);
  - BUILD/output/SUPPLEMENTARY_INFORMATION.pdf (text extracted: Supplementary Tables 2–3);
  - BUILD/output/PAPER_1_npj_Microgravity.pdf (text extracted; page 6 rendered to read Fig. 4b);
  - BUILD/research/F_signatures/F_REPORT.md (after Part 3 was final).
- Directory listings of PKG/inputs and BUILD/src, plus a header-line column count of the three counts
  files.

Environment notes:
- MCP servers fusion (failed to start), cloudflare-api and github (awaiting authorisation) were
  unavailable. None was needed.
- No network calls and no git commands were made. Nothing was written outside this folder, apart from
  temporary text and PNG files in the session scratch directory.
