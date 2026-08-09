# Spatial Transcriptomics of the Mouse Brain under Spaceflight and Antioxidant Treatment

Re-analysis of NASA OSDR accessions **OSD-682, OSD-685, OSD-698 and OSD-699**,
centred on a pre-specified CNS-injury / extracellular-vesicle (EV) mRNA target panel.

> **The analysis lives in
> [`OSD-682-685-698-699-brain-spatial/`](OSD-682-685-698-699-brain-spatial) —
> start with [its README](OSD-682-685-698-699-brain-spatial/README.md).**

The four accessions are not four experiments. They are four brain regions of one
NanoString GeoMx Digital Spatial Profiling run (GEO
[GSE239336](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE239336)): CA1,
dentate gyrus, frontal cortex and cerebral cortex, in a 2 × 2 factorial of
spaceflight and the antioxidant BuOE — 48 ROIs, 15,782 targets, n = 3 per cell.

## Headline findings

**The differential-expression statistics deposited with the study cannot be used
as they stand.** The `Adjusted pvalue` column is smaller than its own raw
`Pvalue` for 77% of targets, which no multiple-testing correction can produce,
and it is rounded to two decimals. Read as deposited, 1,317–2,361 targets look
significant per contrast; correct Benjamini-Hochberg control yields **zero**
across the entire experiment.

**Transcriptome-wide, the design has almost no power — and that is quantifiable.**
At n = 3 a target would need a **6–8-fold change** to survive correction across
15,782 targets. The empty genome-wide result reflects the design, not the biology.

**Restricting to the pre-specified target panel is what makes the data informative.**
Testing the 111 measurable panel members instead of 15,782 targets cuts
multiplicity ~140-fold, and eight results then reach BH < 0.05 — **JUNB**,
**MMP12** and **neurogranin**. All have genome-wide adjusted p-values of
0.32–0.87, so none is visible without the restriction.

**Six Ensembl IDs in the source target spreadsheet point to unrelated genes** and
were corrected against Ensembl (Iba1 → NMRAL1, Neurogranin → DIDO1, SBDP → DCX,
SAA → FLCN, plus two apparent transposed-digit typos). One of the three
significant genes is neurogranin, which the supplied ID would have missed.
# Spatial Transcriptomics of the Mouse Brain under Spaceflight and Antioxidant Treatment — OSD-682, OSD-685, OSD-698, OSD-699

A reproducible re-analysis of four NASA Open Science Data Repository (OSDR)
accessions that together form a single NanoString GeoMx Digital Spatial Profiling
experiment, centred on a pre-specified CNS-injury / extracellular-vesicle (EV)
mRNA target panel.

> The four accessions are **not** four experiments. They are four brain regions of
> one GeoMx run: a 2 × 2 factorial of spaceflight and the antioxidant BuOE.
> This repository analyses them as such, re-derives every statistic from the
> normalised counts, and reports what the design can and cannot support.

Sibling project: [`../`](..) covers the bulk RNA-seq datasets OSD-561 and OSD-562.
Both projects share the same CNS/EV target panel and the same dependency pins.

---

## Three findings up front

**1. The differential-expression statistics deposited with this study cannot be used as they stand.**
In the deposited tables the `Adjusted pvalue` column is *smaller* than the raw
`Pvalue` for 77–78% of targets. No multiple-testing correction can lower a
p-value, so the two columns are not consistent with one another. The column is
also rounded to two decimals (at most 101 distinct values, with 1,244–2,261
targets at exactly 0). Read as deposited, 1,317–2,361 targets look significant per
contrast. Recomputing Benjamini-Hochberg from the deposited raw p-values gives
**zero** significant targets across the entire experiment — the smallest adjusted
p-value anywhere is 0.129.

**2. Transcriptome-wide, this design has almost no power, and that is quantifiable rather than speculative.**
With n = 3 per cell the median standard error of a fold change is 0.21–0.25 log2
units. To be detected at 80% power a target would need |log2FC| ≥ 0.67–0.80
nominally, and ≥ 2.6–3.1 (a **6–8-fold change**) to survive correction across
15,782 targets. Brain tissue does not respond to these perturbations at that
magnitude, so the empty genome-wide result reflects the design, not the biology.

**3. Restricting to the pre-specified CNS/EV panel is what makes the data informative.**
Testing the 111 measurable panel members instead of 15,782 targets reduces the
multiple-testing burden ~140-fold, and eight target × region × contrast results
then reach BH < 0.05 — led by **JUNB**, with **MMP12** and **neurogranin (NRGN)**.
Every one of them has a genome-wide adjusted p-value between 0.32 and 0.87, so
none is visible without the panel restriction. This is legitimate only because
the panel was fixed in advance, from the project's target spreadsheet.

<p align="center">
  <img src="results/figures/panels/Figure7_CNS_EV_target_panel.png" alt="CNS/EV target panel results" width="880">
</p>

---

## Datasets

| OSD | GLDS | Brain region | Label | ROIs |
|---|---|---|---|---|
| [OSD-682](https://osdr.nasa.gov/bio/repo/data/studies/OSD-682) | GLDS-613 | Cornu Ammonis 1 (hippocampus) | CA1 | 12 |
| [OSD-685](https://osdr.nasa.gov/bio/repo/data/studies/OSD-685) | GLDS-616 | Dentate Gyrus (hippocampus) | DG | 12 |
| [OSD-698](https://osdr.nasa.gov/bio/repo/data/studies/OSD-698) | GLDS-627 | Frontal Cortex | FCtx | 12 |
| [OSD-699](https://osdr.nasa.gov/bio/repo/data/studies/OSD-699) | GLDS-628 | Cerebral Cortex | Ctx | 12 |
| **combined** | — | all four | — | **48** |

Study: *Spaceflight-Induced Gene Expression Profiles in the Mouse Brain Are
Attenuated by Treatment with the Antioxidant BuOE*. Assay: NanoString GeoMx DSP
spatial transcriptomics, 15,782 gene targets plus one negative-control probe
(`NegProbe-WTX`). Animals: C57BL/6 male mice, 10 weeks at launch, 34-day mission.

### Where the data actually comes from

OSDR hosts **only raw FASTQ** for these four accessions (7–11 GB per study) plus a
raw-read MultiQC report. There are no GeneLab-processed count matrices and no
differential-expression tables, so the approach used for OSD-561/562 cannot be
pointed at them. The processed layer comes from the submitter deposit in NCBI GEO,
[GSE239336](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE239336):
per-region Q3-normalised expression matrices, per-region DSP Analysis Suite
differential-expression exports, the series matrix carrying per-ROI design and QC
metadata, and the 48 raw DCC files. Total download is about 12 MB.

Re-processing the FASTQ was not attempted: the GeoMx DND pipeline also requires
the instrument configuration files, which are not deposited.

## Design

<p align="center">
  <img src="results/figures/panels/Figure1_design_and_QC.png" alt="Design and QC" width="880">
</p>

A 2 × 2 factorial per region, n = 3:

|  | Saline | BuOE |
|---|---|---|
| **Ground Control** | 3 ROIs | 3 ROIs |
| **Space Flight** | 3 ROIs | 3 ROIs |

Each slide carries three tissue sections, and each section is sampled at four
ROIs — one per brain region. So `rep1`–`rep3` are sections (the biological unit),
and the four regions are repeated measures within a section.

### What the design cannot support

`results/tables/csv/design_confounds.csv` records this explicitly.

- **Group is perfectly confounded with slide.** Each experimental group occupies
  its own GeoMx slide (`GC_SAL_2`, `GC_TX_2`, `FLT_SAL_2`, `FLT_TX_2`). Slide,
  batch, staining and scan effects therefore cannot be separated from treatment
  effects by any statistical method.
- **n = 3 per cell**, giving 8 residual degrees of freedom for the 2 × 2 model.
- **Regions are not independent.** All four come from the same sections, so
  cross-region agreement is weaker evidence than agreement across cohorts.

## Methods

### Reading the deposited data honestly

Two conventions in the GEO files are ambiguous and were resolved from the data
rather than assumed.

**Fold-change orientation.** The exports are named `GCvsFLT` and `SALvsBuOE`, but
nothing in the files states which level is the numerator, and guessing would
invert every biological conclusion. `scripts/check_conventions.py` correlates each
deposited `Log2` column against group-mean differences computed from the
Q3 matrix in both orientations. The result is unambiguous: r = 1.0000 with slope
1.0000 in all 12 region × contrast combinations, and the numerator is the
**Flight** or **BuOE** level in every case. The resolved signs are written to
`config/sign_convention.json` and consumed by the rest of the pipeline.

That exact correlation also establishes that the deposited `Log2` column is
precisely the difference of group means on log2 Q3 values, which is why this
pipeline reproduces it to within 1 × 10⁻⁹ (`deposited_vs_derived.csv`).

**Cross-file comparability.** Each region was normalised in its own file, so
combining them needed checking rather than assuming. Every ROI in all four files
has an identical upper quartile (40.02312), so the matrices are directly
comparable (`normalisation_check.csv`).

### Statistics

- **Model.** Per target, a 2 × 2 fixed-effects model on log2 Q3 values with
  effect coding, giving both main effects, the interaction, and the four simple
  contrasts evaluated against the model's residual variance (df = 8). This also
  supplies the Flight-vs-Ground comparison *within BuOE*, which was never
  deposited, and the interaction term that the attenuation claim actually rests on.
- **Multiplicity.** Benjamini-Hochberg at full precision, computed within region ×
  contrast. The deposited adjusted p-values are carried through for provenance
  but never used to call anything.
- **Sensitivity.** Minimum detectable effect at 80% power from the realised
  standard errors, at nominal and genome-wide thresholds.
- **Panel-restricted testing.** BH within the 111 measurable CNS/EV targets. Both
  the panel-restricted and the genome-wide adjusted p-value are reported side by
  side for every target so the effect of the restriction is always visible.
- **Set enrichment.** Whether the panel is more responsive than the rest of the
  transcriptome, as a rank test and as a label-permutation test (500
  permutations, flight/ground shuffled within drug arm) that respects the
  correlation between targets.
- **Attenuation.** See below.

### Mapping a human target list onto mouse data

The project target list is 150 human analytes with human Ensembl gene IDs; the
data are mouse. `scripts/map_targets.py` resolves each ENSG to its current human
symbol and then to its mouse ortholog through Ensembl Compara, keeping the
orthology type. Title-casing human symbols would not have been safe — human
`CXCL8` has no mouse ortholog at all, and several panel members differ in name
between the genomes. **111 of 150** analytes are measurable on this assay;
39 are not (no mouse ortholog, or not on the GeoMx panel). Every mapping decision
is auditable in `config/cns_ev_targets_mapped.csv`.

#### Corrections to the source spreadsheet

Resolving the IDs surfaced six rows whose ENSG points to a gene unrelated to its
label. Four of them are core CNS-injury biomarkers, so analysing them as supplied
would have reported results for the wrong gene under a recognised biomarker name.
Each was corrected to the intended gene and **verified against Ensembl**:

| Analyte | Supplied ENSID | Resolves to | Corrected to | Resolves to | Measured |
|---|---|---|---|---|---|
| Iba1 | `ENSG00000153406` | NMRAL1 | `ENSG00000204472` | **AIF1** | no |
| Neurogranin | `ENSG00000101191` | DIDO1 | `ENSG00000154146` | **NRGN** | yes |
| SBDP (SNTF) | `ENSG00000077279` | DCX | `ENSG00000197694` | **SPTAN1** | yes |
| Serum amyloid alpha (SAA) | `ENSG00000154803` | FLCN | `ENSG00000173432` | **SAA1** | yes |
| APC-CC1 | `ENSG00000135982` | *(no gene)* | `ENSG00000134982` | **APC** | yes |
| PYCARD | `ENSG00000103483` | *(no gene)* | `ENSG00000103490` | **PYCARD** | yes |

The last two look like transposed digits. This is worth fixing at source: one of
the eight panel-restricted findings below is neurogranin, which would have been
missed entirely under the supplied ID.

A further 13 labels differ from their Ensembl symbol but point to the *right*
gene and were left alone — `Amyloid-beta` → APP, `GLT-1` → SLC1A2, `ICE` → CASP1,
`Tau` → MAPT, `TDP43` → TARDBP, `NfL` → NEFL, `ICAM`/`VCAM` → ICAM1/VCAM1, and so
on. All are recorded in the mapping table.

---

## Results

### Global structure: region dominates, treatment barely registers

<p align="center">
  <img src="results/figures/panels/Figure2_global_structure.png" alt="Global structure" width="880">
</p>

PCA over all 48 ROIs separates brain regions first (PC1 = 29.1%, PC2 = 15.2%).
Decomposing each target's variance:

| Factor | Median share of per-target variance |
|---|---|
| Region | 13.0% |
| Section (rep) | 3.2% |
| Spaceflight | 1.3% |
| Treatment | 1.0% |

Section-to-section variability exceeds either experimental factor. That ordering
sets expectations for everything that follows.

### The deposited statistics

<p align="center">
  <img src="results/figures/panels/Figure3_deposited_statistics_audit.png" alt="Deposited statistics audit" width="880">
</p>

Left: 77% of targets fall below the identity line, meaning the "adjusted"
p-value is smaller than the raw one. Middle: the column takes only 101 distinct
values. Right: deposited significance versus correct BH control — every blue bar
is zero.

### Transcriptome-wide differential expression and sensitivity

<p align="center">
  <img src="results/figures/panels/Figure4_DE_landscape_and_sensitivity.png" alt="DE landscape and sensitivity" width="880">
</p>

Across all 28 region × contrast comparisons: **0** targets at BH < 0.05 and
|log2FC| ≥ 1; best adjusted p-value anywhere 0.063. Nominal p < 0.01 counts run
0.80× to 2.56× the null expectation, so there is a modest excess in some
contrasts but nothing that survives correction.

| Region | Median SE (log2) | Detectable log2FC, nominal | Detectable log2FC, genome-wide | As fold change |
|---|---|---|---|---|
| CA1 | 0.250 | 0.80 | 3.07 | 8.4× |
| DG | 0.223 | 0.71 | 2.74 | 6.7× |
| FCtx | 0.215 | 0.69 | 2.64 | 6.2× |
| Ctx | 0.209 | 0.67 | 2.56 | 5.9× |

<p align="center">
  <img src="results/figures/panels/Figure5_volcanoes.png" alt="Volcano plots" width="920">
</p>

Exploratory structure is nonetheless biologically coherent: cold-inducible
`Rbm3` and `Cirbp` up, prostaglandin D synthase `Ptgds` up, metallothionein `Mt1`
up, and the activity-dependent gene `Arc` down in flight. Six targets recur in at
least three regions for Flight-vs-Ground in saline (`Mt1`, `Ptgds`, `Rbm3`,
`Cirbp`, `Arc`, `4930578G10Rik`), all with consistent direction — but a
label-permutation null that preserves the shared-section structure produces the
same recurrence count, so cross-region recurrence here is not evidence beyond
the per-region statistics (`recurrence_null_test.csv`).

### Does BuOE attenuate the spaceflight response?

<p align="center">
  <img src="results/figures/panels/Figure6_BuOE_attenuation.png" alt="BuOE attenuation" width="880">
</p>

The formal per-target test is the Flight × BuOE interaction, and no target reaches
FDR significance for it transcriptome-wide. So the question was also asked as an
aggregate effect-size comparison. Because the sampling variance of every fold
change is known from the model, the mean squared *true* effect can be estimated
without bias as mean(log2FC² − SE²) and compared between arms.

| Region | RMS effect, saline | RMS effect, BuOE | saline − BuOE (95% CI) | bootstrap p | Verdict |
|---|---|---|---|---|---|
| CA1 | 0.101 | 0.081 | +0.0037 (+0.0008, +0.0067) | 0.013 | attenuation (borderline) |
| DG | 0.119 | 0.071 | +0.0091 (+0.0066, +0.0117) | < 0.001 | **attenuation** |
| FCtx | 0.081 | 0.066 | +0.0022 (+0.0000, +0.0044) | 0.047 | attenuation (borderline) |
| Ctx | 0.040 | 0.064 | −0.0025 (−0.0045, −0.0006) | 0.010 | amplification (borderline) |

Directionally this supports the published claim in the hippocampus, most clearly
in dentate gyrus, and not in cerebral cortex. Two caveats are load-bearing:

- **The conclusion depends on the variance assumption.** The BuOE arm carries
  1.09–1.42× the within-arm variance of the saline arm. Using one pooled variance
  for all four cells understates BuOE noise and overstates saline noise, and flips
  the sign in **all four regions**. The primary analysis therefore estimates each
  arm's variance from its own two cells; both versions are reported
  (Figure 6d–e, `attenuation_summary.csv`).
- **Selecting targets on the saline response manufactures attenuation.** Taking
  the top 500 targets by saline p-value gives median |log2FC| ≈ 0.5 in saline
  versus ≈ 0.18 in BuOE — but 500 *random* targets give ≈ 0.18 in both arms
  (Figure 6f). Any analysis that conditions on the saline response will show
  apparent attenuation through regression to the mean alone.

Effect magnitudes are small in absolute terms throughout: an RMS true effect of
0.10–0.12 log2 units is a 7–9% expression change.

### CNS / EV target panel

This is where the data becomes usable. Restricting to the pre-specified panel,
eight results reach BH < 0.05 within the panel:

| Analyte | Gene | Region | Contrast | log2FC | Panel BH | Genome-wide BH |
|---|---|---|---|---|---|---|
| JUNB | *Junb* | FCtx | BuOE vs saline (flight) | +1.238 | 0.0053 | 0.38 |
| JUNB | *Junb* | DG | Flight × BuOE interaction | −0.534 | 0.0122 | 0.53 |
| JUNB | *Junb* | FCtx | Flight vs Ground (BuOE) | +1.039 | 0.0183 | 0.38 |
| JUNB | *Junb* | FCtx | Flight × BuOE interaction | +0.720 | 0.0211 | 0.60 |
| MMP12 | *Mmp12* | CA1 | Flight × BuOE interaction | −0.453 | 0.0315 | 0.87 |
| Neurogranin | *Nrgn* | DG | Flight vs Ground (saline) | +0.640 | 0.0400 | 0.32 |
| JUNB | *Junb* | DG | Flight vs Ground (saline) | +0.571 | 0.0400 | 0.34 |
| MMP12 | *Mmp12* | CA1 | Flight vs Ground (BuOE) | −0.602 | 0.0480 | 0.36 |

Reading these:

- **JUNB** (AP-1 immediate-early transcription factor) is the strongest signal in
  the panel. It rises in flight in dentate gyrus under saline (+0.57), and rises
  markedly in frontal cortex specifically in the BuOE-treated flight animals
  (+1.24 vs saline in flight; +1.04 flight vs ground within BuOE). It carries a
  significant Flight × BuOE interaction in both DG (−0.53) and FCtx (+0.72) — with
  **opposite signs**, i.e. the drug modifies the flight response in opposite
  directions in the two regions.
- **MMP12** (macrophage metalloelastase, matrix and barrier remodelling) falls in
  flight under BuOE in CA1 (−0.60) with a significant interaction (−0.45): a
  hippocampal flight response that the antioxidant suppresses.
- **Neurogranin (NRGN)**, a postsynaptic calmodulin-binding protein used clinically
  as a CSF marker of synaptic injury, rises in flight in dentate gyrus (+0.64).
  This result exists only because the spreadsheet's ENSG for "Neurogranin" was
  corrected from DIDO1 to NRGN.

The panel is also more responsive than the transcriptome as a whole in seven
region × contrast combinations (label-permutation p < 0.05), most strongly for
BuOE-vs-saline in flight in dentate gyrus (p = 0.004) and for flight-vs-ground in
saline in dentate gyrus (p = 0.026).

At the functional-category level, two results survive BH correction across
categories, both in frontal cortex under BuOE in flight: the
**blood-brain-barrier / endothelial** group (mean t = −0.99, q = 0.028) and the
**neuroinflammation / cytokine** group (mean t = −0.58, q = 0.028) both shift
downward, while the neuronal / synaptic group shifts upward (mean t = +1.39,
p = 0.007, q = 0.082). Under flight in saline, the hippocampal pattern runs the
other way: neuronal/synaptic, glial/myelin and injury-biomarker categories all
trend upward in CA1 and DG (Figure 7d).

Panel-restricted attenuation testing is inconclusive in all four regions: with
111 targets rather than 15,782 the bootstrap interval widens beyond a verdict.

**These are exploratory findings from n = 3 with treatment confounded by slide.**
They are the right candidates to carry into a targeted follow-up — JUNB, MMP12
and neurogranin, in dentate gyrus and frontal cortex — not established effects.

---

## Repository layout

```
.
├── README.md
├── Makefile                          # one-command pipeline
├── requirements.txt                  # pinned, matches the OSD-561/562 project
├── config/
│   ├── cns_ev_targets_mapped.csv     # human -> mouse mapping audit trail
│   ├── focus_genes.txt               # generated: measurable mouse targets
│   └── sign_convention.json          # generated: resolved log2 orientation
├── scripts/
│   ├── common.py                     # IO, design, 2x2 model, BH
│   ├── plotting.py                   # figure style and volcano component
│   ├── download_data.py              # GEO processed layer + OSDR ISA metadata
│   ├── check_conventions.py          # resolves the deposited log2 sign empirically
│   ├── map_targets.py                # target panel -> mouse GeoMx targets
│   ├── analyze_qc.py                 # cohort, QC, confound audit
│   ├── analyze_expression.py         # normalisation check, PCA, variance
│   ├── analyze_de.py                 # DE, deposited-stats audit, sensitivity, attenuation
│   ├── analyze_targets.py            # panel-restricted deep analysis
│   ├── analyze_recurrence.py         # cross-region recurrence vs permutation null
│   ├── make_panels.py                # Figures 1-7
│   ├── export_results.py             # CSV + Excel workbook
│   └── dataset_summary.py            # summary md + json
├── data/                             # gitignored, ~12 MB, regenerate on demand
└── results/
    ├── figures/panels/               # Figure1-Figure7 (PDF + PNG)
    ├── figures/                      # individual source figures
    └── tables/
        ├── OSD-682_685_698_699_results.xlsx
        ├── dataset_summary.{md,json}
        └── csv/
```

## Quick start

```bash
make setup     # uv venv on Python 3.12 + pinned requirements (falls back to venv/pip)
make all       # data -> conventions -> analyze -> panels -> tables -> summary
```

Or step by step:

```bash
uv venv --python 3.12 .venv && uv pip install -r requirements.txt
./.venv/bin/python scripts/download_data.py       # ~12 MB from GEO + OSDR
./.venv/bin/python scripts/check_conventions.py
./.venv/bin/python scripts/map_targets.py         # needs network once, then cached
./.venv/bin/python scripts/analyze_qc.py
./.venv/bin/python scripts/analyze_expression.py
./.venv/bin/python scripts/analyze_de.py
./.venv/bin/python scripts/analyze_targets.py
./.venv/bin/python scripts/analyze_recurrence.py
./.venv/bin/python scripts/make_panels.py
./.venv/bin/python scripts/export_results.py
./.venv/bin/python scripts/dataset_summary.py
```

`scripts/download_data.py --with-dcc` also unpacks the 48 raw DCC files.
`scripts/map_targets.py --offline` reuses the cached Ensembl responses.

## Result tables

`results/tables/OSD-682_685_698_699_results.xlsx` collects 27 sheets, opening with
a README sheet that carries the caveats above. The same tables are CSVs under
`results/tables/csv/`.

| Table | Contents |
|---|---|
| `cohort_summary`, `sample_metadata`, `qc_metrics` | 48 ROIs with design, geometry and sequencing QC |
| `design_confounds` | The slide confound and replication structure, stated explicitly |
| `normalisation_check`, `log2FC_orientation` | Evidence that the deposited files were read correctly |
| `deposited_stats_audit` | Deposited adjusted p-values vs correct BH |
| `deposited_vs_derived` | Reproduction of the deposited fold changes (max diff 1e-9) |
| `PCA_*`, `variance_partition` | Global structure |
| `DEG_summary_all`, `sensitivity_MDE` | Counts under each definition, and what was detectable |
| `attenuation_summary` | BuOE attenuation including the variance-assumption sensitivity |
| `target_mapping`, `target_coverage` | Human → mouse mapping with corrections |
| `target_stats`, `target_sig_panelFDR` | Panel results, panel-restricted and genome-wide FDR |
| `target_setenrichment`, `target_categories` | Panel vs transcriptome; category-level effects |
| `recurrent_targets`, `recurrence_null` | Cross-region recurrence vs its permutation null |

Full gene-level DE for all 28 comparisons is written as
`csv/DE_<region>_<contrast>.csv` (~440k rows total). These are regenerable and
gitignored.

## Limitations

- **Secondary re-analysis.** Built from the submitter's Q3-normalised matrices.
  No re-processing from raw reads; the GeoMx DND configuration files needed for
  that are not deposited.
- **Treatment is confounded with slide.** Unavoidable in this deposit, and it
  bounds every causal reading.
- **n = 3, with regions nested in shared sections.** Cross-region agreement is
  not independent replication.
- **Nothing survives transcriptome-wide correction.** All gene-level findings here
  are either panel-restricted or explicitly labelled exploratory.
- **Panel coverage is partial.** 111 of 150 analytes are measurable; notable
  absences include IL1B, IL6, TLR4, TGFB1, VCAM1, PECAM1, CASP1 and AIF1/Iba1.
  Their absence is an assay-content limitation, not a negative result.
- **The panel-restricted analysis is only valid because the panel is
  pre-specified.** Do not extend it to genes chosen after seeing these results.
- **Category assignments are editable judgement calls** (`CATEGORIES` in
  `scripts/map_targets.py`); they affect aggregation and plotting, never which
  targets are tested.

## Reproducibility

Python 3.12.13 with the exact pins in `requirements.txt`, matching the sibling
OSD-561/562 project. Bootstraps and permutations use fixed, region-seeded
generators (`zlib.crc32`, not the salted built-in `hash`), so intervals are
reproducible across runs and independent of loop order. Figures embed TrueType
fonts (`pdf.fonttype 42`).

## Authors

**Isarar Siddique** — [@isararsiddique](https://github.com/isararsiddique) ·
<siddisrar786@gmail.com>

To cite this analysis:

```
Siddique, I. Spatial Transcriptomics of the Mouse Brain under Spaceflight and
Antioxidant Treatment: a re-analysis of OSD-682, OSD-685, OSD-698 and OSD-699.
GitHub repository.
https://github.com/isararsiddique/OSD-561-562-brain-transcriptomics
```

## Data availability and attribution

Inputs are openly available from the NASA Open Science Data Repository
([OSD-682](https://osdr.nasa.gov/bio/repo/data/studies/OSD-682),
[OSD-685](https://osdr.nasa.gov/bio/repo/data/studies/OSD-685),
[OSD-698](https://osdr.nasa.gov/bio/repo/data/studies/OSD-698),
[OSD-699](https://osdr.nasa.gov/bio/repo/data/studies/OSD-699)) and from NCBI GEO
([GSE239336](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE239336)).
Ortholog mapping uses the [Ensembl REST API](https://rest.ensembl.org). No raw
data is redistributed here; `scripts/download_data.py` fetches it on demand.
Please cite OSDR, GEO and the original data contributors when reusing.

## License

Not yet set. The analysis code may be released under a permissive license of your
choice; the underlying OSDR and GEO data remain subject to their own terms.
Full detail, figures, caveats and limitations:
[`OSD-682-685-698-699-brain-spatial/README.md`](OSD-682-685-698-699-brain-spatial/README.md).

## Layout

```
.
├── CNS mRNA targets & Others in EV ... .csv   # source target panel (150 analytes)
└── OSD-682-685-698-699-brain-spatial/         # the analysis
    ├── README.md                              # methods, results, limitations
    ├── Makefile                               # make setup && make all
    ├── scripts/                                # 13 pipeline stages
    ├── config/                                 # target mapping, sign convention
    └── results/                                # figures (PDF + PNG) and tables
```

## Quick start

```bash
cd OSD-682-685-698-699-brain-spatial
make setup     # uv venv on Python 3.12 + pinned requirements
make all       # data -> conventions -> analyze -> panels -> tables -> summary
```

Input data (~12 MB) is downloaded on demand and is not versioned here.

## The OSD-561/562 bulk RNA-seq project

That project previously occupied this repository root. It is preserved on the
remote and is not lost:

| Ref | Contents |
|---|---|
| `backup/osd-561-562-main` branch | full project as it stood, 112 files |
| `osd-561-562-main-backup` tag | same commit, `81f83ef` |

To recover it:

```bash
git checkout backup/osd-561-562-main
```

## Data availability

Inputs are openly available from the NASA Open Science Data Repository
([OSD-682](https://osdr.nasa.gov/bio/repo/data/studies/OSD-682),
[OSD-685](https://osdr.nasa.gov/bio/repo/data/studies/OSD-685),
[OSD-698](https://osdr.nasa.gov/bio/repo/data/studies/OSD-698),
[OSD-699](https://osdr.nasa.gov/bio/repo/data/studies/OSD-699)) and from NCBI GEO
([GSE239336](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE239336)).
Ortholog mapping uses the [Ensembl REST API](https://rest.ensembl.org). No raw
data is redistributed; the download scripts fetch it on demand. Please cite OSDR,
GEO and the original data contributors when reusing.

## Authors

**Isarar Siddique** — [@isararsiddique](https://github.com/isararsiddique) ·
<siddisrar786@gmail.com>

## License

Not yet set. The analysis code may be released under a permissive license of your
choice; the underlying OSDR and GEO data remain subject to their own terms.
