# Spaceflight Transcriptomics of the Mouse Brain — OSD-561 & OSD-562

A reproducible, full-data re-analysis of two NASA Open Science Data Repository
(OSDR) bulk RNA-seq datasets from the RRRM-2 mission, focused on the cerebellum
and hippocampus, with an extracellular-vesicle (EV) / neuroinflammation gene panel.

> Every unique comparison in the differential-expression tables is analysed and
> categorised (Spaceflight / Age / Environment / Confounded) — not just the
> spaceflight contrasts — and every downloaded table is used.

<p align="center">
  <img src="results/graphical_abstract/graphical_abstract.png" alt="Graphical abstract" width="720">
</p>

---

## Overview

This repository takes the GeneLab-processed bulk RNA-seq tables for two brain
regions flown on the ISS, and produces a complete, honest analysis:

- **Global structure** — PCA on variance-stabilised counts (all samples, all genes).
- **Differential expression** — every unique group comparison in the DE tables,
  categorised by which experimental factor differs.
- **Focus-gene panel** — an EV-biogenesis + neuroinflammation gene set tracked
  across all clean contrasts.
- **Quality control** — sequencing depth, genes detected, and count-table
  concordance, built from the raw (RSEM) and normalized count tables.

All figures are journal-styled (Arial, embedded fonts, no overlapping labels) and
composed into eight multi-panel figures. All results are exported to CSV and a
single Excel workbook.

## Key result

The dominant axis of transcriptional change is **not** the Flight-vs-Ground
(Spaceflight) contrast. When all comparisons are analysed, the **On-ISS vs
On-Earth (Environment)** axis drives far more differential expression, while the
clean spaceflight effect is small.

| Region | Spaceflight (Flight vs Ground) | Environment (On ISS vs On Earth) |
|---|---|---|
| Cerebellum (OSD-561) | 2–6 DEGs | **41** (36 down, 5 up) |
| Hippocampus (OSD-562) | 1–2 DEGs | **42** (22 up, 20 down) |

The largest gene counts come from **confounded** contrasts (up to 1,081 DEGs in
cerebellum) that mix two or more factors; these are reported for completeness but
are flagged as not cleanly interpretable. Thresholds: FDR < 0.05 and
|log2 fold change| ≥ 1.

<p align="center">
  <img src="results/figures/panels/Figure2_DEG_landscape.png" alt="DEG landscape across all contrasts" width="820">
</p>

## Datasets

| OSD | GLDS | Region | Samples | Genes measured |
|---|---|---|---|---|
| [OSD-561](https://osdr.nasa.gov/bio/repo/data/studies/OSD-561) | GLDS-556 | Cerebellum | 26 | 26,005 |
| [OSD-562](https://osdr.nasa.gov/bio/repo/data/studies/OSD-562) | GLDS-557 | Hippocampus | 27 | 26,067 |
| Combined | — | Cb + HPC | 53 | 24,769 shared |

Design factors per sample: **age** (12 vs 29 week), **treatment**
(Ground Control vs Space Flight), and **environment / collection condition**
(On Earth vs On ISS). Cohort composition (combined): 29 Space Flight / 24 Ground
Control, 27 old / 26 young, 30 ISS-terminal / 23 live-animal-return.

Only the small GeneLab-processed tables are used (a few MB each) — **not** the
raw FASTQ/BAM (hundreds of GB) and **not** the spatial-transcriptomics layer.

## Repository structure

```
.
├── README.md
├── Makefile                      # one-command pipeline
├── requirements.txt              # pinned dependencies (Python 3.12)
├── config/
│   └── focus_genes.txt           # EV / neuroinflammation gene panel (editable)
├── scripts/
│   ├── download_data.py          # fetch GeneLab tables from NASA OSDR
│   ├── analyze.py                # PCA, spaceflight volcanoes, focus heatmaps, contrast helpers
│   ├── analyze_full.py           # DEG landscape, all clean-contrast volcanoes, focus log2FC, QC
│   ├── export_results.py         # all result tables -> CSV + Excel
│   ├── make_panels.py            # compose multi-panel figures
│   └── dataset_summary.py        # dataset overview (md + json)
├── data/                         # gitignored — regenerate with download_data.py (~196 MB)
│   ├── OSD-561/
│   └── OSD-562/
└── results/
    ├── figures/
    │   ├── panels/               # Figure1–Figure8 (PDF + PNG) — the main deliverables
    │   └── *.pdf / *.png         # individual source figures
    ├── tables/
    │   ├── OSD-561_562_results.xlsx   # consolidated workbook (all summary sheets)
    │   ├── all_contrasts_DEG_summary.csv
    │   ├── dataset_summary.{md,json}
    │   └── csv/                  # per-table CSVs (large DE_*.csv are gitignored)
    └── graphical_abstract/
        ├── graphical_abstract.svg
        └── graphical_abstract.png
```

## Requirements

- Python 3.12 (tested on 3.12.13)
- Dependencies pinned in [`requirements.txt`](requirements.txt): pandas, numpy,
  matplotlib, seaborn, scikit-learn, scipy, adjustText, openpyxl.

## Quick start

```bash
# 1. Environment + dependencies
make setup                      # creates .venv and installs requirements.txt

# 2. Download the GeneLab tables from NASA OSDR (~196 MB)
make data

# 3. Run everything: analysis -> panels -> tables -> summary
make all
```

Prefer to run the steps yourself:

```bash
python3 -m venv .venv
./.venv/bin/pip install -r requirements.txt
./.venv/bin/python scripts/download_data.py
./.venv/bin/python scripts/analyze.py
./.venv/bin/python scripts/analyze_full.py
./.venv/bin/python scripts/make_panels.py
./.venv/bin/python scripts/export_results.py
./.venv/bin/python scripts/dataset_summary.py
```

## Figures

All panels live in [`results/figures/panels/`](results/figures/panels) as PDF
(print) and PNG (preview). Individual source figures are in
[`results/figures/`](results/figures).

### Figure 1 — Principal component analysis

<p align="center">
  <img src="results/figures/panels/Figure1_PCA.png" alt="Figure 1: PCA" width="80%">
</p>

PCA on variance-stabilised counts. Combined datasets separate first by brain
region (PC1 = 49.7%); within each region, samples separate by collection
condition more than by flight treatment.

### Figure 2 — Differential-expression landscape (all contrasts)

<p align="center">
  <img src="results/figures/panels/Figure2_DEG_landscape.png" alt="Figure 2: DEG landscape" width="100%">
</p>

DEG counts for all 31 unique contrasts per region (symlog scale, coloured by
category). The Environment (On ISS vs On Earth) and confounded contrasts carry
almost all of the signal; clean Flight-vs-Ground effects are small.

### Figure 3 — Volcano plots: Spaceflight (Flight vs Ground)

<p align="center">
  <img src="results/figures/panels/Figure3_Volcanoes_Spaceflight.png" alt="Figure 3: Spaceflight volcanoes" width="92%">
</p>

All five age- and condition-matched Flight-vs-Ground comparisons.

### Figure 4 — Volcano plots: Age (Old vs Young)

<p align="center">
  <img src="results/figures/panels/Figure4_Volcanoes_Age.png" alt="Figure 4: Age volcanoes" width="92%">
</p>

All five matched 29-week-vs-12-week comparisons.

### Figure 5 — Volcano plots: Environment (On ISS vs On Earth)

<p align="center">
  <img src="results/figures/panels/Figure5_Volcanoes_Environment.png" alt="Figure 5: Environment volcanoes" width="72%">
</p>

All four matched On-ISS-vs-On-Earth comparisons — the largest clean effects.

### Figure 6 — Focus-gene log2 fold-change matrix

<p align="center">
  <img src="results/figures/panels/Figure6_FocusGene_log2FC.png" alt="Figure 6: Focus-gene log2FC" width="100%">
</p>

EV / neuroinflammation focus panel across every clean contrast, per region.

### Figure 7 — Focus-gene heat maps

<p align="center">
  <img src="results/figures/panels/Figure7_FocusGene_heatmaps.png" alt="Figure 7: Focus-gene heatmaps" width="82%">
</p>

Focus-gene VST z-scores for the combined cohort and each region.

### Figure 8 — Quality control

<p align="center">
  <img src="results/figures/panels/Figure8_QC.png" alt="Figure 8: QC" width="100%">
</p>

Library size and genes detected per sample (from RSEM raw counts), and
VST-vs-Normalized count concordance (r ≈ 0.97).

## Result tables

`results/tables/OSD-561_562_results.xlsx` collects every summary table; the same
tables are also written as CSVs under `results/tables/csv/`.

| Table | Description |
|---|---|
| `cohort_summary` | Sample and gene counts per study and combined |
| `sample_metadata` | Every sample with region, treatment, age, environment, animal ID |
| `PCA_coordinates` / `PCA_variance` | PC1/PC2 per sample; variance explained |
| `DEG_summary_all` | DEG counts for all 31 unique contrasts, tagged by category |
| `significant_DEGs_all` | Gene-level significant DEGs across every contrast |
| `focus_gene_stats_all` | log2FC / p / adj-p for the focus panel across every contrast |
| `focus_log2FC_matrix` | Focus genes × clean single-factor contrasts |

Full gene-level DE for every contrast is written to
`results/tables/csv/DE_all_contrasts_long.csv` plus one `DE_<study>_<contrast>.csv`
per comparison. These are large (~200 MB total) and regenerable, so they are
gitignored — run `export_results.py` to rebuild them.

## What "all the data" means here

| Data | Used |
|---|---|
| All 53 samples | Yes |
| All measured genes (~26k per study; 24,769 shared) | Yes |
| All 62 DE columns → 31 unique comparisons (reverse duplicates collapsed) | Yes |
| VST, differential_expression, RSEM (raw), Normalized count tables | Yes |
| SampleTable / contrasts files | Information reconstructed from column and sample names |
| Raw FASTQ / BAM (hundreds of GB) | No — GeneLab-processed tables used instead |
| Spatial-transcriptomics layer | No |

## Analysis details

- **Counts.** PCA and heatmaps use variance-stabilised (VST) counts on the top
  variable genes. QC concordance compares VST against log2(Normalized + 1).
- **Contrasts.** GeneLab's precomputed DESeq2 statistics are read directly. Each
  `Log2fc_(A)v(B)` column is parsed; reverse-direction duplicates (A-vs-B and
  B-vs-A) are collapsed. A comparison is *clean* when exactly one design factor
  differs (Spaceflight, Age, or Environment) and *confounded* otherwise.
- **Orientation.** Clean contrasts are oriented consistently: Space Flight − Ground
  Control, 29 week − 12 week, On ISS − On Earth (positive = up in the first term).
- **DEG calling.** FDR < 0.05 and |log2FC| ≥ 1.

## Limitations

- This is a **secondary re-analysis** of GeneLab's precomputed DESeq2 results; no
  re-alignment or re-modelling from raw reads is performed.
- `config/focus_genes.txt` currently holds a **placeholder** EV/neuroinflammation
  panel (Cd9, Cd63, Cd81, Pdcd6ip, Tsg101, Nfkb1, Il1b, Tnf, Gfap, Aif1). Replace
  it with your finalised list and re-run — the focus figures regenerate automatically.
- **Confounded** contrasts are included for completeness but are not directly
  interpretable because more than one factor differs.
- The "On ISS" vs "On Earth" factor reflects the collection/housing condition as
  encoded by GeneLab; its biological interpretation should be made in context.

## Reproducibility

Pinned to Python 3.12.13 with the exact package versions in `requirements.txt`.
Figures use fixed random seeds where relevant and embedded fonts (matplotlib
`pdf.fonttype 42`).

## Data availability and attribution

Input data are from the NASA Open Science Data Repository (OSDR / GeneLab):
[OSD-561](https://osdr.nasa.gov/bio/repo/data/studies/OSD-561) and
[OSD-562](https://osdr.nasa.gov/bio/repo/data/studies/OSD-562). OSDR data are
openly available; please cite the repository and the original data contributors
when reusing. This repository redistributes **no** raw data — it downloads the
processed tables on demand via `scripts/download_data.py`.

## License

No license is set yet. The analysis code in this repository may be released under
a permissive license of your choice (e.g. MIT); the underlying OSDR data remain
subject to NASA's open-data terms. Add a `LICENSE` file before publishing.
