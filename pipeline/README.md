# Portable CNS transcriptomics pipeline

`run_pipeline.py` is the single entry point for the complete deposited-count study. It runs 62 ordered stages covering the corrected 151-target panel, all six original datasets, independent cohorts, single-cell and spatial context, sensitivity analyses, figures and independent checks. Statistical models are fitted to whole-gene counts before extracting the target panel. No external review application or model service is needed.

## Run the complete analysis

From the folder containing `run_pipeline.py`, with Python and R dependencies installed:

```sh
python run_pipeline.py run --inputs pipeline_inputs --output pipeline_results --threads 4
```

Acquire `pipeline_inputs/` using `python run_pipeline.py download-inputs --destination pipeline_inputs`. It is a self-contained, checksum-recorded snapshot: 303 files, approximately 1.69 GB. Assay matrices are not committed to Git. [Direct sources and every file](../data/README.md). The runner checks every input hash before fitting, copies the needed cached sources into the output folder, runs Python and R in dependency order, and writes separate logs. It uses the current Python interpreter and `Rscript` on PATH. If necessary, supply `--rscript /path/to/Rscript`; paths containing spaces must be quoted. Output folders are separate from input folders and existing analysis archives.

This is a pipeline for deposited count matrices, not raw FASTQ alignment. The large mouse-deletion and single-cell stages take substantial time and memory; use a machine with at least 16 GB RAM and several GB of free output space. Linux, macOS, or Windows through WSL can use the same code. Native Windows also needs a compatible R environment. No workstation path, local virtual environment, previous paper image, email or collaborator account is required by the analysis.

## Install dependencies

The recorded environment uses Python 3.12, R 4.5.3, limma 3.66.0, edgeR 4.8.2 and statmod 1.5.2. Python package pins are in `pipeline/requirements.txt`.

With conda or micromamba:

```sh
conda env create -f pipeline/environment.yml
conda activate osdr-cns
python -m pip install -r pipeline/requirements.txt
Rscript pipeline/install_R.R
```

Alternatively, install Python 3.12 and R 4.5.x using your normal system tools, install the same Python requirements, and run the R installer. Building the pinned Bioconductor source packages requires the platform's R build tools. The pipeline's `check` command reports missing dependencies and records the versions actually used; changed environments require a new output folder. Small floating-point differences between platforms can occur. The complete-family and source checks remain mandatory.

```sh
python run_pipeline.py check --inputs pipeline_inputs
python run_pipeline.py list
```

## Resume, verify, or redraw

```sh
python run_pipeline.py run --inputs pipeline_inputs --output pipeline_results --resume
python run_pipeline.py verify --inputs pipeline_inputs --output pipeline_results
python run_pipeline.py figures --inputs pipeline_inputs --output pipeline_results
python run_pipeline.py report --inputs pipeline_inputs --output pipeline_results
```

Resume skips completed stages with unchanged code, input hashes and package versions. Changed stage code invalidates that stage and every subsequent dependency. Changed input bundles or package versions require a new output folder. An interrupted process can leave `.pipeline.lock`; remove that file only after confirming that its recorded process has stopped. Do not manually edit fitted tables and then resume; use `verify` to check saved numerical outputs or rerun into a fresh folder.

`--through STAGE_NAME` deliberately ends at a named stage and records a partial run. Omit it to execute the whole study. Figures and verification commands require the corresponding analysis outputs. The default full run prohibits network requests. `--allow-network` explicitly permits primary repository requests, but live sources may differ from the recorded snapshot; a new source version should be archived as a new input bundle rather than silently replacing this study's definition.

## Included analyses and biological units

| Data | Analysis |
|---|---|
| OSD-561/562 | Endogenous STAR DESeq2, rounded RSEM sensitivity, edgeR/voom, full age/cohort contrasts, paired-region and exact tests |
| OSD-682/685/698/699 / GSE239336 | GeoMx log-expression contrasts, Welch sensitivity; unresolved animal links retained |
| OSD-352 | Donor-level RNA/ATAC pseudobulk, held-donor reference-mixture validation, promoter context, Visium local structure and mouse-level exposure tests |
| OSD-563/564 | Independent RR10 bulk and paired-region context |
| OSD-536 | RR9 nCounter normalization, fixed detection gate and exact animal tests |
| OSD-738 | Submitted-workbook agreement and protein-source scale audit; no invented value repair or pooling |
| OSD-927–932 / GSE295428 | Eight shared RR8 mice: six enriched populations, three primary sensitivity gates, article 5% gate, two count engines, animal deletions, marker-relative and depth endpoints |
| OSD-525/457 | RR23 habitat/vivarium contrasts and MHU-3 genotype/flight models; all 50 actual mouse-deletion refits |
| OSD-612/613 | Shared versus different animal identity checks and cerebral-region models; missing deposited samples retained |
| RR8 forebrain states | Exposure-blind non-panel PCA and K=4/6 clustering, mouse-state counts, technical-PC sensitivity, composition checks |
| All eligible branches | Complete 151-target families, 50 Hallmark sets where available, correlation/df sensitivity, exact animal tests, evidence distributions, publication figures |

Cells, spots, technical libraries, sections, shared regions and reused modalities do not add independent exposure animals. Unavailable effects remain missing with hypothesis p=1 in the planned correction family. Zero corrected discoveries remain reported. State learning/testing shares mice and is exploratory. These RNA analyses do not establish protein activation, EV cargo or a causal mechanism.

## Outputs

`pipeline_results/ANALYSIS_REPORT.md` summarizes the current run. `core/`, `extended/` and `deep/` contain `prepared/`, `tables/`, `sources/` and `figures/`. Every stage has a log in `logs/`. `metadata/` records the input manifest, environment, stage status/fingerprints, analysis family summaries, PDF checks and output hashes. Generated source provenance uses output-relative paths so results can move with their sources.

Figures use black integrated headings and embedded fonts. Arial is used when installed, with Liberation Sans or DejaVu Sans as portable alternatives. Workflow diagrams are drawn directly from vector elements and have no external image dependency. Numerical figure values are independently checked. PDF checks verify embedding and page bounds; final journal layout still needs visual inspection at the intended print size.

## Download controls

```sh
python run_pipeline.py download-inputs --destination pipeline_inputs
python run_pipeline.py download-inputs --destination pipeline_inputs --cache /path/to/verified/local/bundle --offline
python run_pipeline.py download-inputs --destination partial_inputs --only visium_data
```

The downloader uses only the published primary-source URLs and frozen study definitions. It verifies every accepted file, resumes only matching files, extracts checked ZIP byte ranges and writes a complete input manifest only when all 303 inputs exist and match. `--only` intentionally acquires a partial subset and does not make it eligible for a full run. `--cache` reuses local files only after matching the published hash. A changed source is an error, not an automatic update to the study.

The source acquisition command uses Python's standard library and can run before scientific package installation. The `check`/`run` commands require the pinned Python/R environment. The active analysis needs no historical workspace or external review service.
