# Methods

[Study overview](../README.md) · [Source files](../data/FILE_DOWNLOADS.md) · [Findings](FINDINGS.md) · [Figure captions](FIGURES.md)

## Study definition and source integrity

The supplied corrected panel contains 151 human analytes and 143 unique primary mouse symbols after audited orthology mapping. Protein-style labels are mapped to their transcript-producing genes: ICE→Casp1, IL-18→Il18 and Amyloid-beta→App are examples. Analyte labels are preserved for traceability; those labels do not turn RNA counts into protein or EV measurements. Ambiguous mappings and unavailable assays are retained in the complete panel. The input and source manifests record exact filenames, lengths and SHA256 hashes. Public metadata and orthology snapshots are frozen in the repository; assay counts are acquired from primary repositories. OSD numbers and GLDS filename prefixes can differ.

This workflow begins with deposited count matrices and submitted assay tables. It does not realign FASTQ files, reproduce every original authors' preprocessing step, or treat the deposited differential-expression statistics as new independently fitted results. Deposited genome tables support annotation; the new count models use the deposited STAR counts. The rounded RSEM analysis is explicitly a source sensitivity.

## Biological identities and eligibility

ISA Source Name identifies animals where available. Technical libraries are collapsed or treated according to assay structure, rather than counted as independent animals. Original OSD-561/562 have 53 libraries, 29 animals and 24 paired-region animals. RR8 has 47 libraries across six enriched populations from eight shared mice; cerebellar myeloid has seven represented mice. OSD-612 measures nine already represented RR10 animals, while OSD-613 has 18 different RRRM-2 identities without overlap with the original 29. Missing count columns, absent donor/states and unresolved GeoMx links remain recorded.

RR23 has 29 technical libraries from 26 mice; technical libraries are summed by verified mouse. Flight, habitat-ground and vivarium controls contain 9, 9 and 8 mice and remain separate. MHU-3 has 24 deposited mice, six per genotype/exposure group. A published caption reports five flight KO mice; the deposited six are retained and each KO omission is evaluated in the complete deletion analysis.

## Count models and complete families

DESeq2 uses whole-gene endogenous STAR counts for the original bulk arm. Independent edgeR quasi-likelihood and limma/voom models use TMM normalization and whole-gene filtering before extracting the panel. Designs preserve eligible age, collection-cohort, exposure and interaction contrasts. Region-specific and mission-specific effects remain separate; no pooled pan-brain effect is estimated across incompatible designs.

Benjamini–Hochberg adjustment is recomputed over each complete planned family, including unavailable target hypotheses with p=1. An unavailable effect remains missing, rather than being imputed as zero. Repeated human labels mapping to the same mouse transcript remain traceable panel rows; counts of distinct genes are deduplicated when stated. Families are method/assay specific and do not correct across every exploratory analysis choice.

| Branch | Planned family |
|---|---:|
| Original bulk | 151 targets × 28 contrasts = 4,228 per engine |
| Original GeoMx | 151 × 28 = 4,228 |
| Earlier independent bulk | 151 × three contrasts = 453 per engine |
| RR8 enriched populations | 151 × six populations × four age/exposure contrasts = 3,624 per engine/QC |
| Additional/regional bulk | 151 × 12 contrasts = 1,812 per engine |
| RR3 promoter context | 151 × three eligible broad types = 453 |
| Visium local structure | 151 × 12 sections = 1,812 |
| Visium mouse exposure | 151 |
| State models | 151 × K × four contrasts per partition/engine |

Nominal 95% model intervals and corrected q-values describe different uncertainty rules. Mouse-deletion ranges are influence summaries, not confidence intervals. Filled plot markers use the relevant full-family q<0.05 threshold. Full-family outputs, including uncorrected estimates and failures, accompany focused figures.

## RNA-mixture and cellular context

The RR3 reference uses author-assigned broad cellular annotations. Non-panel reference transcripts support NNLS and ridge RNA-mixture estimation. Held-donor synthetic mixtures test numerical estimation without using the held donor to define the reference; 500 mixtures are evaluated per method. These are synthetic reference-validation samples, not additional experimental animals. Estimated fractions describe RNA mixtures, not absolute cell abundance. Adjustment for a highly collinear inferred fraction cannot establish causal mediation.

RR3 snRNA/ATAC uses 21,178 author-annotated nuclei from five donors, with modalities measured on the same nuclei. Donor-level pseudobulk and promoter context preserve that dependence. Promoter windows use mm10/GRCm38 Ensembl release 102 TSS ±2 kb. Three broad types meet eligibility. Promoter and RNA effects are contextual evidence, not independent replication of exposure.

## RR8 cell filtering and donor pseudobulk

The 47 attached matrices contain 156,746 cells before the studied gates. The primary gate follows the current GEO processing description: 200–2,500 detected features and ≤50% mitochondrial RNA, retaining 122,084 cells. Alternatives include a ≤20%/200–6,000-feature gate and a non-panel marker-positive gate. A separate article-defined ≤5%/200–2,500-feature sensitivity retains 62,351 cells across all 47 libraries. The article/GEO discrepancy is explicit. This analysis does not claim exact reproduction of author ambient-RNA or doublet correction.

Counts are aggregated to actual mouse/population/QC units. Models estimate age-specific, balanced-age exposure and exposure-by-age contrasts with two count engines. The eight engine/QC combinations are two engines × four filtering definitions in the same mice. Marker-relative models use five non-panel myeloid reference transcripts and HC3 uncertainty. They are biological comparators, not assumed invariant housekeeping genes. Analytical 500-UMI detection probabilities assess depth sensitivity; mouse means remain the exposure units. Age-blocked exact tests enumerate 36 labelings, or 18 where seven donors restrict the design.

## Exposure-blind expression states

The primary forebrain state input contains 23,036 cells from eight mice. Balanced training uses up to 1,500 cells per mouse, at most 12,000 total. Two thousand non-panel, non-mitochondrial, non-ribosomal variable features, 30 principal components and K-means K=4/K=6 define numerical partitions. Exposure and age labels are excluded from state learning. A technical sensitivity residualizes UMI and mitochondrial effects before clustering. Seeds and parameters are recorded in the preparation tables and scripts.

Within-state models require at least 50 cells in every mouse of a state. Missing donor/state coverage is not filled with zero expression. All eligible state hypotheses retain their specified correction family. Fixed-composition summaries compare measured state mixtures with equal state weights, but some states are sparse and primary K=6 has missing coverage. State learning and testing share the same eight mice, so partitions and fitted effects are exploratory rather than an independently validated taxonomy. State labels are not matched across alternative partitions.

## Spatial localization and exposure

The valid RR3 Visium archive supplies 36 checked matrix/coordinate/scaling members, 12 sections, six mice and 29,839 retained spots. The invalid standalone RDS is excluded. Local Moran analyses use 499 spatial permutations per section and the complete target-by-section family. Mouse-level expression averages technical sections before exact exposure tests over 20 group labelings. Local spatial structure and between-mouse exposure are different hypotheses.

Displayed sections use a predetermined first section per donor, raw full-resolution pixel coordinates and a common Aif1 expression scale. No unverified histological compartment or physical distance is assigned. Original GeoMx analyses use the deposited Q3-normalized matrices, log-expression factorial contrasts and Welch sensitivity. Unresolved animal/ROI links and slide/group confounding limit inference and prevent a causal antioxidant-rescue interpretation.

## Pathway and orthogonal assay checks

All 50 mouse Hallmark sets (MSigDB 2025.1.Mm) are retained in the eligible branch families: 1,400 core, 150 independent-bulk, 250 RR3 cell-type and 1,200 RR8 tests. CAMERA is evaluated at fixed intergene correlations 0.01 and 0.10, and with estimated correlation. Estimated CAMERA also changes residual sample degrees of freedom (3–19), so differences cannot be attributed to correlation alone. Animal-level singscores receive exact label tests and complete-family adjustment.

RR9 nCounter uses submitted RCC files, the stated normalization/detection gate and exact animal tests. Of 757 endogenous targets, 31 panel analytes meet detection coverage. AIF1 is absent and PYCARD is below the detection gate. OSD-738 is a submitted-workbook/source-scale audit: four extreme ASC values are present in the source and retained with flags. No invented rescaling or pooled RNA/protein inference is applied.

## Independent verification and figure construction

The pipeline records environment versions, stage fingerprints, input/output hashes and logs. Numerical checks independently recompute source agreement, family adjustment, donor counts, coverage, selected effects, deletion summaries and displayed values. Results matched 18 earlier complete tables (258,677 rows; rtol 1e-8, atol 1e-10). The complete 62-stage execution and relocation checks were run on macOS arm64.

Scientific figures are drawn from executed tables, with black integrated headings and embedded publication fonts. Arial is preferred; Liberation Sans and DejaVu Sans are portable alternatives. Workflow boxes/arrows are native vector elements. PNGs are previews, PDFs are publication exports and SVGs are editable where provided. The 19 PDFs passed single-page/font/text-bound checks; journal-specific column sizing still requires final layout inspection.
