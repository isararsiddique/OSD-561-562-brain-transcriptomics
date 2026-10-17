# CNS transcriptomics: reproducible analysis report

This run analyzes the supplied 151-target CNS panel against deposited count matrices and verified biological identities. The complete planned testing families are retained, with unavailable effects left missing and their hypothesis p-values set to one. Every numerical summary below comes from this output folder.

44 distinct mapped RR8 forebrain genes pass q < 0.05 in all eight method/QC combinations. Those combinations reuse the same eight mice and do not add independent biological replication.

## Complete testing families

| Analysis | Specification | Planned tests | q < 0.05 | Minimum q |
| --- | --- | --- | --- | --- |
| Original bulk / DESeq2 |  | 4228 | 17 | 3.48013e-06 |
| Original GeoMx |  | 4228 | 0 | 0.20092 |
| Bulk method comparison | ('edgeR_QL_TMM',) | 4228 | 10 | 0.0261817 |
| Bulk method comparison | ('limma_voom_TMM',) | 4228 | 0 | 0.0799437 |
| Independent bulk | ('edgeR_QL_TMM',) | 453 | 1 | 0.0452009 |
| Independent bulk | ('limma_voom_TMM',) | 453 | 0 | 0.307083 |
| RR8 enriched populations | ('alternative_gate', 'edgeR_QL') | 3624 | 226 | 1.38133e-06 |
| RR8 enriched populations | ('alternative_gate', 'limma_voom') | 3624 | 224 | 5.63634e-08 |
| RR8 enriched populations | ('marker_positive_gate', 'edgeR_QL') | 3624 | 206 | 1.524e-07 |
| RR8 enriched populations | ('marker_positive_gate', 'limma_voom') | 3624 | 220 | 3.09548e-09 |
| RR8 enriched populations | ('published_gate', 'edgeR_QL') | 3624 | 268 | 2.03161e-07 |
| RR8 enriched populations | ('published_gate', 'limma_voom') | 3624 | 264 | 3.28131e-09 |
| Additional and regional bulk | ('edgeR_QL',) | 1812 | 10 | 0.000221043 |
| Additional and regional bulk | ('limma_voom',) | 1812 | 9 | 0.000118036 |
| Article 5% mitochondrial gate | ('edgeR_QL',) | 3624 | 235 | 1.29403e-07 |
| Article 5% mitochondrial gate | ('limma_voom',) | 3624 | 229 | 3.37575e-09 |

## Selected biological contrasts

| OSD | analyte | method | contrast | log2FC | q_expanded_1812_family |
| --- | --- | --- | --- | --- | --- |
| OSD-525 | AIM2 | limma_voom | flight_vs_habitat_ground | +0.741 | 0.014624 |
| OSD-525 | PYCARD | limma_voom | flight_vs_habitat_ground | -0.176 | 0.98734 |
| OSD-525 | S100B | limma_voom | flight_vs_habitat_ground | -0.162 | 0.97925 |
| OSD-525 | SLC1A3 | limma_voom | flight_vs_habitat_ground | +0.064 | 1 |
| OSD-525 | AIM2 | edgeR_QL | flight_vs_habitat_ground | +0.733 | 0.010846 |
| OSD-525 | PYCARD | edgeR_QL | flight_vs_habitat_ground | -0.179 | 0.98298 |
| OSD-525 | S100B | edgeR_QL | flight_vs_habitat_ground | -0.173 | 0.9315 |
| OSD-525 | SLC1A3 | edgeR_QL | flight_vs_habitat_ground | +0.055 | 1 |
| OSD-525 | AIM2 | limma_voom | flight_vs_vivarium | +0.208 | 0.84283 |
| OSD-525 | PYCARD | limma_voom | flight_vs_vivarium | -0.637 | 0.070448 |
| OSD-525 | S100B | limma_voom | flight_vs_vivarium | -0.104 | 1 |
| OSD-525 | SLC1A3 | limma_voom | flight_vs_vivarium | +0.088 | 1 |
| OSD-525 | AIM2 | edgeR_QL | flight_vs_vivarium | +0.174 | 0.89783 |
| OSD-525 | PYCARD | edgeR_QL | flight_vs_vivarium | -0.630 | 0.079357 |
| OSD-525 | S100B | edgeR_QL | flight_vs_vivarium | -0.124 | 1 |
| OSD-525 | SLC1A3 | edgeR_QL | flight_vs_vivarium | +0.090 | 1 |
| OSD-457 | AIM2 | limma_voom | flight_WT | +0.375 | 0.45326 |
| OSD-457 | PYCARD | limma_voom | flight_WT | +0.440 | 0.49469 |
| OSD-457 | S100B | limma_voom | flight_WT | -0.154 | 0.44415 |
| OSD-457 | SLC1A3 | limma_voom | flight_WT | -0.214 | 0.017233 |
| OSD-457 | AIM2 | edgeR_QL | flight_WT | +0.325 | 0.51848 |
| OSD-457 | PYCARD | edgeR_QL | flight_WT | +0.466 | 0.46001 |
| OSD-457 | S100B | edgeR_QL | flight_WT | -0.153 | 0.40947 |
| OSD-457 | SLC1A3 | edgeR_QL | flight_WT | -0.214 | 0.020735 |
| OSD-457 | AIM2 | limma_voom | flight_balanced_genotype | +0.316 | 0.35308 |
| OSD-457 | PYCARD | limma_voom | flight_balanced_genotype | +0.244 | 0.73269 |
| OSD-457 | S100B | limma_voom | flight_balanced_genotype | -0.189 | 0.044207 |
| OSD-457 | SLC1A3 | limma_voom | flight_balanced_genotype | -0.175 | 0.0079351 |
| OSD-457 | AIM2 | edgeR_QL | flight_balanced_genotype | +0.305 | 0.29792 |
| OSD-457 | PYCARD | edgeR_QL | flight_balanced_genotype | +0.272 | 0.59458 |
| OSD-457 | S100B | edgeR_QL | flight_balanced_genotype | -0.193 | 0.036288 |
| OSD-457 | SLC1A3 | edgeR_QL | flight_balanced_genotype | -0.175 | 0.0082406 |
| OSD-612 | AIM2 | limma_voom | flight_vs_ground | -0.029 | 1 |
| OSD-612 | PYCARD | limma_voom | flight_vs_ground | +0.613 | 0.6597 |
| OSD-612 | S100B | limma_voom | flight_vs_ground | -0.804 | 0.45791 |
| OSD-612 | SLC1A3 | limma_voom | flight_vs_ground | -0.409 | 0.67037 |
| OSD-612 | AIM2 | edgeR_QL | flight_vs_ground | -0.017 | 1 |
| OSD-612 | PYCARD | edgeR_QL | flight_vs_ground | +0.436 | 0.93309 |
| OSD-612 | S100B | edgeR_QL | flight_vs_ground | -0.981 | 0.37066 |
| OSD-612 | SLC1A3 | edgeR_QL | flight_vs_ground | -0.449 | 0.60309 |
| OSD-613 | AIM2 | limma_voom | flight_balanced_age | -0.339 | 0.51301 |
| OSD-613 | PYCARD | limma_voom | flight_balanced_age | -0.265 | 0.45326 |
| OSD-613 | S100B | limma_voom | flight_balanced_age | -0.030 | 1 |
| OSD-613 | SLC1A3 | limma_voom | flight_balanced_age | -0.055 | 1 |
| OSD-613 | AIM2 | edgeR_QL | flight_balanced_age | -0.331 | 0.523 |
| OSD-613 | PYCARD | edgeR_QL | flight_balanced_age | -0.264 | 0.43723 |
| OSD-613 | S100B | edgeR_QL | flight_balanced_age | -0.032 | 1 |
| OSD-613 | SLC1A3 | edgeR_QL | flight_balanced_age | -0.057 | 1 |

All target effects, including uncorrected and unavailable estimates, are in the full CSV tables. Mouse-deletion results are in [the expanded deletion summary](deep/tables/reported_expanded_deletion_summary.csv). Expression-state, technical-state, fixed-composition, depth, paired-region, pathway, spatial and promoter checks retain their original families and experimental units.

## Interpretation limits

- Original terminal bulk preservation differs by flight group. RR8 exposure and processing day are confounded. RR23 collection dates and habitat/vivarium comparator differences matter.
- Cells, nuclei, technical libraries, spatial spots, sections and shared modalities do not increase the number of exposure animals. Original GeoMx animal links remain unresolved.
- Expression states are numerical exploratory clusters learned and tested in the same mice. RNA fractions are not absolute cell fractions; adjusted associations are not causal mediation.
- CAMERA results depend on intergene correlation and residual degrees of freedom. An uncorrected test or zero discovery count is not evidence of equivalence.
- These data measure transcripts. They do not establish EV cargo, protein activation, inflammasome function, protection or a first-ever mechanism. Protein-source scale flags remain recorded.

## Verification and outputs

The independent core, extension, figure-value and deeper audits completed 931 checks (56 + 134 + 91 + 650). Source integrity is recorded in [the input manifest](../data/input_manifest.json), package versions in [the environment](verification/environment.json), and file hashes in [the output manifest](publication_manifest.json).

Figures are in [core](core/figures/), [integrated analyses](extended/figures/) and [deeper validation](deep/figures/). Publication figures use black headings and the installed Arial font, embedded in PDFs. Every stage has a separate log. The complete script order is available with `python run_pipeline.py list`.
