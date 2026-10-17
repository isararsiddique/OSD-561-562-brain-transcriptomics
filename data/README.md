# Data availability and direct primary-source downloads

[Study overview](../README.md) · [Every file with a direct link](FILE_DOWNLOADS.md) · [Machine-readable source table](SOURCE_FILES.tsv) · [Download manifest](download_manifest.json) · [Input SHA256 manifest](input_manifest.json)

## What is stored here

The current repository contains **no assay count matrices, FASTQ files, Cell Ranger matrices, H5 expression objects or fitted model pickles**. Assay inputs are acquired into the Git-ignored `pipeline_inputs/` folder. Derived statistics, QC/audit tables and figures are in `results/`.

The supplied corrected target definition, audited mapping and frozen public metadata/orthology records are retained in `reference/`. These small supporting definitions make the same study reproducible even if live metadata/orthology services change. The MGI annotation snapshot is compressed. Provider metadata files are not new experimental data.

## Acquire and verify the full input bundle

```sh
python run_pipeline.py download-inputs --destination pipeline_inputs
python run_pipeline.py check --inputs pipeline_inputs
python run_pipeline.py run --inputs pipeline_inputs --output pipeline_results --threads 4
```

The full bundle contains **303 files, approximately 1.69 GB**, including some duplicated source copies required by the study's provenance checks. The downloader validates recorded byte lengths and SHA256 hashes; archived members also receive ZIP CRC32 checks. It downloads only the required ZIP members through verified byte ranges, including the author cell metadata inside the 3.93 GB RR3 archive. It does not silently replace a file when a source has changed. An interrupted acquisition can be resumed with the same command; mismatched existing files require a fresh destination. See [pipeline instructions](../pipeline/README.md).

## All NASA accessions

| Accession | Role | Primary study | Metadata | File inventory |
|---|---|---|---|---|
| OSD-561 | Original cerebellum bulk | [Study](https://osdr.nasa.gov/bio/repo/data/studies/OSD-561) | [JSON](https://osdr.nasa.gov/osdr/data/osd/meta/561) | [Files](https://osdr.nasa.gov/osdr/data/osd/files/561) |
| OSD-562 | Original hippocampus bulk | [Study](https://osdr.nasa.gov/bio/repo/data/studies/OSD-562) | [JSON](https://osdr.nasa.gov/osdr/data/osd/meta/562) | [Files](https://osdr.nasa.gov/osdr/data/osd/files/562) |
| OSD-682 | Original GeoMx CA1 | [Study](https://osdr.nasa.gov/bio/repo/data/studies/OSD-682) | [JSON](https://osdr.nasa.gov/osdr/data/osd/meta/682) | [Files](https://osdr.nasa.gov/osdr/data/osd/files/682) |
| OSD-685 | Original GeoMx DG | [Study](https://osdr.nasa.gov/bio/repo/data/studies/OSD-685) | [JSON](https://osdr.nasa.gov/osdr/data/osd/meta/685) | [Files](https://osdr.nasa.gov/osdr/data/osd/files/685) |
| OSD-698 | Original GeoMx frontal cortex | [Study](https://osdr.nasa.gov/bio/repo/data/studies/OSD-698) | [JSON](https://osdr.nasa.gov/osdr/data/osd/meta/698) | [Files](https://osdr.nasa.gov/osdr/data/osd/files/698) |
| OSD-699 | Original GeoMx cortex | [Study](https://osdr.nasa.gov/bio/repo/data/studies/OSD-699) | [JSON](https://osdr.nasa.gov/osdr/data/osd/meta/699) | [Files](https://osdr.nasa.gov/osdr/data/osd/files/699) |
| OSD-352 | RR3 bulk, RNA/ATAC and Visium | [Study](https://osdr.nasa.gov/bio/repo/data/studies/OSD-352) | [JSON](https://osdr.nasa.gov/osdr/data/osd/meta/352) | [Files](https://osdr.nasa.gov/osdr/data/osd/files/352) |
| OSD-563 | RR10 cerebellum bulk | [Study](https://osdr.nasa.gov/bio/repo/data/studies/OSD-563) | [JSON](https://osdr.nasa.gov/osdr/data/osd/meta/563) | [Files](https://osdr.nasa.gov/osdr/data/osd/files/563) |
| OSD-564 | RR10 hippocampus bulk | [Study](https://osdr.nasa.gov/bio/repo/data/studies/OSD-564) | [JSON](https://osdr.nasa.gov/osdr/data/osd/meta/564) | [Files](https://osdr.nasa.gov/osdr/data/osd/files/564) |
| OSD-536 | RR9 brain nCounter | [Study](https://osdr.nasa.gov/bio/repo/data/studies/OSD-536) | [JSON](https://osdr.nasa.gov/osdr/data/osd/meta/536) | [Files](https://osdr.nasa.gov/osdr/data/osd/files/536) |
| OSD-738 | RR1 submitted protein context | [Study](https://osdr.nasa.gov/bio/repo/data/studies/OSD-738) | [JSON](https://osdr.nasa.gov/osdr/data/osd/meta/738) | [Files](https://osdr.nasa.gov/osdr/data/osd/files/738) |
| OSD-927 | RR8 enriched population | [Study](https://osdr.nasa.gov/bio/repo/data/studies/OSD-927) | [JSON](https://osdr.nasa.gov/osdr/data/osd/meta/927) | [Files](https://osdr.nasa.gov/osdr/data/osd/files/927) |
| OSD-928 | RR8 enriched population | [Study](https://osdr.nasa.gov/bio/repo/data/studies/OSD-928) | [JSON](https://osdr.nasa.gov/osdr/data/osd/meta/928) | [Files](https://osdr.nasa.gov/osdr/data/osd/files/928) |
| OSD-929 | RR8 enriched population | [Study](https://osdr.nasa.gov/bio/repo/data/studies/OSD-929) | [JSON](https://osdr.nasa.gov/osdr/data/osd/meta/929) | [Files](https://osdr.nasa.gov/osdr/data/osd/files/929) |
| OSD-930 | RR8 enriched population | [Study](https://osdr.nasa.gov/bio/repo/data/studies/OSD-930) | [JSON](https://osdr.nasa.gov/osdr/data/osd/meta/930) | [Files](https://osdr.nasa.gov/osdr/data/osd/files/930) |
| OSD-931 | RR8 enriched population | [Study](https://osdr.nasa.gov/bio/repo/data/studies/OSD-931) | [JSON](https://osdr.nasa.gov/osdr/data/osd/meta/931) | [Files](https://osdr.nasa.gov/osdr/data/osd/files/931) |
| OSD-932 | RR8 enriched population | [Study](https://osdr.nasa.gov/bio/repo/data/studies/OSD-932) | [JSON](https://osdr.nasa.gov/osdr/data/osd/meta/932) | [Files](https://osdr.nasa.gov/osdr/data/osd/files/932) |
| OSD-525 | RR23 cerebellum / habitat and vivarium | [Study](https://osdr.nasa.gov/bio/repo/data/studies/OSD-525) | [JSON](https://osdr.nasa.gov/osdr/data/osd/meta/525) | [Files](https://osdr.nasa.gov/osdr/data/osd/files/525) |
| OSD-457 | MHU-3 cerebrum / Nrf2 genotype | [Study](https://osdr.nasa.gov/bio/repo/data/studies/OSD-457) | [JSON](https://osdr.nasa.gov/osdr/data/osd/meta/457) | [Files](https://osdr.nasa.gov/osdr/data/osd/files/457) |
| OSD-612 | RR10 added cerebral region, shared animals | [Study](https://osdr.nasa.gov/bio/repo/data/studies/OSD-612) | [JSON](https://osdr.nasa.gov/osdr/data/osd/meta/612) | [Files](https://osdr.nasa.gov/osdr/data/osd/files/612) |
| OSD-613 | RRRM-2 cerebral region, different animals | [Study](https://osdr.nasa.gov/bio/repo/data/studies/OSD-613) | [JSON](https://osdr.nasa.gov/osdr/data/osd/meta/613) | [Files](https://osdr.nasa.gov/osdr/data/osd/files/613) |
| OSD-611 | Metadata comparison only; no outcome model | [Study](https://osdr.nasa.gov/bio/repo/data/studies/OSD-611) | [JSON](https://osdr.nasa.gov/osdr/data/osd/meta/611) | [Files](https://osdr.nasa.gov/osdr/data/osd/files/611) |

OSD-611 is used for a metadata comparison only. OSD-738 supplies a protein-source audit, not a pooled RNA/protein outcome. OSD-612 reuses RR10 animals. OSD-927–932 share eight mice across enriched populations. Accession count is not an independent-mission or independent-animal count.

OSD-682/685/698/699 regional assignments follow the project/deposited registry; the corresponding four GEO matrices are CA, DG, FCT and CT. OSD numbers and GLDS count-file prefixes may differ: for example, OSD-561 uses GLDS-556 and OSD-562 uses GLDS-557 deposited bulk filenames.

## GEO and versioned annotation sources

| Resource | Landing/source page | Direct files |
|---|---|---|
| Original GeoMx | [GSE239336](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE239336) | [CA](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE239nnn/GSE239336/suppl/GSE239336_CA_GeneExpression_Q3norm.txt.gz), [DG](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE239nnn/GSE239336/suppl/GSE239336_DG_GeneExpression_Q3norm.txt.gz), [FCT](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE239nnn/GSE239336/suppl/GSE239336_FCT_GeneExpression_Q3norm.txt.gz), [CT](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE239nnn/GSE239336/suppl/GSE239336_CT_GeneExpression_Q3norm.txt.gz) |
| RR8 matrices | [GSE295428](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE295428) | All 141 barcode/feature/matrix URLs in [FILE_DOWNLOADS.md](FILE_DOWNLOADS.md) |
| Promoter annotation, mm10/GRCm38 | Ensembl release 102 | [GTF.gz](https://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz) |
| Mouse Hallmark sets | MSigDB 2025.1.Mm, all 50 sets | [GMT](https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Mm/mh.all.v2025.1.Mm.symbols.gmt) |
| Human/mouse orthology | MGI Homology report | [Current report](https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt); frozen version included for this run |
| Ensembl lookup/homology | [Ensembl REST](https://rest.ensembl.org/) | Frozen response records in `reference/sources/core/` |

All individual bulk counts, deposited annotation tables, RCC files, submitted protein tables and archive members have direct primary-source links in [FILE_DOWNLOADS.md](FILE_DOWNLOADS.md). ZIP-member rows identify both the source archive and required member; clicking the archive link downloads the archive, while the pipeline extracts only required members by byte range.

## Source/version policy

The analysis snapshot is dated **17 October 2026**. Exact hashes define the version analyzed. Current landing pages or metadata endpoints may change; the pipeline uses the recorded snapshots and refuses different assay contents. One generated download ledger was converted to portable relative paths for publication; assay bytes and biological inputs are unchanged.

GEO describes a 50% mitochondrial gate for RR8 while the article Methods describe 5%; both versions are separately analyzed. An invalid standalone Visium RDS is excluded; the verified archive supplies the analyzed spatial inputs. Raw FASTQ alignment, unsupported sample reconstruction and imputed missing animal/states are outside this pipeline.

Data providers and primary publications retain their terms and attribution. The panel is a collaboration-supplied study definition without a separate public accession. No private email text is included. See [references](../docs/REFERENCES.md).
