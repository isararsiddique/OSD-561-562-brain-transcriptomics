# Dataset links, experimental units and provenance

All source datasets are public NASA Open Science Data Repository (OSDR) deposits. Obtain the processed count/junction files, ISA metadata and associated checksum sidecars from each study's file inventory. The archived frozen inputs provide the exact versions used in these analyses; newly downloaded files should be checked against their recorded checksums before use.

| Accession | Paper 1 use | Paper 2 use | NASA study page |
|---|---|---|---|
| OSD-525 | RR-23 cerebellum; flight, habitat-ground and vivarium contrasts | RR23 discovery; 26 animals, 9 FLT / 9 GC / 8 VIV | [OSD-525](https://osdr.nasa.gov/bio/repo/data/studies/OSD-525) |
| OSD-561 | RRRM-2 cerebellum; protocol and matched-flight comparisons | RRRM2 external validation; 19 included animals, 10 FLT / 9 GC | [OSD-561](https://osdr.nasa.gov/bio/repo/data/studies/OSD-561) |
| OSD-562 | RRRM-2 hippocampus | Not included in Paper 2's cerebellar junction cohorts | [OSD-562](https://osdr.nasa.gov/bio/repo/data/studies/OSD-562) |
| OSD-563 | RR-10 cerebellum | RR10 external validation; 10 animals, 5 FLT / 5 GC | [OSD-563](https://osdr.nasa.gov/bio/repo/data/studies/OSD-563) |
| OSD-564 | RR-10 hippocampus | Not included | [OSD-564](https://osdr.nasa.gov/bio/repo/data/studies/OSD-564) |
| OSD-612 | RR-10 left cerebral hemisphere | Not included | [OSD-612](https://osdr.nasa.gov/bio/repo/data/studies/OSD-612) |
| OSD-613 | RRRM-2 left cerebral hemisphere | Not included | [OSD-613](https://osdr.nasa.gov/bio/repo/data/studies/OSD-613) |
| OSD-202 | External brain expression/prediction evaluation | Not included | [OSD-202](https://osdr.nasa.gov/bio/repo/data/studies/OSD-202) |
| OSD-352 | RR-3 external brain expression/prediction evaluation | Not included | [OSD-352](https://osdr.nasa.gov/bio/repo/data/studies/OSD-352) |
| OSD-374 | External evaluation; two-donor pools remain pool-level units | Not included | [OSD-374](https://osdr.nasa.gov/bio/repo/data/studies/OSD-374) |
| OSD-611 | External evaluation; primary NASA SampleTable pool design and full-pool sensitivity | Not included | [OSD-611](https://osdr.nasa.gov/bio/repo/data/studies/OSD-611) |

## Identity and cohort boundaries

Paper 1 uses seven primary deposits from three missions, comprising 129 libraries from 83 distinct mice. Multiple brain regions or technical libraries do not add independent animals. Its external evaluation includes individual animals and pools; pool donors are not treated as separately measured expression units.

Paper 2 uses 55 distinct cerebellar animals across three deposited cohorts. RR23 and RR10 share SpaceX-21, so the three cohorts span two launches. Distinct source names/RFIDs establish separate animals; they do not provide genotype fingerprinting or make the shared-launch cohorts independent launch replications. RRRM2 metadata conflicts were handled by the archived full-stratum exclusion, not by selecting on observed junction effects.

## Frozen-input provenance

The complete Paper 1 archive includes `inputs/provenance.json`, `inputs/upstream_count_provenance.json` and the external acquisition/design records. Source identity and technical-library collapse are checked before modelling.

The complete Paper 2 archive includes frozen integer junction counts, reference annotation, technical-library/source-animal maps, OSDR file inventories, official MD5 checks and local SHA-256 records. The original 29 RR23 junction-file MD5 values matched their processed-file sidecar. The GRCm39 release 112 annotation defines exact intron-coordinate matching.

Source manifests and the corresponding NASA inventories document the exact download URLs and filenames. Their measurement content remains in the full archives, not in this GitHub release. [SOURCE_FILE_MANIFEST.tsv](SOURCE_FILE_MANIFEST.tsv) lists all archived paths and hashes without copying the datasets.

## Additional reference material

Paper 1's published-signature analysis stores its dated, hashed plan and signature provenance in `research/F_signatures/`. The plan/hash gate must remain intact. Its exploratory Tyssowski gene-class analysis is explicitly labelled post hoc.

Paper 2's cell-type enrichment uses the supplementary marker tables from Kozareva et al., *Nature* 598, 214–219 (2021), [doi:10.1038/s41586-021-03220-z](https://doi.org/10.1038/s41586-021-03220-z). The external supplementary tables are not redistributed; derived marker sets and the extraction scripts are supplied. Online enrichment/reference steps query g:Profiler or Crossref, while frozen outputs are retained in the full archive.

No new animal work, sequencing or raw-read alignment was performed for either study. Changes in an upstream file, annotation release, dataset version or analysis environment can change regenerated results; use the supplied frozen inputs for exact reproduction.
