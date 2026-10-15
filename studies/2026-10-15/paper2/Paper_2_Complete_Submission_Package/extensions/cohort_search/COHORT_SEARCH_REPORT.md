# Additional cohort search and acquisition

Dated 15 October 2026. The outcome-independent retrospective eligibility plan was written before new data were retrieved; its SHA256 is `a103dce4736d4d5399262945c63bb25b9777dfa57adbefcf7692423345559d90`. The existing RR23 and RRRM2 findings were already known, so this search is not a preregistered study.

## Actionable result

OSD-563 adds a **distinct-animal, same-launch external cerebellar cohort**. Its ten processed STAR junction files were acquired and verified against the NASA publisher MD5 sidecar. Every retained matrix column was independently reconstructed using a second, line-by-line implementation. The frozen matrix contains 590,325 unique junctions and ten animals, with five Flight and five Ground Control. No expression, junction effect, candidate significance or direction was inspected to choose this cohort. The parent analysis will use its existing frozen directions and full feature families.

Inputs are `inputs/RR10/RR10_animal_junction_counts.csv.gz` and `inputs/RR10/RR10_animal_metadata.csv`. The matrix key is chromosome:strand:intron-start:intron-end. Filters exactly follow the existing validation pipeline: canonical motifs 1–6, strand 1/2, STAR annotation status 1, maximum overhang at least 8, positive uniquely mapped read support, and chromosomes 1–19/X/Y/MT. All count values are integer unique-read counts. No outcome-based measurement filter was applied.

The cohort consists of female B6129SF2/J wild-type mice aged 14–15 weeks at launch. These are distinct source names and RFIDs from the male C57BL/6J RR23 discovery animals. Both studies launched on **SpaceX-21**. RR10 is therefore a new biological-animal cohort, but it does not establish replication across independent launches. All RR10 animals received Alizarin Red, Calcein and BrdU marker injections. They are not untreated animals, and their different strain, sex, handling and preservation limit direct comparison. The study includes wild type animals only; the p21-null arm mentioned in the broader mission description is not in this dataset.

The publisher protocol specifies STAR 2.7.11b, GRCm39 with ERCC92 and Ensembl release 112. The publisher-MD5-verified processing-information ZIP contains the STAR commands, confirming two-pass Basic and the GRCm39 reference. All ten biological source animals have one bulk junction file. The downloaded software-versions Markdown has a malformed 33-character publisher-sidecar MD5 digest; a second download returned identical bytes; that exception is preserved in `tables/RR10_publisher_checksum_exception.json`. Do not describe every metadata file as checksum verified. All ten STAR count files and the processing-information ZIP do match publisher checksums.

## Design limitations that the new data do not repair

| Feature | RR10 Flight | RR10 Ground |
|---|---|---|
| Animals | 5 | 5 |
| Euthanasia | 2 January 2021, on ISS | 5 January 2021, on Earth |
| BSP dissection | 12 May 2021 | 14 May 2021 |
| Preserved material | Frozen carcass, Cryochiller | Frozen carcass, Cryochiller |
| Diet | Nutrient Upgraded Rodent Food Bar | Nutrient Upgraded Rodent Food Bar |

Exposure remains perfectly aliased with both euthanasia date and dissection date. This cohort cannot isolate microgravity causality. Its frozen-carcass preparation also differs from RR23 live return and direct tissue preservation. The repeated ALSDA Subject ID field is not a valid individual identifier in this ISA; biological identities are taken from unique source names and RFIDs. Animal identity is verified from deposited metadata, not through independent genotype fingerprinting.

Only 252 unrestricted label assignments exist with 5+5 animals. A prespecified positive-direction global score can resolve a minimum one-sided p of 1/252. An absolute two-sided statistic has a minimum p at least 2/252. With all 923 local clusters, Bonferroni/Holm rejection at 0.05 cannot be reached even at the one-sided minimum. Absence of local corrected discoveries would not demonstrate biological absence.

## Search coverage and complete decisions

NASA OSDR live searches for “cerebellum”, “cerebellar” and “brain” returned 375, 232 and 4,777 records, respectively; all reported records were returned. The OSD accession subsets contained 12, 0 and 64 records, representing 72 unique OSD studies. A prior workspace mouse-brain catalog adds OSD-374, giving **73 accessions with explicit inclusion/exclusion decisions** in `tables/complete_accession_decisions.csv`. Detailed live metadata and file inventories were frozen for 22 plausibly relevant cohorts. The broader search endpoint also returns non-OSD records; these were not all individually screened. This is explicitly a bounded NASA-accession search, supported by literature discovery.

Two corrected date-bounded Europe PMC queries returned all 58 and 2 records, respectively. All 60 have title/abstract decisions in `tables/complete_literature_title_abstract_screen.csv`. Initial queries are also preserved: the unquoted hyphenated mission query was overbroad (11,228 hits, only 1,000 returned) and was replaced rather than misreported as a complete search. The corrected queries include a publication-date bound through 15 October 2026. Literature screening was not an exhaustive full-text systematic review.

This bounded search identified **no further independently launched bulk cerebellar RNA-seq cohort** beyond RRRM2/OSD-561. That statement is limited to the sources and scope searched; it does not prove that none exists elsewhere.

## Secondary evidence and why it does not replace independent cerebellar replication

| Dataset | What it can provide | Why it is not another primary replication |
|---|---|---|
| OSD-222, STS-135 | Independent-launch cerebellar microarray context | Gene-level array measurements do not recover the frozen STAR junctions or Tcaf1 first-exon usage. |
| OSD-352, RR3 | Independent-launch bulk/single-nucleus/spatial cerebral tissue context | Bulk samples are cerebral-hemisphere sections at hippocampal level, not a separately sampled cerebellum; only 3+3 bulk animals. |
| OSD-562/613, RRRM2 | Within-mission hippocampal/cerebral regional comparison with OSD-561 | Same mission and source-animal pool. Matched-animal tissue differences could be useful, but are not independent replication. |
| OSD-564/612, RR10 | Within-mission hippocampal/cerebral regional comparison with OSD-563 | Same RR10 source-animal pool and shared launch with RR23. Tissue contrasts need animal matching and must preserve the marker-injection context. |
| OSD-927/930, RRRM1 SpaceX-16 | Independently launched cerebellar myeloid/non-myeloid cell context | Seven/eight 3′ 10x single-cell samples; no processed junction matrix in NASA inventory. These two tissues share the eight-animal cohort, not two replications. |
| OSD-682/685/698/699, SpaceX-24 | Regional GeoMx gene-expression context | Cortex/hippocampus regions, with BuOE treatment arms; not cerebellar bulk junction sequencing. |
| OSD-536, SpaceX-12 | Independent-launch targeted brain gene-expression context | NanoString panel in the associated publication; no full junction counts. |
| OSD-916, RRRM1 | Brain miRNA context | miRNA sequencing cannot validate local mRNA splice junctions. |

RRRM1/OSD-927/930 uses female BALB/cAnNTac animals, with two animals in each age-by-exposure cell (3 and 8 months; four flight and four ground before the one missing myeloid sample). Flight animals were dissected on 16 January 2019 and ground animals on 17 January 2019. Age entries in individual NASA sample-factor records are incomplete, requiring careful source mapping before any analysis. The 3′ assay is particularly poorly suited to validating a **5′ untranslated/alternative-first-exon** Tcaf1 hypothesis. Deposited gene-count or cell-type expression would be contextual evidence, not proof of that junction mechanism. Full raw-read acquisition was deferred because no immediately usable processed junction matrix was identified, the reads exceed 170 GB combined, and technical reprocessing does not remove the assay/design mismatch. The complete raw inventory is retained for a separately planned study.

The most useful next step from this search is the complete frozen RR10 global transport test, while retaining RRRM2 as the independent-launch result and every shared-launch/design caveat. An additional positive result, if observed, would strengthen cross-animal transport; it would not make the microgravity effect causally identifiable or guarantee journal acceptance.

## Publisher sources

- NASA OSD-563 study: https://osdr.nasa.gov/bio/repo/data/studies/OSD-563
- Live metadata: https://osdr.nasa.gov/osdr/data/osd/meta/563
- Live files: https://osdr.nasa.gov/osdr/data/osd/files/563
- NASA brain search: https://osdr.nasa.gov/osdr/data/search?term=brain&size=10000
- RRRM1 cerebellar myeloid: https://osdr.nasa.gov/bio/repo/data/studies/OSD-927
- RRRM1 cerebellar non-myeloid: https://osdr.nasa.gov/bio/repo/data/studies/OSD-930
- Associated RRRM1 GEO accession: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE295428

All exact query URLs, retrieval times, response hashes, original inventories, successful downloads and checksum exceptions are retained in the accompanying tables and sources.
