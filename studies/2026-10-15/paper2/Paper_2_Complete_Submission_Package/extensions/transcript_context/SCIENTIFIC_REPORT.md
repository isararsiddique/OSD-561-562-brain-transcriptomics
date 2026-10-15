# Exact transcript-context audit of RR23 local junction measurements

Author: Isarar Siddique

Affiliation:

Email:

Completion date: 15 October 2026. Annotation and descriptive audit of existing retrospective outcomes; no new hypothesis-testing family or independent expression assay.

## Full universe and annotation limits

The frozen whole-gene GRCm39.112 context contains 2,381 genes and 404,117 feature records selected for every gene relevant to the measured universe, without significance screening. The union includes 8,825 distinct needed junctions from all 923 multiway clusters and 3,239 uniquely gene-mapped binary choices. Of these, 7,626 have exact transcript-intron matches. Within the 3,248 distinct multiway junctions, 2,049 match and 1,199 do not. All unmatched coordinates remain in multiway_junction_exact_match_status.csv and the full cluster audit; no nearby gene is substituted.

The original STAR canonical/annotation-status gate and exact chromosome/strand GTF transcript-intron membership are different criteria. The mismatch is a real audit limit, not grounds to silently discard original hypotheses or assign full-transcript interpretations. Only 293/923 clusters have every alternative exactly matched; 290 have a unique common-gene intersection, and seven of these contain an alternative individually matching multiple genes while their all-alternative intersection remains one gene. Among the original 41 corrected clusters, 22 are fully matched/uniquely assigned; 19 have limited annotation. The original coordinate-based tests and gene-intersection mappings are preserved.

All 3,239 binary event annotations match their original common gene. All 1,983 measured genes, including three blank symbols, are retained. Gene IDs, not invented symbols, identify those three records:

| gene_id | gene_symbol | measured_events | all_event_annotations_match | primary_BY_gene |
| --- | --- | --- | --- | --- |
| ENSMUSG00000121470 | blank | 2 | True | False |
| ENSMUSG00000121502 | blank | 1 | True | False |
| ENSMUSG00000121596 | blank | 1 | True | False |

## Tcaf1 exact reference architecture

For the minus-strand cluster 6:2:42663893:start, shared low genomic intron start 42,663,893 lies on the acceptor side, while the alternative high intron ends are donor-side boundaries. The three inclusive intron intervals have lengths 5,974, 22,570 and 23,000 bp. The short intron matches Tcaf1-202/-203 (protein_coding) and Tcaf1-206 (retained_intron); the longest matches Tcaf1-201 (protein_coding) and Tcaf1-208 (protein_coding_CDS_not_defined). The middle intron matches Tcaf1-205 (retained_intron), with donor/acceptor exon ranks 2→3. All other matches have ranks 1→2. These are nonexclusive memberships: short reads at a shared intron cannot choose one of the matching complete transcripts.

| junction_id | intron_length_bp | transcript_name | transcript_biotype | donor_exon_rank | acceptor_exon_rank | donor_boundary_region | acceptor_boundary_region |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 6:2:42663893:42686892 | 23000 | Tcaf1-201 | protein_coding | 1 | 2 | 5prime_UTR | 5prime_UTR |
| 6:2:42663893:42669866 | 5974 | Tcaf1-203 | protein_coding | 1 | 2 | 5prime_UTR | 5prime_UTR |
| 6:2:42663893:42669866 | 5974 | Tcaf1-202 | protein_coding | 1 | 2 | 5prime_UTR | 5prime_UTR |
| 6:2:42663893:42669866 | 5974 | Tcaf1-206 | retained_intron | 1 | 2 | no_CDS_annotation | no_CDS_annotation |
| 6:2:42663893:42686892 | 23000 | Tcaf1-208 | protein_coding_CDS_not_defined | 1 | 2 | no_CDS_annotation | no_CDS_annotation |
| 6:2:42663893:42686462 | 22570 | Tcaf1-205 | retained_intron | 2 | 3 | no_CDS_annotation | no_CDS_annotation |

The coding transcripts Tcaf1-201/-202/-203 place both measured exonic junction boundaries in the 5′ untranslated region and share the same seven CDS coordinate/phase segments, 2,766 CDS bases and CCDS20076. Their first-exon intervals differ; their shared annotated coding structure does not demonstrate an altered protein. Tcaf1-205/-206/-208 have no CDS annotation; retained_intron and CDS-not-defined are reference biotypes, not observations of retained introns or noncoding function in these samples.

The CDS coordinate/phase structure, printed in genomic order, is [[42650312, 42650572, "0"], [42651976, 42652313, "2"], [42653725, 42653964, "2"], [42654264, 42654477, "0"], [42655004, 42655101, "2"], [42655360, 42656354, "1"], [42663259, 42663878, "0"]]. The figure joins these segments in reverse order to follow minus-strand transcript direction; introns are removed solely for drawing. Full exons, CDS/UTR features, transcript IDs, protein IDs and CCDS IDs are retained in the tables.

| transcript_name | transcript_biotype | n_exons | CDS_bases | five_prime_UTR_bases | CCDS_ids |
| --- | --- | --- | --- | --- | --- |
| Tcaf1-201 | protein_coding | 9 | 2766 | 144 | CCDS20076 |
| Tcaf1-203 | protein_coding | 9 | 2766 | 169 | CCDS20076 |
| Tcaf1-202 | protein_coding | 8 | 2766 | 141 | CCDS20076 |
| Tcaf1-206 | retained_intron | 6 | 0 | 0 | none |
| Tcaf1-208 | protein_coding_CDS_not_defined | 3 | 0 | 0 | none |
| Tcaf1-205 | retained_intron | 3 | 0 | 0 | none |
| Tcaf1-204 | protein_coding_CDS_not_defined | 3 | 0 | 0 | none |
| Tcaf1-207 | protein_coding_CDS_not_defined | 3 | 0 | 0 | none |

Exact chromosome 6 first-exon and 5′ UTR intervals for the three coding matches are 1-based inclusive and listed below. Each first exon is wholly annotated as 5′ UTR, followed in transcript direction by the shared 14-base second-exon UTR at 42,663,879–42,663,892. Their shared CDS begins immediately downstream in minus-strand transcript direction at genomic base 42,663,878. These reference coordinates do not measure transcription initiation or promoter activity.

| transcript_name | first_exon_start | first_exon_end | first_exon_bases | five_prime_UTR_intervals | five_prime_UTR_bases |
| --- | --- | --- | --- | --- | --- |
| Tcaf1-201 | 42686893 | 42687022 | 130 | 42686893–42687022; 42663879–42663892 | 144 |
| Tcaf1-202 | 42669867 | 42669993 | 127 | 42669867–42669993; 42663879–42663892 | 141 |
| Tcaf1-203 | 42669867 | 42670021 | 155 | 42669867–42670021; 42663879–42663892 | 169 |

## Empirical read support

All 26 biological animals contribute actual deposited count vectors, after inherited technical-library collapsing. The three introns are supported in every animal, with minimum total selected reads 80 in FLT, 79 in GC and 105 in VIV. Mean long-intron fractions are 0.781322, 0.617234 and 0.658661. The companion intron 6:2:42686572:42686892 is present with 20–169 unique reads in all 26 animals. These counts add same-library local support; they do not identify one full RNA molecule or independently validate its architecture.

| group | n | minimum_selected_reads | median_selected_reads | mean_long_fraction | mean_long_vs_short_conditional |
| --- | --- | --- | --- | --- | --- |
| FLT | 9 | 80.0 | 124.0 | 0.781322485 | 0.926503213 |
| GC | 9 | 79.0 | 126.0 | 0.617234055 | 0.863517301 |
| VIV | 8 | 105.0 | 134.5 | 0.658661 | 0.902689594 |

The conditional long/(long+short) fraction excludes the 22,570-bp alternative and is reported descriptively, not as a newly tested binary event. Coverage correlations in the complete table receive no new p values; pooled and within-arm correlations are not independent replication or proof of a sequencing mechanism:

| group | n | spearman_rho_fraction_vs_selected_reads | interpretation |
| --- | --- | --- | --- |
| all26 | 26 | -0.008551394 | Descriptive coverage relation; no gene/animal significance assigned and no independent replication |
| FLT | 9 | -0.366666667 | Descriptive coverage relation; no gene/animal significance assigned and no independent replication |
| GC | 9 | 0.610884008 | Descriptive coverage relation; no gene/animal significance assigned and no independent replication |
| VIV | 8 | 0.238095238 | Descriptive coverage relation; no gene/animal significance assigned and no independent replication |

## Relationship to statistical findings

Tcaf1 was selected for this structural description after the primary complete-family BY tests were known. Its initial q values are 0.032287 for FLT–GC and 0.045051 for FLT–VIV. The stronger 95% gate leaves it measurable but changes full-family q to 0.038744 and 0.053811, respectively; its FLT–VIV finding fails that sensitivity. Both comparisons fail full-family BY in all ten half-depth count perturbations despite remaining measurable. This audit does not change the original findings or rescue failed robustness/transport tests. Reference compatibility narrows a follow-up question toward alternative first-exon/5′ architecture, but does not establish promoter activity, transcript abundance, RNA stability, translation, complete isoforms, intron retention, protein function or injury. No broader literature search or priority claim was added.

## Provenance and independent verification

22 checks pass. An independent parser, without importing analysis functions, reconstructs every needed intron/transcript membership, donor/acceptor rank and CDS/UTR endpoint class; every transcript CDS/phase summary; all 923 common-gene intersections; all 3,239 binary gene annotations; every Tcaf1 animal/alternative count and every available annotated intron count. It verifies the frozen local GTF hash and publisher processing/software MD5/SHA values. The publisher reference names specify GRCm39 primary assembly plus ERCC92, Ensembl112, STAR 2.7.11b and two-pass Basic, consistent with this context. The complete FASTA sequence and original full GTF are not reconstructed by these completion checks; their source hash/URL and initial preparation provenance remain recorded.

Original provenance.json is unchanged. Historical byte hashes of regenerated parent CSV.gz are recorded alongside current values, while semantic checks establish unchanged gene IDs, definitions, memberships and raw animal counts. Gzip headers/timestamps may change without any numerical change. The verifier therefore does not mistake derived gzip byte drift for scientific data drift. Immutable publisher and frozen local GTF bytes remain required to match their recorded hashes.

All local outputs are portable, with --base pointing at the NEW package. README documents offline reproduction from the supplied frozen subset; initial preparation is separate and requires the original full GTF/publisher files. All input provenance, completion plan, code, full tables, verification, report and vector/raster figure are included in MANIFEST.csv.

## Design limits

These measurements reuse the same bulk cerebellar libraries as the binary, gene and multiway analyses. They do not add independent animals, cell-type measurements, long reads, proteins or function. Flight samples were collected 14 January 2021; both ground groups 17 January 2021. The flight/date alias persists through every annotation and robustness check. The annotated reference contexts support coordinate-level hypotheses for subsequent assays while limiting causal flight and protein-consequence claims. No outreach was performed.
