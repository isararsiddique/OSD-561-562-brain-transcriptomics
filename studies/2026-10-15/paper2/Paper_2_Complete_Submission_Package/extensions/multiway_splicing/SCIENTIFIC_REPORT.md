# RR23 multiway local splice-choice extension

Author: Isarar Siddique

Affiliation:

Email:

Analysis date: 15 October 2026. Additional retrospective family; not independent confirmation of the previous binary-junction analysis.

## Scientific outcome

All three global multiway junction-composition comparisons survive Holm correction over the three primary tests. Across the complete local exploratory family of 2,769 tests, BY q≤0.05 identifies 40 flight–ground, three flight–vivarium, and zero ground–vivarium associations (43 contrast–cluster tests, 41 distinct clusters). BH over the same family identifies 155, 28, and 24, respectively. The global ground–vivarium finding therefore coexists with a null conservative local-discovery result; it does not license assigning uncorrected local events. These observations concern local compositions in the same bulk cerebellar animals, and do not establish complete isoforms, proteins, injury, or a causal flight effect.

## Frozen retrospective measurement

The dated plan was saved before loading multiway outcomes, with SHA256 8e22f6071305b131dbba786be9d17deeccf0d0eda0fae991404805aa8851311c. It is an analysis lock, not prospective preregistration. Starting from 805,562 canonical annotated junctions across 26 collapsed animals (flight 9, ground hardware 9, vivarium 8), 173,745 junctions qualify by ≥5 reads in ≥3 animals. Separate chromosome/strand/shared-genomic-start and shared-genomic-end groups require ≥3 qualifying alternatives, ≥20 selected reads in every animal, and ≥90% of all canonical boundary counts in every animal. Retain every qualifying alternative. Label-free pooled mean composition requires every retained alternative ≥1%, at least three ≥5%, and maximum ≤90%. Selection is fixed and is not repeated inside permutations or omissions. The complete 8,695 candidate-cluster table records every gate, including null and excluded candidates.

| boundary_type | qualifying_boundary_groups | at_least_three | coverage_pass | coverage_and_dominance_pass | selected |
| --- | --- | --- | --- | --- | --- |
| start | 148796 | 4322 | 3344 | 3028 | 466 |
| end | 148595 | 4373 | 3336 | 3019 | 457 |

The final 923 clusters comprise 466 shared starts and 457 shared ends. Alternative multiplicities are 3: 542 clusters, 4: 248 clusters, 5: 82 clusters, 6: 26 clusters, 7: 16 clusters, 8: 5 clusters, 9: 2 clusters, 10: 1 clusters, 12: 1 clusters. There are 3,376 alternative occurrences but 3,248 distinct retained junctions; 128 junctions each occur in two clusters, creating 128 intersecting cluster pairs. Complete overlap tables document this dependence. Equal cluster weight therefore does not mean equal gene or nonoverlapping molecular-event weight.

## Exact global inference

For each cluster e, H_e contains square roots of animal-level compositions and K_e=H_e H_eᵀ. K is the average of all K_e, giving every cluster equal weight. Its minimum eigenvalue is 0.00984769, and its diagonal differs from one by at most 3.77e-15. The primary statistic T=dᵀKd is the squared difference of group centroids in this sqrt space, and sqrt(T/2) is the Hellinger centroid distance. This is not the Hellinger distance of normalized arithmetic-mean compositions. Every fixed-size label assignment is enumerated; exact p counts ≥observed, including the observed assignment, with absolute tie tolerance 10⁻¹². The balanced comparison includes its label-complement tie. Holm correction applies to the three global tests. Exact arithmetic does not remove the biological exchangeability assumption.

| contrast | hellinger_centroid_distance | exact_p | n_ge_observed | n_assignments | holm_p_3 |
| --- | --- | --- | --- | --- | --- |
| FLT_vs_GC | 0.083129639 | 4.1135e-05 | 2 | 48620 | 0.000123406 |
| FLT_vs_VIV | 0.066147417 | 0.000658165 | 16 | 24310 | 0.001316331 |
| GC_vs_VIV | 0.065514911 | 0.002180173 | 53 | 24310 | 0.002180173 |

## Full local family

Each selected cluster receives the same complete exact centroid test for each of the three contrasts. BH and BY were calculated over one family of all 923×3=2,769 p values; BY controls under general dependence. Per-contrast corrections are also provided, but the declared conservative local finding uses BY over 3E. All tests, counts, and raw/adjusted p values are in local_exact_results_complete.csv.gz; the positive subset is BY_discoveries_full_family.csv. No pathway enrichment or full-isoform consequence is inferred.

| contrast | n_tests | minimum_raw_p | minimum_BH_q | minimum_BY_q | BH_positive | BY_positive |
| --- | --- | --- | --- | --- | --- | --- |
| FLT_vs_GC | 923 | 4.1135e-05 | 0.003796791 | 0.032286539 | 155 | 40 |
| FLT_vs_VIV | 923 | 4.1135e-05 | 0.003796791 | 0.032286539 | 28 | 3 |
| GC_vs_VIV | 923 | 0.000575895 | 0.014629839 | 0.124406848 | 24 | 0 |

## All three flight–vivarium local findings

Two of the three also survive full-family BY in flight–ground. Only the chromosome-6 cluster has a single gene consistent with every retained alternative. Partial matching of one alternative in the other clusters is insufficient for a gene assignment. All coordinates are 1-based inclusive intron intervals, chromosome and strand specified (STAR strand code 1=positive, 2=negative).

| cluster_id | gene_id | gene_symbol | gene_mapping_count | n_alternatives | exact_p | BY_q_full_3E | hellinger_centroid_distance | native_mean_total_variation | also_BY_significant_FLT_vs_GC |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 11:2:109886437:end | unassigned | unassigned | 0 | 3 | 4.1135e-05 | 0.032286539 | 0.107284776 | 0.145852178 | True |
| 2:1:55311420:start | unassigned | unassigned | 0 | 4 | 4.1135e-05 | 0.032286539 | 0.203467167 | 0.201932084 | False |
| 6:2:42663893:start | ENSMUSG00000036667 | Tcaf1 | 1 | 3 | 8.2271e-05 | 0.045050985 | 0.101115875 | 0.122661485 | True |

The complete intron alternatives and their individual exact gene matches are:

| cluster_id | alternative_index | junction_id | matching_gene_ids |
| --- | --- | --- | --- |
| 11:2:109886437:end | 0 | 11:2:109868700:109886437 | none |
| 11:2:109886437:end | 1 | 11:2:109872606:109886437 | ENSMUSG00000020620 |
| 11:2:109886437:end | 2 | 11:2:109873095:109886437 | none |
| 2:1:55311420:start | 0 | 2:1:55311420:55312239 | ENSMUSG00000087040 |
| 2:1:55311420:start | 1 | 2:1:55311420:55316969 | none |
| 2:1:55311420:start | 2 | 2:1:55311420:55317007 | none |
| 2:1:55311420:start | 3 | 2:1:55311420:55322064 | none |
| 6:2:42663893:start | 0 | 6:2:42663893:42669866 | ENSMUSG00000036667 |
| 6:2:42663893:start | 1 | 6:2:42663893:42686462 | ENSMUSG00000036667 |
| 6:2:42663893:start | 2 | 6:2:42663893:42686892 | ENSMUSG00000036667 |

Three-arm compositions below are means of animal proportions in alternative-index order. Entropy and coverage are descriptive; these were not new hypothesis tests. All animals and all three clusters, not just the displayed Tcaf1 example, are included in the characterization files.

| cluster_id | group | n_animals | mean_proportions_semicolon | mean_entropy_nats | mean_normalized_entropy | mean_selected_reads | minimum_selected_reads | minimum_dominance |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 11:2:109886437:end | FLT | 9 | 0.06122262;0.37596398;0.56281340 | 0.845056 | 0.769203 | 171.777778 | 103 | 0.980952 |
| 11:2:109886437:end | GC | 9 | 0.10580462;0.52691941;0.36727597 | 0.91848 | 0.836036 | 134.555556 | 54 | 0.964286 |
| 11:2:109886437:end | VIV | 8 | 0.10105816;0.48198062;0.41696122 | 0.927445 | 0.844196 | 144.125 | 102 | 0.95 |
| 2:1:55311420:start | FLT | 9 | 0.68194679;0.09632418;0.17446707;0.04726196 | 0.91744 | 0.661793 | 154.333333 | 107 | 0.953333 |
| 2:1:55311420:start | GC | 9 | 0.84980908;0.04360130;0.09083979;0.01574983 | 0.525404 | 0.378999 | 154.0 | 73 | 0.993103 |
| 2:1:55311420:start | VIV | 8 | 0.88387887;0.04627577;0.06009560;0.00974976 | 0.424926 | 0.306519 | 152.0 | 111 | 0.992857 |
| 6:2:42663893:start | FLT | 9 | 0.06177270;0.15690482;0.78132249 | 0.648223 | 0.590038 | 124.777778 | 80 | 0.976378 |
| 6:2:42663893:start | GC | 9 | 0.09808807;0.28467787;0.61723406 | 0.871156 | 0.79296 | 166.666667 | 79 | 0.95082 |
| 6:2:42663893:start | VIV | 8 | 0.06937388;0.27196512;0.65866100 | 0.79386 | 0.722602 | 147.875 | 105 | 0.967949 |

At the uniquely mapped Tcaf1 boundary, the intron ending at 42,686,892 has mean proportion 0.781 in flight, 0.617 in ground hardware, and 0.659 in vivarium; the intron ending at 42,686,462 has means 0.157, 0.285, and 0.272. The shorter intron ending at 42,669,866 has means 0.062, 0.098, and 0.069. Thus its flight–vivarium difference mainly redistributes composition between the two longer introns. Mean normalized entropy is 0.590, 0.793, and 0.723, respectively. These local count relationships do not identify a complete transcript or a protein consequence.

For the unassigned chromosome-11 negative-strand cluster, the intron starting at 109,873,095 increases from vivarium mean 0.417 to flight mean 0.563, while the other two alternatives decrease; ground hardware lies at 0.367 for that intron. For the unassigned chromosome-2 positive-strand cluster, the shortest intron ending at 55,312,239 decreases from vivarium mean 0.884 to flight mean 0.682, with higher shares of all three longer alternatives. Ground hardware has mean 0.850 for this shortest intron, but the flight–ground local test fails the declared full-family BY threshold (q=0.0891). Coordinate-level redistribution is reported for these two clusters because most alternatives lack an exact transcript-intron gene match; proximity-based gene naming would overstate the annotation evidence.

Flight–ground results for those same three clusters, including the negative result, are:

| cluster_id | exact_p | BY_q_full_3E | native_mean_total_variation |
| --- | --- | --- | --- |
| 11:2:109886437:end | 4.1135e-05 | 0.032286539 | 0.195537425 |
| 2:1:55311420:start | 0.000329083 | 0.089066315 | 0.167862297 |
| 6:2:42663893:start | 4.1135e-05 | 0.032286539 | 0.16408843 |

## Annotation and sensitivity

Exact GRCm39 release-112 introns were reconstructed from adjacent transcript exons (upstream exon end+1 through downstream exon start−1) with chromosome and strand matching. The intersection of exact gene sets over every alternative maps 290 clusters to 256 genes; 633 clusters have no common-gene intersection and none are assigned a gene by proximity. The local frozen lookup contains all 3,248 retained junctions, including explicit absent matches, and is sufficient for offline reproduction with these counts. Its raw GTF provenance/hash is retained in annotation/source.json. This mapping coverage limits equal-gene inference to a subset and limits biological naming of local discoveries.

On the 290 mapped clusters, equal-cluster and equal-gene sensitivities are reported separately. Equal-gene weight assigns 1/256 total weight per gene divided equally across its cluster definitions. The 95% dominance sensitivity retains 607 of the frozen clusters. All three global comparisons survive Holm within each explicitly secondary three-test family. These are measurements of robustness using the same animals, not independent replications.

| contrast | sensitivity | n_clusters | hellinger_centroid_distance | exact_p | holm_p_within_sensitivity_3 |
| --- | --- | --- | --- | --- | --- |
| FLT_vs_GC | mapped_subset_equal_cluster | 290 | 0.080260709 | 4.1135e-05 | 0.000123406 |
| FLT_vs_GC | equal_gene_weight | 290 | 0.07866238 | 4.1135e-05 | 0.000123406 |
| FLT_vs_GC | dominance_95 | 607 | 0.07587115 | 4.1135e-05 | 0.000123406 |
| FLT_vs_VIV | mapped_subset_equal_cluster | 290 | 0.05535739 | 0.002180173 | 0.004360346 |
| FLT_vs_VIV | equal_gene_weight | 290 | 0.055102358 | 0.002426985 | 0.00485397 |
| FLT_vs_VIV | dominance_95 | 607 | 0.058023861 | 0.000740436 | 0.001480872 |
| GC_vs_VIV | mapped_subset_equal_cluster | 290 | 0.058753826 | 0.003414233 | 0.004360346 |
| GC_vs_VIV | equal_gene_weight | 290 | 0.058225656 | 0.003743316 | 0.00485397 |
| GC_vs_VIV | dominance_95 | 607 | 0.05968316 | 0.003002879 | 0.003002879 |

Single-animal omission uses the frozen 923 definitions and reports effect ranges only; no repeated omission p-value family is used.

| contrast | n_omissions | minimum_distance | maximum_distance |
| --- | --- | --- | --- |
| FLT_vs_GC | 18 | 0.081491034 | 0.086240289 |
| FLT_vs_VIV | 17 | 0.065919698 | 0.069750509 |
| GC_vs_VIV | 17 | 0.064645955 | 0.07115801 |

## Independent verification

47 checks pass. A separate script reloads original counts, independently reconstructs all candidate and selected clusters, repeats pooled selection with animal order reversed, reconstructs every composition, intersects every alternative gene set, checks kernel PSD and unit diagonal, verifies the full global exact distributions from closed-form kernel sums, and verifies every one of the 2,769 local exact p values separately. SciPy independently checks all BH/BY calculations; Holm and all secondary exact results are also independently checked. Selection never loads group metadata. Detailed checked values are saved in independent_verification.csv and .json. The annotation lookup is checked against its gene-set intersections and source hash; the report does not claim an independent external biological annotation validation.

## Limits and date confounding

Flight animals were euthanized and dissected on 14 January 2021, both ground arms on 17 January 2021. Flight contrasts are perfectly confounded with this three-day collection-date difference. The global ground–vivarium association shares a collection date and shows that control context also differs, but cannot identify a flight mechanism. Exact label enumeration does not adjust away date, cell composition, sequencing measurement, or unrecorded cohort differences. These animals and measurements are reused across the binary and multiway retrospective analyses; Holm over this new three-test family does not provide a combined multiplicity correction for all analyses previously tried. Canonical annotated filtering excludes novel and poorly covered junctions. Boundary clusters can overlap, share genes, and are local choices rather than full transcript reconstructions. Avoid complete-isoform, reading-frame, protein, injury, cell-intrinsic, radiation-mechanism, or causal-flight language.

No outreach was performed.

## Artifacts and offline reproduction

Figure MW1 is supplied as PNG (300 dpi), vector PDF, and editable SVG, with a full caption. Complete animal-level compositions, all-arm summaries, coverages/entropies, candidate gates, definitions, overlap, exact distributions, complete local tests, sensitivities, omissions, and independent checks are retained. Run the scripts in order using the README commands, with --base pointing at the frozen input directory. Outputs resolve to the script-relative multiway folder. The original 32 MB GTF is unnecessary for offline reproduction because the exact lookup and source hash are bundled. No network call is needed. The manifest hashes every deliverable except the manifest itself.

Python 3.12.14, NumPy 2.3.5, pandas 2.2.3, SciPy 1.16.3, Matplotlib 3.10.7.
