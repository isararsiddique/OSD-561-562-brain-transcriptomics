# Results and interpretation

[Overview](../README.md) · [Methods](METHODS.md) · [Complete tables](../results/README.md) · [Limitations](LIMITATIONS.md)

## Original panel: effects agree more than significance decisions

The original bulk family contains 4,228 planned target/contrast hypotheses per engine: 3,584 estimable and 644 unavailable. DESeq2 identifies 17 corrected contrasts, comprising two flight, 14 collection-cohort and one age result. The flight results are lower young-terminal cerebellar AIF1/PYCARD RNA. EdgeR has 10 corrected contrasts; voom has none, with minimum q=0.07994. Young-terminal PYCARD is corrected by DESeq2 and edgeR (edgeR q≈0.031) but not voom. Therefore the association is method-sensitive, despite similar effect directions.

Across the 3,584 shared estimable panel contrasts, DESeq2–edgeR Spearman correlation is 0.968 and same-direction fraction is 0.941; DESeq2–voom values are 0.894 and 0.904. The Figure 2 whole-genome young-terminal display has a separate Spearman value of about 0.92. These compare different scopes and must not be substituted for one another. Original terminal preservation is associated with exposure and remains an alternative explanation.

The complete original GeoMx family has zero corrected discoveries under its primary specification; Welch sensitivity does not resolve the absence of corrected support. That does not establish equivalence or the absence of treatment effects. ROI/animal identity and slide confounding constrain a rescue claim.

Sources: [bulk family](../results/core/tables/bulk_panel_complete.csv), [method family](../results/extended/tables/bulk_method_comparison_panel.csv), [effect agreement](../results/extended/tables/bulk_effect_agreement_summary.csv), [GeoMx family](../results/core/tables/spatial_panel_complete.csv).

## RR8: reproducible RNA directions within a confounded design

Forty-four distinct mapped forebrain genes pass correction in all eight combinations of two engines and four cell-QC definitions. This is sensitivity agreement in eight shared mice, not eight replications. The article 5% gate gives forebrain voom effects of AIF1 −1.092 (q≈1.12e−8), PYCARD −1.256 (q≈1.47e−8), AIM2 +0.571 (q≈1.54e−5), NINJ1 −1.130 (q≈4.14e−7), and IL18 −0.776 (q≈7.21e−6).

The opposing sensor/adaptor RNA pattern is a candidate for further experiments. AIF1/PYCARD/AIM2 marker-relative HC3 tests remain corrected at approximately q=0.0296, while PYCARD's 500-UMI detection sensitivity is narrowly uncorrected (q≈0.0598). Non-panel reference genes may themselves change biologically. Exact age-blocked animal tests do not yield corrected discoveries in these complete families. Flight and ground libraries were processed on different days with perfect exposure/day confounding; computational sensitivity cannot disentangle that factor.

Sources: [eight-combination consensus](../results/deep/tables/forebrain_all_eight_QC_method_consensus.csv), [article-gate family](../results/deep/tables/article_5pct_panel_complete.csv), [Figure 7 source ledger](../results/extended/figures/figure7_source_ledger.json).

## Expression states: positive support and a meaningful failure

Primary shared-state PYCARD effects are negative and corrected: K=4 −0.867/−0.893 (voom/edgeR), and K=6 −1.192/−1.187. All eight donors meet the ≥50-cell eligibility threshold in the illustrated primary shared states. Some other states have poor overlap, including a K=6 state absent from one ground mouse.

The technical K=6 eligible state has effects −0.430/−0.451 and q=0.280/0.147, failing correction. Equal-state composition tests do not produce corrected confirmation, and primary K=6 all-state equal weighting is unavailable. This prevents a universal within-state or composition-independent claim. Labels are numerical partitions, not independently established cell identities.

Sources: [primary state family](../results/deep/tables/state_panel_complete.csv), [technical state family](../results/deep/technical_state_tables/state_panel_complete.csv), [displayed effects](../results/deep/tables/figure8_displayed_count_effects.csv).

## Added bulk cohorts: habitat/vivarium dependence

RR23 AIM2 flight versus habitat-ground has log2FC +0.740942, q=0.014624 (voom), and +0.732579, q=0.010846 (edgeR) in the expanded 1,812-test family. All 26 mouse-deletion refits remain positive and corrected in each engine; maximum q values are approximately 0.0270/0.0239. This is the strongest influence-stable additional-cohort target result.

Against vivarium, AIM2 has +0.208/+0.174 and q=0.843/0.898. RR23 PYCARD flight versus vivarium has −0.637/−0.630 and expanded q=0.0704/0.0794. Its smaller earlier 1,057-test family passed correction; that earlier value must not be presented as corrected in the expanded analysis. Every deletion is negative, but only 2/26 voom and 1/26 edgeR deletion fits pass expanded correction. Collection dates and habitat/vivarium differences qualify interpretation.

In MHU-3, SLC1A3 WT flight effects are about −0.214 with q=0.0172/0.0207; 22/24 voom and 23/24 edgeR deletions pass correction, all negative. Balanced-genotype S100B passes the full fit but only 4/24 and 7/24 deletions remain corrected. No corrected Nrf2 interaction establishes protection. All 24 deposited mice are used, with the published flight-KO count discrepancy documented.

Sources: [expanded family](../results/deep/tables/additional_panel_all_regions.csv), [all expanded deletions](../results/deep/tables/additional_mouse_deletion_expanded_complete.csv.gz), [deletion summary](../results/deep/tables/reported_expanded_deletion_summary.csv), [Figure 9 effects](../results/deep/tables/figure9_displayed_contrasts.csv).

## Regional transfer remains incomplete

OSD-613 has 18 different RRRM-2 animals without identity overlap with the original 29, but shares mission context. Its PYCARD balanced-age effects are −0.265/−0.264 and q=0.453/0.437. OSD-612 supplies nine cerebral measurements from existing RR10 animals: PYCARD +0.613/+0.436, q=0.660/0.933. One additional ISA flight sample lacks a deposited count column and is not imputed. These are unconfirmed regional tests, not demonstrations of absence or a pooled contradiction.

Sources: [animal-independence audit](../results/deep/tables/terminal_region_animal_independence.json), [missing count sample](../results/deep/tables/terminal_region_missing_bulk_samples.csv), [regional figure effects](../results/deep/tables/figure10_displayed_effects.csv).

## Spatial, promoter and pathway support has limits

RR3 Visium shows 883 corrected gene-by-section local spatial tests across 12 sections from six mice. The count concerns repeated gene/section hypotheses; it is not a count of independent genes or exposure discoveries. Exact six-mouse exposure tests over 20 assignments have zero corrected panel discoveries. The three eligible RR3 promoter-context types also have zero corrected discoveries across their complete 453-test family. Shared RNA/ATAC donors do not create independent replication.

At CAMERA fixed correlation 0.01, corrected counts are 156 core, 30 independent-bulk, 15 RR3 cell-type and 209 RR8 tests. Fixed 0.10 and estimated-correlation specifications have zero corrected discoveries. Estimated CAMERA changes residual degrees of freedom as well as correlation. Exact animal singscores do not confirm corrected pathway effects. The interferon-alpha score differences in RR3 and RR10 cerebellum are descriptive (exact p≈0.300/0.302; full-family q≈0.854).

Sources: [Visium summary](../results/extended/tables/Visium_summary.json), [promoter summary](../results/extended/tables/promoter_context_summary.json), [complete pathway distributions](../results/extended/tables/pathway_evidence_distribution_summary.csv).

## Orthogonal assays do not close the mechanistic gap

RR9 nCounter has 31 detectable panel targets and zero corrected exact discoveries; AIF1 is absent and PYCARD is below the detection gate. Submitted OSD-738 protein records contain four extreme ASC values; source agreement is confirmed, values remain flagged, and no invented repair or pooled RNA/protein effect is applied. These branches neither directly confirm nor falsify the full RR8 program.

Sources: [nCounter summary](../results/extended/tables/OSD536_analysis_summary.json), [protein audit](../results/extended/tables/OSD738_protein_source_audit.json).

## Working conclusion

The supported result is a context-dependent innate-immune RNA pattern whose evidence changes with cell QC, numerical state, control definition and tissue region. Stable RR23 AIM2 and multi-QC RR8 forebrain directions justify focused functional experiments. Preservation, processing day, small donor numbers, incomplete transfer and failed conservative checks prevent a causal or universal mechanism claim.
