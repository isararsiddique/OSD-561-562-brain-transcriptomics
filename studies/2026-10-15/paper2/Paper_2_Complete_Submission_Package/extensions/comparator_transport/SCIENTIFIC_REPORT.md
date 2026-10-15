# Comparator-specific transport across three cerebellar cohorts

This extension adds real publisher-verified RR10data from10distinct source animals and tests both discovery ground-comparator directions in both external cohorts. It follows a plan frozen before the new outcomes, while acknowledging that the previous RRRM2 FLT–GC-trained result was already known. RR10 was selected on tissue, source-animal and processed-data criteria, not an observed effect. The total55animals comprise26RR23,19eligibleRRRM2 and10RR10 animals across two launches, not three independent missions.

## Complete global outcomes

| cohort | direction | panel | n_clusters | n_animals | n_FLT | n_GC | n_assignments | effect_pp | conditional_CI95_low_pp | conditional_CI95_high_pp | two_cohort_CI95_low_pp | two_cohort_CI95_high_pp | exact_p | exact_two_sided_p | Holm_four_global_p |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| RRRM2 | FLT_minus_GC | cohort_measured | 281 | 19 | 10 | 9 | 11200 | -0.2770039 | -0.9635373 | 0.3777646 | -1.0630885 | 0.7392302 | 0.73625 | 0.521875 | 1.0 |
| RRRM2 | FLT_minus_GC | intersection_both_cohorts | 270 | 19 | 10 | 9 | 11200 | -0.407354 | -1.0711686 | 0.2158823 | -1.1029356 | 0.6645014 | 0.8396429 | 0.3225893 | 0.9166667 |
| RRRM2 | FLT_minus_VIV | cohort_measured | 281 | 19 | 10 | 9 | 11200 | 0.1823613 | -0.5070325 | 0.8072152 | -0.7197008 | 1.4658688 | 0.3495536 | 0.6982143 | 1.0 |
| RRRM2 | FLT_minus_VIV | intersection_both_cohorts | 270 | 19 | 10 | 9 | 11200 | 0.0728082 | -0.5911246 | 0.7045436 | -0.7707445 | 1.4362167 | 0.4400893 | 0.8757143 | 0.9166667 |
| RR10 | FLT_minus_GC | cohort_measured | 702 | 10 | 5 | 5 | 252 | 0.4357856 | -1.0133282 | 2.0329407 | -0.9991941 | 1.890736 | 0.3134921 | 0.6269841 | 1.0 |
| RR10 | FLT_minus_GC | intersection_both_cohorts | 270 | 10 | 5 | 5 | 252 | 0.5026652 | -0.9936677 | 2.1455625 | -0.9493972 | 1.8940904 | 0.3055556 | 0.6111111 | 0.9166667 |
| RR10 | FLT_minus_VIV | cohort_measured | 702 | 10 | 5 | 5 | 252 | 0.0596201 | -1.1054978 | 1.2389795 | -0.9561099 | 1.3182304 | 0.4642857 | 0.9285714 | 1.0 |
| RR10 | FLT_minus_VIV | intersection_both_cohorts | 270 | 10 | 5 | 5 | 252 | 0.3865709 | -0.5373376 | 1.3464868 | -0.6823687 | 1.5145795 | 0.2261905 | 0.452381 | 0.9047619 |

No primary or intersection-panel global test passes its four-test Holm family. RR10 measures702/923clusters (76.06%), compared with281/923 (30.44%) in RRRM2. Thus unsuccessful transfer is also observed in a cohort with substantially larger measurement overlap. The RR10 FLT–GC-trained effect is+0.4358pp,95%conditional interval[-1.0133,2.0329],exact positive-transferp=.31349. FLT–VIV-trained direction does not rescue transfer in either cohort. These nonsignificant results do not prove a universal effect is absent or identify why cohorts differ. The intervals preserve uncertainty.

## Design and causal limits

RR10/OSD-563 is female wild-typeB6129SF2/J,5FLT/5GC,14–15weeks at launch, with marker injections shared across arms. RR23 is maleC57BL/6J. Both were launched on SpaceX21; their source animals/RFIDs do not overlap. RR10 flight animals were euthanized02January2021 and controls05January; carcass dissections were12May and14May. These context differences and exposure/date aliasing prevent causal interpretation. All10STARdata files and the processing ZIP match publisher MD5; software-version sidecar has a malformed33-character digest that is preserved as a provenance exception. Verified processing commands and frozen reference metadata establish coordinate compatibility. See the complete cohort-search record for inclusions/exclusions.

RRRM2 is an independent-launch cohort with19eligible animals and the previously documented full-stratum exclusion. Its eligibility and stratification were unchanged. RR10 provides distinct-animal external evidence, not new independent-launch replication.

## Exact frozen analysis

All923original clusters and selected alternatives remain frozen. Label-blind measurement gates are evaluated separately in each external cohort. For comparatorC, d_c is the RR23 FLT−C mean native composition contrast, u_c=d_c/||d_c||. Each validation animal receives u_c'*(p_i−mean_RR23_C); the common reference centering cancels from group contrasts. The score averages equally over the cohort's measurable clusters. No validation outcome enters training weights or comparator choice. Four global one-sided tests form one new retrospective Holm family; the270cluster intersection is a separate complete four-test sensitivity.

RR10 uses252exact allocations. RRRM2 uses11200allocations preserving three route-by-age strata and weights each stratum by its fraction of all19animals. Ties use1e−12tolerance; absolute two-sided tests are also retained. The new local family contains4×923=3692slots, assigningp=1 to unmeasurable clusters, with jointBH/BY/Holm. NeitherBHnorBY yields a local discovery. For Holm, even the finest p minimum1/11200 exceeds.05/3692; the first adjusted minimum is.32964. Holm rejection is therefore impossible throughout the joint family, and that zero count is not biological evidence.

Four thousand bootstrap draws resample biological animals within group/stratum; primary intervals hold RR23directions fixed. A separate two-cohort sensitivity re-estimates directions from RR23animal resamples without fitting to validation outcomes. All bootstrap weights and draws are saved. Intervals are approximate/pointwise and nonsignificance is not equivalence.

## Tcaf1 and complete outputs

| cohort | direction | effect_pp | conditional_low_pp | conditional_high_pp | exact_p | joint_BY |
| --- | --- | --- | --- | --- | --- | --- |
| RRRM2 | FLT_minus_GC | -5.3892882 | -14.4306871 | 4.0930831 | 0.8594643 | 1.0 |
| RRRM2 | FLT_minus_VIV | -5.1251633 | -14.0482445 | 4.2790756 | 0.851875 | 1.0 |
| RR10 | FLT_minus_GC | -0.3372012 | -6.7297755 | 5.8459456 | 0.5357143 | 1.0 |
| RR10 | FLT_minus_VIV | -0.632568 | -6.7851672 | 5.0526719 | 0.5595238 | 1.0 |

Tcaf1 is measurable in both external cohorts but neither frozen direction supports a transferred association. Its complete fourtest row family remains part of all3692tests, not a selected standalone test. Candidate intervals are descriptive and unadjusted for selection.

## Verification

218independent checks reconstruct every measurement, native direction, global and local exact allocation result, all64000bootstrap effect draws, intervals and multiplicity arithmetic from frozen count matrices/metadata. The verifier does not import analysis functions. All sources, unfavourable outcomes, animal scores, gates and preanalysis decisions are retained. Run scripts/analyse.py, scripts/verify_independently.py, scripts/build_report_figure.py. These checks establish numerical consistency, not causal identification or successful replication.
