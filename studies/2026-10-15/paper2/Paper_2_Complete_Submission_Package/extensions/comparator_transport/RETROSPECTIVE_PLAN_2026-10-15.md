# Comparator-specific frozen transport and RR10 external cohort

Frozen15October2026 before reading RR10junction outcomes or calculating new comparator-specific validation scores. Previously known: RR23 FLT–GC-trained binary and multialternative scores did not support positive transport in RRRM2. The RR10 metadata show distinct source animals, 5FLT+5GC, but a shared SpaceX21 launch with RR23, femaleB6129SF2/J rather than maleC57BL/6J, shared marker injections and different collection/dissection procedures. RR10 exposure is again aliased with collection dates. It is an external distinct-animal cohort, not independent-launch replication or causal isolation.

## Fixed question

Does the discovery ground comparator change the direction of the frozen junction profile's external transfer? Existing validation trained only on FLT–GC. Test FLT–GC and FLT–VIV native-composition directions separately in both RRRM2 and RR10, retaining the previously unsuccessful result. No choice of comparator or blend based on validation outcomes. Do not optimize weights against either validation cohort.

## Measurements and inference

Use the original923cluster definitions and every frozen selected alternative, counts and canonical boundary denominators under identical annotation/coordinate filters. Each validation cluster requires at least20selectedreads and90%dominance in every included animal. No gate tuning. Cluster membership selection uses counts without exposure labels.

For comparator C in GC,VIV, freeze d_c=mean_RR23_FLT,p−mean_RR23_C,p and u_c=d_c/||d_c|| from RR23 alone; zero direction retained if norm<=1e−12. Score each validation animal as equal-cluster mean of u_c'*(p_validation−mean_RR23_C). Frozen reference centering cancels from contrasts. Primary estimands remain RR23trained directional shifts in each cohort's measurable overlap, not shared biological features across all three cohorts.

RR10: all10eligible source animals, one5/5allocation block, all252assignments. RRRM2: previous19eligible animals, three route-by-age strata and11200assignments, same exclusions. Whole-animal exact positive one-sided tests, plus absolute two-sided values. The four directional global tests (2comparators ×2externalcohorts) constitute ONE new retrospective Holm family atalpha.05. Report all four raw and adjusted values even if none passes. This newfamily does not create a project-wide confirmatory guarantee or retrospectively prespecify earlier work.

Bootstrap4000resamples within biological group/stratum; primary intervals keep frozen RR23directions. Include a two-cohort sensitivity that re-estimates directions from RR23resampled animals. Intervals approximate, pointwise, not proof of equivalence on failure. Seed2026101541+cohortindex. Local tests across all923clusters, p=1 for unmeasurable; jointBY/BH/Holm across ALL4×923=3692newlocal slots. Original tests remain unchanged. RR10's minimum one-sidedp1/252 implies local Holm rejection impossible; do not use zeroHolm as biological evidence.

Export source hashes, all animal scores, all four global null distributions, all gate failures, all3692local results and preselectedTcaf1 measurements. Tcaf1 is descriptive within that family and not a separately rescued candidate. Nativeeffects in RR10 andRRRM2 have different animal/mission contexts; no averaging into a universal causal effect. If no common measurable clusters, report measurement failure without an invented score.

## Additional fixed diagnostics

Report both comparator score effects on the intersection of clusters measurable in BOTH validation cohorts, as a separate complete four-test Holm sensitivity. This uses label-free measurement overlap and tests comparability without replacing the primary each-cohort overlap. Report magnitude/direction differences without selecting the more favorable panel. No new benchmark or newly chosen feature set based on these results.
