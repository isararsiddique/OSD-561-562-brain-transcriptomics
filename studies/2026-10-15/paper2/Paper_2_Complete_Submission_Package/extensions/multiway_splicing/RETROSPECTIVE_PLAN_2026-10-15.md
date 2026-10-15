# RR23 multiway junction-composition analysis: retrospective plan

Written 15 October 2026 before loading or testing multiway cluster outcomes. This additional analysis family extends an existing retrospective binary-junction study of the same RR23 cerebellar animals. It is not preregistered, independent confirmation, or a new cohort.

Author: Isarar Siddique. Affiliation and email: blank.

## Frozen measurement rules

Use the existing canonical, annotated all-animal junction-count matrix (805,562 genomic junctions and 26 collapsed animals) and its animal metadata. Keep the animals as the independent units. No library is treated as a biological replicate.

Separately group genomic junctions by chromosome, strand, and shared intron start or shared intron end. Strand is retained explicitly; these are genomic start/end boundaries, not claims about transcript isoforms. A junction qualifies if it has at least five reads in at least three of all 26 animals. A multiway cluster must have at least three qualifying alternatives. Retain all qualifying alternatives in that cluster, never the most associated subset. Require the sum of retained alternative counts to be at least 20 in every animal and at least 90% of ALL canonical annotated junction counts at the same boundary in every animal. The denominator includes nonqualifying alternatives. No group labels are used for these steps.

For the retained alternatives, compute each animal's composition. Using the unweighted mean of compositions across all 26 animals, require every retained alternative to have mean proportion at least 0.01, at least three alternatives to have mean proportion at least 0.05, and the largest mean proportion to be at most 0.90. These pooled rules remove effectively binary or degenerate multiway clusters without using group contrasts. Selection is frozen once; it is not repeated inside permutations or animal-omission checks. Report the full selection funnel and all exclusions.

## Primary global comparisons

Use sqrt(composition) locally. For cluster e, H_e is the 26 by k_e matrix of square roots of proportions; K_e = H_e H_e^T is PSD. The equal-cluster animal kernel is K = mean_e K_e, so a gene with several clusters can contribute repeatedly. Document every reused junction and intersecting cluster pair.

For FLT versus GC (9 versus 9), FLT versus VIV (9 versus 8), and GC versus VIV (9 versus 8), use the squared centroid distance in this kernel, T = d^T K d, where d is +1/n_A for group A and -1/n_B for group B. Report sqrt(T/2) as the mean-cluster Hellinger centroid distance (centroids in sqrt space, not compositions renormalized after averaging). Enumerate every fixed-size assignment: 48,620 for 9 versus 9 and 24,310 for each 9 versus 8. Count greater or equal statistics, including the observed assignment; exact p = count / total, with a documented floating-point tie tolerance. Holm controls the three primary global p values at 0.05. A plot based on centered K is descriptive and will show explained variance and every animal.

## Complete local exploratory family and characterization

For every selected cluster and all three comparisons use the same exact centroid statistic and complete fixed-size label enumeration. Report raw p, BH and BY q over the single full family of 3E tests, plus per-contrast BH/BY for transparency. BY over 3E is the conservative local discovery criterion at 0.05. No pathway or candidate selection is promoted to independent confirmation. Include native-composition group means, differences, each animal's composition and coverage, alternative count, dominance, entropy, exact coordinates, and local effect sizes for every cluster, including null results.

## Annotation and sensitivities

If the cached Mus_musculus.GRCm39.112.gtf.gz is available, reconstruct transcript introns from adjacent exons (upstream exon end+1 to downstream exon start-1, chromosome and strand matching). Intersect matching gene sets over EVERY retained alternative in a cluster. Assign a gene only if this intersection contains exactly one gene; report absent and ambiguous mappings without invented names. On this uniquely mapped subset, assign each gene equal total weight, dividing each gene's weight across all its clusters. Perform the same exact global comparisons and Holm correction within this explicitly secondary three-test family. Report the mapped-subset equal-cluster comparison to distinguish subset from weight changes. A 95% dominance sensitivity and leave-one-animal-out effect ranges use the original frozen clusters. They describe robustness in the same data and do not supply independent validation; no additional hypothesis family is formed from omission results.

## Independent verification and interpretation

Check PSD eigenvalues, equal diagonal entries, identities connecting feature and kernel statistics, direct permutation enumeration for global statistics, selection invariance to animal-order/group-label changes, permutation counts and observed-assignment inclusion, and local multiplicity calculations. Save input/output hashes, scripts, complete tables, verification report, scientific report, and a legible four-panel figure.

Flight animals were euthanized and dissected 14 January 2021, both ground groups 17 January 2021; flight contrasts are perfectly confounded with this date. Junction composition describes local splice choices in bulk cerebellum; it cannot resolve complete isoforms, proteins, cell-intrinsic changes, injury, a radiation mechanism, or a causal flight effect. Shared boundaries can overlap and are not independent molecular events. Read-based coverage and dominance filters restrict inference to well-covered canonical annotated clusters and exclude rare/novel alternatives. Exchangeability of nonrandom cohort labels remains an assumption even with exact arithmetic.
