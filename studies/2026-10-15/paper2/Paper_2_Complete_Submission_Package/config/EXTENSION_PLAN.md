# RR23 three-arm RNA-processing extension

Written 14 October 2026 before fitting the new vivarium junction contrasts. This is a retrospective extension, not a prospective preregistration or an independent validation of previously inspected RR23 results.

Objective: test whether the pattern of local junction usage distinguishes flight, habitat ground control and vivarium control samples, and quantify how the control comparison affects interpretation of earlier flight–ground findings.

The earlier paired-module draft is preserved. Stronger prose is contingent on stronger evidence; significant findings and unsuccessful hypotheses will both be reported.

## Fixed design and measurements

- OSD-525: 26 deposited source animals, 9 flight, 9 habitat ground control, 8 vivarium control; 29 technical libraries collapsed to source animals.
- All animals are male C57BL/6J. Flight specimens were collected on 14 January 2021 and both ground groups on 17 January. Date and flight status are aliased. Shared husbandry/preservation fields do not eliminate that confounding.
- Use the historical 6,277 RR23 local binary junction definitions, selected using pooled measurement coverage rather than differential-expression significance. No positive-event subset is used to choose the primary features.
- Retain annotated canonical unique STAR junction counts, strand 1/2, intron motif 1–6, maximum overhang at least 8. Sum technical libraries before forming animal-level fractions.
- Primary common measurement gate: both selected junctions have total at least 20 unique reads in every one of the 26 animals, and contribute at least 90% of all retained junctions sharing the boundary in every animal. Failed definitions remain p=1 in event correction families. Report a 95% dominance sensitivity.
- Local junction fractions are not full isoforms, intron retention, productive protein, cell composition or functional injury.

## Analysis families

1. Three primary global tests: FLT–GC, GC–VIV and FLT–VIV. Statistic: mean squared difference in animal-mean fractions across all common measured definitions. Null: animal vectors are exchangeable within each pairwise comparison. Enumerate all 48,620 / 24,310 / 24,310 binary allocations. Holm across three tests. Report raw squared contrast and an unbiased squared-distance estimate; negative unbiased estimates are retained.
2. Event-level characterization: the same three contrasts and exact label spaces. BH across all 3 × 6,277 = 18,831 tests, with within-contrast BH and BY as explicit sensitivities. Magnitude threshold of five percentage points is descriptive, not a measured functional threshold.
3. Secondary contrast classification: require FLT–GC and FLT–VIV differences in the same direction, and GC–VIV equivalence within an operational ±5-percentage-point bound using approximate Welch TOST. Conjunction p is the maximum of the three required p values, followed by BH across 6,277 definitions. This is exploratory and its distributional/equivalence assumptions must be disclosed. Non-significant GC–VIV differences alone do not establish flight specificity.
4. Shared-control dependence diagnostic: compare the FLT–GC and GC–VIV effect profiles over the complete common measured set, then repeat across all 126 disjoint 4/5 partitions of GC animals. Do not test correlated junctions as independent observations or call anticorrelation a mechanism. Animal bootstrap intervals preserve within-animal feature dependence.
5. Technical sensitivity: 10 independently seeded binomial thinnings of the complete animal/junction count matrix to the lowest observed total retained-junction depth. All events sharing a junction use the same thinned count. Recompute fractions and global tests on the fixed primary measured set; report every replicate. This is a measurement-depth sensitivity, not biological replication.
6. Existing cross-cohort validation is separately labeled as previously executed evidence. Verify saved multiplicity/count calculations and input provenance; do not rename same-launch RR10 as an independent-launch replication.

The validity of permutation inference depends on exchangeability. Neither exact tests nor adjustment for RIN or depth can identify a flight effect separately from collection date in the observed design.

## Decision

If a reproducible global processing association and interpretable control comparison emerge, develop a paper around context-dependent local RNA processing and its transportability. If not, report the failed extension and retain the corrected baseline. Do not claim a common causal spaceflight mechanism or a new treatment target from this extension alone.
