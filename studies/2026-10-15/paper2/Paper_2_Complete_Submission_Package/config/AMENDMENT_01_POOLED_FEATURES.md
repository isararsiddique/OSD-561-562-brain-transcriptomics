# Amendment 1 — label-invariant all-animal feature construction

14 October 2026. The initial extension used the historical 6,277 definitions constructed from the 18 flight/habitat samples. Although that construction did not select on differential p values, it used a subset whose membership depends on the groups being compared with vivarium. Conditional randomization of those comparisons therefore cannot assume a label-invariant selection procedure. Treating the historical panel as prospectively independent would be incorrect.

The first-pass results are preserved in `results/legacy6277_first_pass/` and are not the primary evidence of the beta manuscript. This amendment is retrospective and is explicitly recorded after inspecting those preliminary comparisons.

The repaired primary definitions are constructed from pooled counts of **all 26 animals**, without using group labels:

1. Retain a canonical annotated unique junction with at least five reads in at least three of the 26 animals.
2. Form groups by chromosome, strand and shared start or end boundary. Require exactly two qualifying alternatives and sort their identifiers lexically.
3. Require at least 20 pair reads in every animal and pooled mean first-junction fraction strictly between 0.05 and 0.95.
4. Apply the 90% dominance gate against all canonical annotated junctions sharing the boundary in every animal. Failed dominance definitions remain p=1 in the event family. A 95% sensitivity is retained.

All feature-selection operations are invariant under every allocation of the exposure labels. The number of qualifying definitions is determined by these fixed measurement rules, not a chosen outcome. Primary global and event tests then use the same all-animal panel. Family sizes are M definitions per contrast, 3M jointly and M for the exploratory conjunction, where M is reported from the prepared inputs.

Gene assignments use the intersection of GRCm39.112 transcript-intron gene mappings for the two exact junctions. Ambiguous or absent gene mappings remain explicitly unassigned.

The remaining extension plan is unchanged. Historical 6,277-result counts and earlier independent-launch validations are separately described as prior work, using their original definitions. They are not fresh validations of the newly formed all-26 panel. Exact test validity still requires exchangeability, and flight remains confounded with collection date.
