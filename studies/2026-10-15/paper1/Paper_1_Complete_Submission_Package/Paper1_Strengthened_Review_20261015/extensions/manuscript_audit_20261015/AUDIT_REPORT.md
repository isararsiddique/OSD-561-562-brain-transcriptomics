# Independent manuscript and figure audit

15 October 2026. Read-only review of the saved scientific outputs and current manuscript, supplement, content builder, document builder, and numbered publication figure index. Shared scientific files were not edited. Final render QA is being performed by the parent task.

The independent checker passed **157 checks with no failures**. It imports no production statistics code. Its machine-readable report contains result summaries, each check, and SHA-256 hashes of the inspected CSVs. This certifies the checked internal numerical consistency, not causal identification or journal acceptance.

## Verified numerical claims

| Claim | Independently checked result |
| --- | --- |
| Data accounting | 129 libraries; 126 collapsed animal–region observations; 83 animals; 29 core animals, of whom 24 are regionally paired; 10 RR-10 paired animals; 7 location-conflict entries from 6 animals. |
| Direct hippocampal Circadian interaction | −0.405321754 log2; bootstrap interval [−0.590650834, −0.215082007]; wild p=0.000838; complete 32-test BH q=0.026816. Count-model panel mean −0.367146698. |
| Circadian metadata sensitivity | −0.459705032 log2; q=0.786272. The manuscript appropriately retains the loss of corrected support and sparse older-flight design support. |
| Circadian omission checks | Every full-data one-animal and one-gene omission remains negative; reported rounded ranges match saved outputs. |
| Primary gene interactions | 109 cerebellar and 7 hippocampal primary discoveries; 90 and 2 retain original global significance after metadata exclusions. |
| Holdout result | 362 overlapping splits from 30 unique ground source animals. Core sensor AUC=1 in every split. RR-23 mean sensor AUC=0.963889. |
| Original-score outside prediction | Mean MSE ratios are 1.090902, 0.961891 and 1.034881 for cerebellum, hippocampus and RR-23. The stated 9.1% and 3.5% deterioration follows these ratios. |
| Operator geometry | 15 aligned geometries; unbounded gains 8.129052–15.309384; bounded gains ≤2 within tolerance. |
| Joint score attenuation | All three outside-sensor bootstrap intervals cross zero; manuscript endpoints match the saved summaries. |
| Multiplicity | Independently recomputed BH for all 232 primary/sensitivity panels, 29 gradients, 104 paired panels, 32 direct-interaction panels, 54,989 reproducibility-screen entries and 72,023 finite gene-interaction tests. |
| Recorded numeric lookups | Every recorded unique row selector and numeric field agrees with its source CSV. |

The manuscript's main interpretation is supported: reliable within-resource collection-score discrimination does not guarantee transfer of genome-wide loadings. It correctly distinguishes expression attenuation, bounded numerical gain, biological recovery and independent replication. Circadian and beta-globin associations remain observational candidates. The absence of corrected paired or cross-cohort discoveries is not presented as equivalence or absence of biological effects. No new substantive biological overclaim was identified in the current text.

## Concrete editorial corrections

1. **Figure 2's word “Unadjusted” is ambiguous and scientifically inaccurate if read as covariate-unadjusted.** The displayed raw effects are estimated with age adjustment and pooled tissue adjustment. Use “Panel effects before ground-derived transformation,” and state that design nuisance adjustment is retained. Its probabilities also refer to a corrected family.
2. **Remove four obsolete figure pairs from the portable final artwork directory.** They are outside the authoritative final manifest but reuse earlier numbers: `figure6_ground_anchor_transfer`, `figure7_disjoint_anchor_partitions`, `figure8_known_truth_validation`, and `figure9_high_gain_domain_shift`. The current fourteen-item manifest and manuscript caption numbers are already consistent.
3. **Correct the supplementary electronic figure-source index.** It currently names `figures/*.svg`, `figures/*.png` and `results/figure_manifest.json`; the final numbered artwork is in `publication_figures/` and indexed by `results/publication_figure_manifest.json`. Keep the original source folder only as an explicitly labelled intermediate source.
4. **Make the supplementary reproduction order agree with the README.** The supplement's current sequence ends at `prepare_final_content.py` and tests. Add the final document builder and the supplied-equation/regeneration note. The README already distinguishes the superseded draft builder with `--figures-only` and records the final document build.
5. **Publisher-specific layout edits are in progress in the parent task.** The inspected pre-edit title is 17 words and abstract exceeds the supplied 150-word limit; the parent is also moving supplementary methodological content and the funding declaration according to verified journal guidance. Recheck the final rebuild after those edits.

## Authorship and declarations

The manuscript, current content builder, `authors.json`, cover letter and PDF metadata consistently use the five user-confirmed authors: Isarar Siddique, Raiyan Siddique, Arshi Niyaz, Jamyang Ugyen Tshomo and Atish Raza. The affiliation and correspondence agree with the pasted user record: Convolity AI Lab, Lucknow, India; isarar@convolity.com. ORCID 0009-0005-2305-9579 and the user-confirmed Atish Raza Data curation role are retained. The current authoritative files do not reproduce the old conflicting single-author identity.

No external funding, no competing interests, secondary-analysis ethics and substantive AI assistance are expressly declared. `all_author_approval` remains unresolved, as it should until confirmed. The permanent public revision archive and related-Paper-2 status are also openly unresolved release items. The review package should preserve these factual limits and not assert submission-ready status prematurely.

## Figure structure

The fourteen final figure entries are ordered 1–10 and S1–S4, and each caption matches its filename number. PNG metadata report at least 599 dpi; SVGs retain editable Arial text. Multipanel letter sequences match caption references: main Fig. 1 a–b, Fig. 3 a–b, Fig. 4 a–d, Fig. 5 a–c, Fig. 6 a–d, Fig. 8 a–b, Figs. 9–10 a–d; supplementary Fig. S1 a–b, S2 a–i, S3–S4 a–d. Figures 2 and 7 are single plots. Caption lengths stay below 350 words. These structural checks do not substitute for the parent's final visual inspection.

## Audit files

- `independent_checks.py`: standalone saved-result and figure-index checker.
- `independent_checks.json`: all 157 outcomes, numerical summaries, panel inventory and inspected input hashes.

Rerun from the analysis directory with the established scientific environment after final content changes. The checker only writes this audit directory.
