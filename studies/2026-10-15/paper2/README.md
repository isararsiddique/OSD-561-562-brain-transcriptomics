# Cerebellar RNA junction profiles show comparator dependence and limited transport across spaceflight missions

**Authors:** Isarar Siddique, Raiyan Siddique, Arshi Niyaz, Jamyang Ugyen Tshomo, M. Windy McNerney and Juan Pablo de Rivero Vaccari.

This study examines how comparator choice, finite exact-test resolution and cross-cohort transport change the interpretation of cerebellar junction profiles. It uses 55 distinct mice in three cohorts: RR23 (26), RRRM2 (19) and RR10 (10). RR23 and RR10 share SpaceX-21; the three cohorts span two launches.

## What we did and found

We constructed a label-invariant RR23 panel of binary junction choices and multialternative clusters, tested all three FLT/GC/VIV comparisons using complete animal-vector allocations, and assessed gene localization, measurement robustness, relative abundance and frozen transport.

The 5,249 primary binary choices and 923 multialternative clusters resolve all three group comparisons. The ground controls differ substantially, showing that the selected comparator changes the flight contrast. Dependence-robust analyses identify 47 gene associations and 41 clusters in the reported primary families. Technical read halves can preserve effect directions while changing discovery decisions; the exact BY arithmetic prevents a non-empty cluster discovery list shorter than 20 tests in the stated family.

Frozen directions reproduce in held-out RR23 animals, including a same-date housing axis, but do not establish positive transfer to either external cohort. Within RR23, flight is inseparable from collection date. Within-cohort prediction does not identify a causal spaceflight mechanism, and zero corrected local discoveries can reflect an unattainable testing threshold.

## Navigate the analysis

| Component | Directory |
|---|---|
| Complete primary offline analysis driver | `Paper_2_Complete_Submission_Package/scripts/reproduce_analyses.py` |
| Binary panel, primary exact tests and independent verification | primary package `scripts/` |
| Gene-level omnibus localization | primary package `extensions/gene_omnibus/` |
| Multialternative cluster construction/inference | primary package `extensions/multiway_splicing/` |
| Gene abundance and explicit equivalence | primary package `extensions/abundance_integration/` |
| RRRM2 binary and multialternative transport | primary package `extensions/independent_validation/`, `extensions/multiway_transport/` |
| Read-half stability, influence and dispersion | primary package `extensions/multiway_robustness/`, `extensions/effect_stability/` |
| Count-model and confounding sensitivity | primary package `extensions/count_model_benchmark/`, `extensions/confounding_sensitivity/` |
| Cohort eligibility, RR10 acquisition/provenance | primary package `extensions/cohort_search/` |
| Two-comparator, two-cohort frozen transport | primary package `extensions/comparator_transport/` |
| Transcript architecture | primary package `extensions/transcript_context/` |
| Exact BY floor simulation | `Paper_2_Revision_Build/research/D_simulation/` |
| Graded transport and five views | `Paper_2_Revision_Build/research/G_graded_transport/` |
| Held-out internal benchmark and same-date control | `Paper_2_Revision_Build/research/K_internal_transport/` |
| Enrichment and experimental-design calculator | `Paper_2_Revision_Build/research/I_enrichment/`, `research/H_design_calculator/` |
| Numerical audit, reference/novelty verification and figures | `Paper_2_Revision_Build/research/`, `refs/` and `figures/` |

[METHODS_AND_WORKFLOW.md](METHODS_AND_WORKFLOW.md) explains the calculation order, statistical families, resolution limits and interpretation. [Dataset links](../DATA_SOURCES.md) and the [source manifest](../SOURCE_FILE_MANIFEST.tsv) identify frozen data and complete outputs.

## Reproduce

Restore the complete ZIP to `.local-reproduction/paper2/`, create a Python 3.12 environment using `requirements.lock.txt`, and run `scripts/reproduce_analyses.py` from the primary package. The driver runs numerical analyses and independent verification without rerunning cohort searches or source downloads. Some added enrichment/reference stages require online services or separately obtained marker tables; the workflow distinguishes those from the offline primary reproduction.

All published code copied here is unchanged from the analysis archive. Full-data validators require the restored archive rather than the intentionally incomplete GitHub-only folder. The DOI `10.5281/zenodo.23387329` is reserved, with publication pending at release preparation. The existing archive specifies MIT for code and CC BY 4.0 for derived tables and documentation; upstream and third-party terms remain applicable.
