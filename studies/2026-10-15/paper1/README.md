# Collection protocol outweighs spaceflight effects on brain immediate-early genes and resists post hoc correction

**Authors:** Isarar Siddique, Raiyan Siddique, Arshi Niyaz, Jamyang Ugyen Tshomo, Atish Raza, Juan Pablo de Rivero Vaccari and M. Windy McNerney.

The study asks whether collection-protocol effects can rival matched spaceflight effects in mouse brain expression, and whether removing a strong collection-response score can worsen prediction through noise amplification. It combines animal-level exact inference, metadata/protocol decomposition, correction geometry, holdouts and gain-constrained regression. The seven primary deposits represent 83 mice from three missions; four additional datasets evaluate the prediction method.

## What we found

In ground mice that never flew, the RRRM-2 ISS-terminal protocol exceeded live-animal-return collection on the rapid primary-response panel by approximately 0.70 log2 units in cerebellum and 1.12 in hippocampus. Matched flight effects varied across mission/comparator settings. Housing can also contribute: collection protocol is not the only source of the signal.

A collection-response sensor can discriminate cohorts while its genome-wide loading remains too noisy to estimate reliably. Score subtraction had spectral gain around 8–15 and worsened held-out prediction in two of three primary settings. A convex estimator constrained to gain at most two improved held-out error by about 1.2–7.2% across seven settings. This is prediction evidence, not a demonstration that correction recovers causal biology or universally outperforms every comparator.

## Navigate the analysis

| Component | Directory |
|---|---|
| Core loaders, normalization, panel contrasts, exact inference and score correction | `Paper_1_Complete_Submission_Package/Paper1_Strengthened_Review_20261015/scripts/` |
| Protocol/flight decomposition and dominance family | `Paper_1_npj_LaTeX_Build_20261015/research/D_protocol/` |
| Leakage-controlled ground holdouts and exact discrimination reference | primary package `extensions/ground_holdout/` |
| Convex constrained estimator, matched-grid benchmarks and independent optimizer checks | primary package `extensions/strengthening_20261015/method_benchmark/` |
| Real-count-calibrated correction benchmark | primary package `extensions/strengthening_20261015/correction_benchmark/` |
| Count models, paired biology, interactions and biological sensitivity | primary package `scripts/` and `extensions/` |
| Identical-loading noise-ceiling diagnostic | `Paper_1_npj_LaTeX_Build_20261015/research/B_novelty/analysis/` |
| Hashed published-signature plan and separately labelled post hoc analysis | `Paper_1_npj_LaTeX_Build_20261015/research/F_signatures/` |
| Independent numerical audit and reimplementation | `Paper_1_npj_LaTeX_Build_20261015/research/A_audit/` and `research/G_independent/` |
| Reference verification and dated literature searches | `Paper_1_npj_LaTeX_Build_20261015/refs/` and `research/H_novelty_systematic/` |
| Figure builders and small outputs | `Paper_1_npj_LaTeX_Build_20261015/figures/`; primary package `publication_figures/`, `results/` and `supplementary_data/` |

[METHODS_AND_WORKFLOW.md](METHODS_AND_WORKFLOW.md) explains each calculation, its inferential unit, multiplicity family, validation and command order. [Dataset links](../DATA_SOURCES.md) and the [distribution manifest](../SOURCE_FILE_MANIFEST.tsv) identify the data and complete outputs stored outside GitHub.

## Reproduce

Restore the full archive using the release helper, then follow the commands in the workflow. The restored archive root is `.local-reproduction/paper1/`. Use Python 3.12 and the supplied dependency/runtime records. Full validators require the complete archive, including its frozen inputs and result tables.

Original component READMEs and historical reports can refer to earlier manuscript titles. They are retained as analysis history; this guide identifies the current study. The source scripts have not been rewritten to change the numerical results.

The Zenodo DOI `10.5281/zenodo.23388691` is reserved; publication is pending at release preparation. Paper 1's redistribution licence is awaiting author confirmation. NASA and third-party materials retain their original terms.
