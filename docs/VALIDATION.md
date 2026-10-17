# Executed validation: 17 October 2026

[Study overview](../README.md) · [Executed run summary](../results/verification/publication_run_summary.json) · [Result comparisons](../results/verification/portable_result_comparison.json) · [Source checks](../results/verification/live_source_validation.json)

## Publication checkout and fresh execution

The cleaned repository's single entry point reconstructed all 303 input files in a separate folder using a verified local cache plus the published frozen definitions. The complete 62-stage analysis then ran in another separate output folder, both with spaces in their names. This tests the data-free checkout without dependence on the historical analysis folders. The default analysis run prohibited network access.

| Check | Executed result |
|---|---:|
| Complete pipeline stages | 62 / 62 |
| Input SHA256 checks | 303 / 303 |
| Original numerical/source checks | 56 / 56 |
| Integrated numerical checks | 134 / 134 |
| Figure-value checks | 91 / 91 |
| Deeper model/source checks | 650 / 650 |
| Total numerical checks | 931 / 931 |
| PDF single-page / embedded-font / text-bound checks | 19 / 19 |
| Input/acquisition safeguard tests | 8 / 8 |
| Representative live primary-source input downloads | 7 / 7 |
| Complete prior tables compared | 18 |
| Rows in those comparisons | 258,677 |

The prior and current complete tables match numeric values at relative tolerance 1e-8 and absolute tolerance 1e-10; categorical identities also match. Comparisons include all 50 actual additional-cohort deletion refits, both engines and expanded correction families. This is computational verification, not an extra biological replication.

## Acquisition and source integrity

Representative live checks included a STAR bulk count file, the CA GeoMx matrix and its duplicate source copy, RR8 barcodes, mouse Hallmark sets, author cell metadata extracted from the 3.93 GB RR3 archive, and Visium coordinates extracted through verified ZIP local headers/byte ranges. Every downloaded file matched its published length and SHA256; extracted members passed CRC32. These seven checks do not assert that every endpoint was retested live.

Eight safeguard tests verify valid paths with spaces, rejection of same-length count edits, duplicate/escaping input paths, changed downloaded contents, ignored ZIP byte ranges, corrupt member CRCs and source-path escape. A changed source does not become an accepted input or a complete bundle.

One generated download ledger was converted to portable relative paths; all assay bytes and biological inputs are unchanged. The public input manifest defines the published bundle. Frozen original analyses remained preserved during the update.

## Figures and document integrity

The workflow was visually inspected, with clear question–data/method–evidence lanes and an explicitly proposed laboratory validation step. The 19 PDFs have embedded fonts and valid text bounds; integrated headings are black. The complete gallery has 21 PNG previews and 13 editable SVGs. Displayed numbers are verified against executed tables. Final journal column sizing remains a separate layout task.

The repository integrity check validates published artifact hashes, script hashes, local Markdown/image links, all 303 source records, absence of assay/model objects and GitHub per-file size limits. It is available with `python tools/verify_repository.py` and runs in the repository CI with the safeguard tests.

[All 62 executed stage logs](../results/verification/stage_logs/) are included with portable paths for inspection.

## Practical limits

The full scientific execution was tested on macOS arm64 using Python 3.12, R 4.5.3, limma 3.66.0 and edgeR 4.8.2. Linux/macOS/WSL instructions use configurable paths and standard dependencies, but a complete Linux or Windows run was not directly tested. Repository CI checks standard-library safeguards and publication integrity; it does not rerun the entire biological analysis.

Correct computations do not resolve preservation/day confounding, identify a causal mechanism or supply independent experimental replication. Complete families, uncorrected results and missing coverage remain visible. See [limitations](LIMITATIONS.md).
