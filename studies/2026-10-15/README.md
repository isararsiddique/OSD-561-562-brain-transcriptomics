# Two spaceflight brain studies: analysis code and reproducibility

This release contains the analyses supporting two current studies, with separate scientific questions and inferential units. It provides the source code, configuration, methods, small result tables, figures and provenance needed to understand the work. NASA measurement datasets and large generated outputs are distributed through the study archives rather than committed to GitHub.

| Study | Question | Start here |
|---|---|---|
| Paper 1: **Collection protocol outweighs spaceflight effects on brain immediate-early genes and resists post hoc correction** | How large are collection-protocol effects relative to matched flight effects, and when can correction amplify error? | [Study guide](paper1/README.md), [methods and workflow](paper1/METHODS_AND_WORKFLOW.md) |
| Paper 2: **Cerebellar RNA junction profiles show comparator dependence and limited transport across spaceflight missions** | How do comparator choice, exact-test resolution and cross-cohort transport affect junction-profile interpretation? | [Study guide](paper2/README.md), [methods and workflow](paper2/METHODS_AND_WORKFLOW.md) |

[Dataset links and provenance](DATA_SOURCES.md) identify every primary and external NASA accession. The studies reuse some deposits but answer different questions; their analyses and conclusions should not be pooled.

## What is included

Scientific source files are copied byte-for-byte from the 15 October 2026 analysis archives. The directory layout is preserved so scripts retain their original relative paths. Supporting figure builders, independent verifiers, mathematical checks, dated plans, sensitivity reports and small aggregate outputs accompany the core pipelines.

The [source manifest](SOURCE_FILE_MANIFEST.tsv) records the size and SHA-256 of every archive member, whether it is included here or remains in Zenodo, and the reason for that decision. Input directories, acquired datasets, cached source material, compressed matrices, simulation/bootstrap arrays, ZIP archives and generated files larger than 2 MiB are omitted from this GitHub release. All scientific scripts are retained regardless of size.

Older component reports can retain earlier working titles, analysis counts and dated conclusions. They document the analysis history. The two study guides above describe the current articles. Superseded manuscripts and cover letters are not part of this release.

## Data and archive access

| Study | Archive | Reserved DOI |
|---|---|---|
| Paper 1 | `Paper1_Analysis_Code_and_Results_Zenodo_20261015.zip` | `10.5281/zenodo.23388691` |
| Paper 2 | `Paper2_Data_Code_Results_20261015.zip` | `10.5281/zenodo.23387329` |

**Deposit status:** both Zenodo records are drafts with uploads in progress at release preparation. Reserved DOIs become resolvable when the records are published. See [ARCHIVES.json](ARCHIVES.json) for archive hashes and current recorded status. Public source datasets are already accessible through the NASA links in [DATA_SOURCES.md](DATA_SOURCES.md).

## Reproduce locally

1. Clone this repository and obtain the required complete study ZIP from the corresponding Zenodo record when public, or use the local verified release ZIP.
2. Run the restoration helper from the repository root. It verifies the complete ZIP SHA-256, rejects unsafe paths and restores the exact archive layout into a local directory excluded from Git.

```sh
python studies/2026-10-15/restore_archive.py --paper paper1 --archive /path/to/Paper1_Analysis_Code_and_Results_Zenodo_20261015.zip
python studies/2026-10-15/restore_archive.py --paper paper2 --archive /path/to/Paper2_Data_Code_Results_20261015.zip
```

When a deposit is public, `--download` can obtain its ZIP through the public record API. No login token or browser cookie is used.

```sh
python studies/2026-10-15/restore_archive.py --paper paper1 --download
python studies/2026-10-15/restore_archive.py --paper paper2 --download
```

3. Follow each study's workflow using a separate Python 3.12 environment. Full reproductions use the restored package under `.local-reproduction/paper1/` or `.local-reproduction/paper2/`. The copied GitHub scripts have the same scientific content as those archive scripts.

Do not run full-data validators against the intentionally incomplete GitHub-only folders. They require the complete archive. Re-running analyses can replace result files and invalidate the archive seal; preserve the supplied outputs and perform new fits in the local reproduction copy.

## Verify this GitHub release

```sh
python studies/2026-10-15/validate_release.py
```

This checks every distributed source file against its archive SHA-256, compiles Python syntax without executing analyses, and verifies that the excluded input directories and binary dataset formats are absent. These packaging checks are separate from the scientific checks recorded in each study archive. They do not constitute a new end-to-end scientific reanalysis.

See [RELEASE_VALIDATION.md](RELEASE_VALIDATION.md) for the checks actually run, including the 11 passing existing statistical regression tests and the exact environment used.

## Interpretation and authorship

Both studies are retrospective secondary analyses. They generate no new animal or sequencing data and do not repeat upstream raw-read alignment. Animals, rather than technical libraries, junctions, genes or overlapping holdout splits, define biological replication. Collection and exposure confounding, finite permutation resolution and measurement coverage are substantive limits, documented in the study workflows.

The current manuscript disclosures state that AI tools were used only to improve writing and grammar. Authors reviewed and edited the resulting text and take responsibility for the scientific content. The original code, results and dated audit records are preserved as research provenance.

Paper 2's existing archive specifies MIT for code and CC BY 4.0 for derived tables and documentation. Paper 1's licence is awaiting author confirmation; proposed Zenodo licence metadata does not itself constitute publication. Upstream datasets and third-party source code retain their original terms.
