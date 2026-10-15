# Bounded cohort search and frozen RR10 data

Read `COHORT_SEARCH_REPORT.md`. This extension identifies and acquires one additional distinct-animal, same-launch cerebellar cohort (RR10/OSD-563), without examining outcomes. It does not establish another independent-launch replication or repair exposure/date confounding.

The parent can directly read `inputs/RR10/RR10_animal_junction_counts.csv.gz` and `inputs/RR10/RR10_animal_metadata.csv`. All ten source STAR files match publisher MD5 and have independently reconstructed matrix columns. See the retained software-metadata checksum exception before making provenance claims.

Python dependencies: Python 3.12, pandas and NumPy; acquisition uses the standard library. Run scripts in this sequence to reproduce the online search/acquisition and offline count reconstruction:

```sh
python scripts/acquire_search.py
python scripts/search_literature_supplement.py
python scripts/acquire_rr10.py
python scripts/freeze_rr10_counts.py
python scripts/summarize_search.py
python scripts/verify_and_manifest.py
```

Online results may evolve. The frozen `sources` responses and complete accession/literature decisions preserve this run. The outcome-independent plan precedes acquisition and should not be rewritten. The matrix reconstruction additionally compares individual source identifiers with `../../inputs/animal_metadata.csv`; it does not inspect RR23 outcomes. `tables/RR10_matrix_checks.csv` gives every count-file check; `tables/RR10_freeze_summary.json` records the final matrix and metadata hashes. `tables/FILE_MANIFEST.csv` hashes the extension inputs, sources, scripts and outputs present at handoff.
