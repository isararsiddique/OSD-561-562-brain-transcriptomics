# Frozen-family multiway robustness extension

Read `RETROSPECTIVE_PLAN_2026-10-15.md` and `SCIENTIFIC_REPORT.md`. This extension uses the NEW npj package inputs and its sibling `extensions/multiway_splicing`; it does not depend on the predecessor research directory. `results/frozen_input_links.csv` lists the package-relative required inputs and hashes.

From this directory:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
MWR_PYTHON=/path/to/python
MWR_PACKAGE=../..
"$MWR_PYTHON" scripts/analyse_robustness.py --base "$MWR_PACKAGE"
"$MWR_PYTHON" scripts/verify_robustness.py --base "$MWR_PACKAGE"
"$MWR_PYTHON" scripts/build_report_figure.py
```

NumPy, pandas, SciPy and Matplotlib are required. Arial or Helvetica must be available for journal figure generation. Versions are recorded in `results/software_versions.json`. All outputs resolve to this script-relative extension, and no network is needed. The analysis reproduces ten complete 805,562×26 binomial-thinned matrices and all inferential/bootstrap results; verification is intentionally comprehensive and may take several minutes. All ten matrices are retained in `thinned_counts` as NPZ files in the original frozen matrix row order, with source/hash/seed metadata in `results/thinned_count_manifest.csv`. They are read perturbations of the same animals, not independent biological datasets.

Core outputs:

- `results/local_dominance95_complete_fixed_family.csv.gz`: all 2,769 sensitivity hypotheses, invalid p=1.
- `results/local_thinning_complete_fixed_families.csv.gz`: all 27,690 local thinning tests, no dropped hypotheses.
- `results/global_thinning_all_replicates.csv`: all 30 global results.
- `results/original_hit_stability_summary.csv` and complete tracking: every original positive test, including failures.
- `results/all_cluster_animal_bootstrap_intervals.csv.gz`, signed alternative intervals, and `animal_bootstrap_replicates.npz`: complete pointwise uncertainty and resampling weights.
- `results/dispersion_and_studentization_diagnostics.csv` and all-arm dispersions: carefully limited reference diagnostics.
- `results/independent_verification.json`/`.csv`: all independent checks.
- `figures/Figure_Multiway_Robustness.pdf`/`.svg`/`.png` and caption: fixed 183 mm width, 8 pt font, ≥1 pt lines, 600 dpi PNG.
- `MANIFEST.csv`: SHA256 of every deliverable except the manifest itself.

Do not promote a bootstrap lower bound to simultaneous local significance, interpret nonsignificant dispersion as equal variances, select favorable thinnings, or treat stable aggregate association as validation of individual depth-sensitive local discoveries. Flight/date confounding remains unchanged.
