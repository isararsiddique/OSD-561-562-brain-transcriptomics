# Independent-launch validation extension

This dated retrospective extension tests the complete RR23 all26-defined local-junction profile in independently launched RRRM2/OSD-561 cerebellum. It does not select RR23 significance-positive events. The primary null result and measurement limits are in `INDEPENDENT_VALIDATION_REPORT.md`; numerical values are in `tables/summary.json`.

The extension is portable. Its output paths are script-relative. The only external input is a beta analysis/package root passed with `--base`, containing the frozen RR23 files below. All required RRRM2 counts, animal metadata, official file inventory, publisher MD5 list, publisher URLs and checksum/aggregation verification are frozen in `inputs/`. Raw STAR junction files are not required to rerun the analysis.

## Reproduction

Python 3.11+ with numpy, pandas, scipy and matplotlib is required. Run from this extension directory. If the extension is installed as `extensions/independent_validation` inside a beta package, these commands work directly:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python scripts/analyze_validation.py --base ../..
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python scripts/verify_independently.py --base ../..
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python scripts/build_report_and_figure.py
```

For another layout, replace `../..` with the beta analysis/package root. No path to the source research repository is needed. The analysis reads these beta files:

- `inputs/primary_all26_event_definitions.csv.gz`
- `inputs/all_animal_junction_counts.csv.gz`
- `inputs/animal_metadata.csv`
- `results/all26_animal_fractions.csv.gz`
- `results/measurement_gates.csv.gz`

The input manifest uses `scope=beta` or `scope=extension`, with relative paths and before/after SHA256 hashes. The frozen retrospective plan is preserved byte-for-byte. It records the original research paths as provenance, not as execution dependencies.

`scripts/freeze_source_inputs.py` was a one-time provenance step. It checked all 26 original OSD-561 junction files against NASA's publisher MD5 hashes and independently rebuilt every retained count in all 19 included animal columns (668,237 junction rows) before freezing the matrix. It is not a required reproduction step; its results and original file URLs are preserved in `inputs/frozen_input_provenance.json` and `inputs/RRRM2_publisher_provenance.csv`.

## Main outputs

- `figures/Figure_Independent_Validation.png`, `.pdf`, `.svg` and `Figure_Independent_Validation_caption.txt`: four-panel figure for the manuscript. The PDF/SVG preserve vector lines, text and points.
- `tables/complete_feature_mapping_and_gates.csv.gz`: all 5,724 definitions and measurement decisions.
- `tables/complete_common_panel_effects.csv.gz`: all 1,512 common events, effects, fixed learned directions and full-family Holm values.
- `tables/RRRM2_heldout_animal_scores.csv`, `RRRM2_exact_score_null.csv.gz` and `animal_bootstrap_draws.csv.gz`: all biological animal scores, 11,200 exact assignments and 2,000 bootstrap draws.
- `tables/independent_checks.json`: separate implementation verification of complete fractions/counts/dominance, projection, all score assignments, deterministic event checks and all Holm values.

The long descriptive filenames `independent_launch_validation_four_panel.png` and `.svg` contain the same figure as the manuscript aliases. All source files are preserved. Author: Isarar Siddique; affiliation and email are blank.
