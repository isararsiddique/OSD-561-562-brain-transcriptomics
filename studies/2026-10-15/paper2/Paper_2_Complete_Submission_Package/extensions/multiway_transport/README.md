# Frozen multiway transport extension

This retrospective analysis tests all 923 frozen RR23 multiway junction clusters in independently launched RRRM2/OSD-561 cerebellum. It learns unit Euclidean native-composition directions from RR23 only and projects complete measurable compositions into held-out animals. No RR23 significance subset is selected.

The extension is self-contained. Frozen count matrices, metadata, complete alternative memberships/features, annotation lookup, prior measurement context and publisher provenance are under `inputs/`; no source-repository path or other extension is needed for analysis. The source-freezing script is a one-time provenance step, not a reproduction dependency.

## Reproduce

Python 3.11+ with numpy, pandas, scipy and matplotlib is required. Arial or Helvetica must be installed for the publication figure. Run these commands from this extension directory:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python scripts/analyse.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python scripts/verify_independently.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python scripts/build_outputs.py
```

Outputs are script-relative, so the scripts may also be called by absolute path from another working directory. Optional self-contained copied-folder verification is `python scripts/check_offline_reproduction.py` after all three commands. It reruns all three commands in a fresh copied folder, checks an identical complete summary, scientific tables, bootstrap arrays, report and caption, and passing independent checks. Its temporary copy is removed afterward. Tested dependency versions are recorded in `requirements-tested.txt` and `tables/software_versions.json`.

## Main files

- `RETROSPECTIVE_PLAN_2026-10-15.md`: outcome-preceding dated retrospective plan and fixed normalization.
- `SCIENTIFIC_REPORT.md`: primary null findings, measurement limits, complete-family correction and finite-label-space Holm bound.
- `figures/Figure_Multiway_Transport.pdf`, `.svg`, `.png` and `_caption.txt`: 183 mm four-panel figure, ≥8 pt Arial/Helvetica, ≥1 pt lines, 600 dpi PNG and separate caption. Group encodings use FLT blue `#0072B2` and GC orange `#D55E00`; cluster/uncertainty encodings are identified separately.
- `tables/all_923_cluster_directional_results.csv`: all 923 family members, gates, learned norms, held-out directional estimates and BH/BY/Holm adjustments; failed measurement p values remain 1.
- `tables/all_3376_alternative_group_composition_summaries.csv.gz`: every frozen alternative, all RR23 group means and held-out route×age-weighted means for common measurable clusters; no component-wise hypothesis tests.
- `tables/animal_bootstrap_cluster_effects.npz` and `animal_bootstrap_global_draws.csv.gz`: all 2,000 conditional fixed-direction and secondary two-cohort retrained-direction draws.
- `tables/exact_global_null_all_11200.csv.gz`: every exact stratified global allocation.
- `tables/independent_numerical_checks.json`: separate full raw-count reconstruction, direction/projection, allocation, correction and bootstrap verification.
- `tables/offline_reproduction_check.json`: PASS for all three commands from copied inputs/scripts; 28 scientific CSV tables, all saved bootstrap arrays, summary, report and caption match exactly.
- `tables/figure_visual_QA.json`: hash-bound PASS after inspecting all four panels of the final PDF render; `artifact_sha256_manifest.csv` inventories the completed extension.

Global and local families are separate retrospective analyses. Local Holm discovery at 0.05 is structurally unattainable for 923 tests and 11,200 assignments; zero Holm calls must not be treated as biological evidence of absence. The primary one-score exact test and BH/BY map remain informative under their stated conditions. RR23 flight and collection date are aliased. No causal flight, injury or conserved universal signature is claimed.

Author: Isarar Siddique. Affiliation and email are blank.
