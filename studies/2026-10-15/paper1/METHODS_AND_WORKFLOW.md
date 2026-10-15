# Paper 1: methods, command order and interpretation

This guide summarizes the final manuscript Methods and connects them to the preserved source scripts. The work is a retrospective secondary analysis of processed NASA data. Plans, previous results and diagnostic amendments are dated research provenance rather than prospective preregistration.

## 1. Establish source identity and experimental units

Use OSD-525, 561, 562, 563, 564, 612 and 613. Resolve animals from deposited Source Name and cross-check RFID. Sum technical libraries within animal and deposit before filtering. A mouse contributes one expression profile per region; multiple tissues from one animal remain linked during inference.

Cross-tabulate euthanasia method, sample/carcass preservation, collection and dissection dates, age, exposure/recovery and habitat against each contrast. A protocol contrast combines the deposited protocol differences and does not identify the separate causal effect of one anaesthetic or preservation step.

## 2. Normalize separately within deposits

Remove ERCC rows. Retain Ensembl genes with mean raw count at least 10. Calculate median-ratio size factors using genes positive in every sample, then use `log2(count / size_factor + 1)`. Normalize deposits separately. Pooled regional analyses use the gene intersection and tissue/age nuisance terms.

The source loaders and checksum gates are in `scripts/run_revision.py`; mathematical primitives are in `scripts/statistics_core.py`. STAR is primary and matched RSEM inputs support the estimator sensitivity.

## 3. Estimate panel effects and exact probabilities

Eight fixed panels cover rapid/delayed primary response, secondary response, heat shock, glia/complement, endothelial/BBB, circadian genes and the 31-gene collection sensor. Frisch–Waugh–Lovell residualization removes nuisance terms; a panel effect averages its gene coefficients.

Enumerate whole-animal allocations preserving age-block group sizes and keeping regions from one animal together. Include the observed allocation and ties. The primary BH family contains all 232 panel tests across 17 contrasts and metadata sensitivities. Percentile intervals use 1,000 whole-animal resamples within label and age. Gene-level sandwich uncertainties are exploratory.

## 4. Decompose protocol and flight effects

Within RRRM-2 estimate ground ISS-terminal minus ground live-animal-return, the corresponding flight protocol contrast, and flight-minus-ground within each protocol. Before calculating new quantities, the protocol script reproduces 80 existing panel effects and exact probabilities.

The direct dominance contrast compares ground ISS-terminal with flight live-animal-return relative to their shared ground live-animal-return reference. Its complete family contains 8 panels × 3 analyses = 24 tests, with BH correction. Cerebellar, hippocampal and pooled exact spaces contain 200, 400 and 700 assignments. Summarize protocol-minus-flight differences with 2,000 stratified whole-animal bootstraps, seed 15102026. All family members, including negative findings, are retained.

## 5. Check count models, interactions and regional pairing

Seven PyDESeq2 0.5.4 models use collapsed integer counts, parametric dispersions, Cook refitting and no coefficient shrinkage; 15 regional contrasts are compared with linear panel estimates. Flight-by-collection interactions use age-adjusted full designs, 4,000 stratified bootstraps and 499,999 restricted-model Rademacher wild draws with HC3 studentization; the declared panel family has 32 tests.

The cross-cohort screen uses the maximum component Wald probability as an intersection–union test. Discordant or unavailable effects receive p=1; all 54,989 tests enter one BH family with BY sensitivity. Paired hippocampus-minus-cerebellum analyses retain within-animal linkage and use a separate 104-test family. Lack of corrected replication is not proof of no biological effect.

## 6. Derive and test correction geometry

With unit sensor weights `w`, residual expression `Y0` and score `s = w.T @ Y0`, estimate loading `gamma = Y0 @ s.T / (s @ s.T)`. Score removal is `A = I - gamma @ w.T`. It is idempotent and removes the sensor score, but its spectral norm equals `norm(gamma)` for the stated multidimensional setting.

Writing `gamma = w + h`, with `h` orthogonal to `w`, gives gain `sqrt(1 + norm(h)**2)`. A bounded loading shrinks the perpendicular component to cap gain at `L=2`. Tests verify idempotence, the score-null constraint and spectral gain. Removing score energy does not quantify the fraction of false biological discoveries removed.

## 7. Evaluate ground-only holdouts without leakage

Each RRRM-2 region uses all 81 splits holding out one animal from each of four age-by-cohort ground cells: 8 training and 4 test animals. RR23 uses 200 fixed ground-control splits. Refit gene filtering, normalization references, nuisance coefficients and sensor weights using training animals only.

Report age-stratified discrimination and outside-sensor prediction relative to the nuisance baseline. Rebuild the entire pipeline under every age-restricted label allocation; overlapping splits are not independent replicates. The exact discrimination family has 18 tests with BH/BY adjustment.

The identical-loading noise-ceiling analysis asks whether low loading agreement could arise even under one shared true loading. It estimates contrast/noise energies, checks the approximation with 400 Gaussian simulations and explicitly ignores gene correlations. This is a diagnostic, not a transport hypothesis test.

## 8. Fit gain-constrained regression

Minimize the training squared error plus sensor-scale ridge penalty subject to `norm(C, 2) <= sqrt(L**2 - 1)`, with `L=2`. The induced correction operator has gain at most two. Accelerated projected-gradient steps use singular-value clipping, monotone restart and a Frank–Wolfe gap stopping threshold of `1e-7` relative to the stated objective scale.

All 16,084 non-baseline fits reached the recorded criterion. Twelve full-coefficient CVXPY/Clarabel solutions provide 96 independent objective, feasibility, gain and identity checks. Nested leave-one-training-unit-out selection considers the baseline and alpha values 10, 1, 0.1 and 0.01. Matched-grid ridge/clipped-ridge comparisons use identical folds, preprocessing and candidate grids. Post-fit clipping and constrained loss minimization answer different optimization questions.

## 9. Evaluate external datasets and simulations

External OSD-202, 352, 374 and 611 evaluations refit the algorithm within each dataset. They do not transport frozen coefficients across experiments. Retain the recorded individual-animal versus two-donor-pool units; use the OSD-611 primary SampleTable design and complete-pool sensitivity.

The count-calibrated benchmark samples 600 genes, estimates within-cell negative-binomial dispersion, evaluates seven scenarios at training sizes 8 and 24 with 100 datasets per cell, and compares 14 pipelines. BH operates over 600 genes within each method run. The 1,400 simulated datasets and 19,600 method runs are computational experiments, not additional animals. Report both favorable and unfavorable comparisons.

## Command order

Restore the complete ZIP first. From the repository root:

```sh
cd .local-reproduction/paper1
python3.12 -m venv .venv
.venv/bin/python -m pip install -r Paper_1_Complete_Submission_Package/Paper1_Strengthened_Review_20261015/requirements_strengthening.txt
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
```

The archive also records historical core requirements and observed runtime versions. Install the packages needed for the selected stage in its dedicated environment; rendering-only Node/Word dependencies are not required for numerical analyses.

| Stage | Command from the restored archive root |
|---|---|
| Core primary analysis | `.venv/bin/python Paper_1_Complete_Submission_Package/Paper1_Strengthened_Review_20261015/scripts/run_revision.py` |
| Protocol decomposition/dominance | `.venv/bin/python Paper_1_npj_LaTeX_Build_20261015/research/D_protocol/D_analysis.py` |
| Published signature (hash-gated) | `.venv/bin/python Paper_1_npj_LaTeX_Build_20261015/research/F_signatures/F_analysis.py` |
| Labelled post hoc gene classes | `.venv/bin/python Paper_1_npj_LaTeX_Build_20261015/research/F_signatures/posthoc_tyssowski.py` |
| Independent numeric recomputation | `.venv/bin/python Paper_1_npj_LaTeX_Build_20261015/research/A_audit/A_recompute.py a b c d e f g h i j k l m` |
| Noise ceiling | `.venv/bin/python Paper_1_npj_LaTeX_Build_20261015/research/B_novelty/analysis/noise_ceiling.py` |
| Independent reimplementation | `.venv/bin/python Paper_1_npj_LaTeX_Build_20261015/research/G_independent/G_code/g_part13.py`, then `g_part2.py` |

For strengthened method reconstruction, change to the primary package root and run the scripts in `extensions/strengthening_20261015/` in order:

1. `method_benchmark/run_method_benchmark.py`
2. `method_benchmark/run_strengthened_methods.py`
3. `method_benchmark/run_matched_grid_comparator.py`
4. `correction_benchmark/run_count_benchmark.py --repetitions 100`
5. `biological_validation/prepare_external_metadata.py`
6. `biological_validation/run_biological_validation.py`

The corresponding `verify_*.py` scripts independently check the generated results. Full statistical regression tests are:

```sh
cd Paper_1_Complete_Submission_Package/Paper1_Strengthened_Review_20261015
../../.venv/bin/python -m unittest discover -s tests -v
```

`scripts/validate_revision.py` additionally checks supplied aggregate families and historical document consistency. Publication-document checks refer to superseded documents omitted from the analysis-only ZIP; distinguish those missing-document checks from computational checks. Do not report a full document-validator pass for this stripped release.

## Results and limits

Use `research/D_protocol/D_tables/` for protocol quantities, the primary package's `results/`, `extensions/` and `supplementary_data/` for panel/correction/benchmark outputs, and `research/A_audit/` plus `research/G_independent/` for independent recomputation. Larger tables remain in the complete Zenodo ZIP, with paths and hashes listed in the release manifest.

The main conclusion supports matching collection protocols by design and reporting correction gain. Observational confounding, small animal groups, overlapping holdouts, pooled donors and cross-experiment differences prevent claims of causal artefact removal, universal prediction superiority or independently replicated molecular mechanisms.
