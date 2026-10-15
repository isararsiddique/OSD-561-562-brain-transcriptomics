# Spaceflight brain transcriptomics review package

This is the authoritative 15 October 2026 deep-analysis revision for **npj Microgravity**, prepared as an Article for author review. It contains the full manuscript in Word, PDF and Markdown, supplementary information, ten main and five supplementary panel figures, an updated optional graphical overview, exact processed NASA inputs, complete result tables and executable code. Historical drafts are not current scientific outputs.

The completion audit refitted supervised ground-panel weights and training-only preprocessing under all 400 age-restricted label assignments in each core region. Sensor31 had AUC 1, exact p=0.005, BH q=0.01125 and BY q=0.03932 across all 18 inherited panel/comparator-by-region tests. Other panels also discriminate the cohorts; this is not Sensor31-specific evidence. All label nulls and failed tests are retained. The age-saturated Circadian estimates remain negative, but no saturated estimand survives its separate 64-test BH family (hippocampal equal-age q=0.113). This strengthens conditional collection association and qualifies the biological interpretation.

Publisher guidance was checked directly on 15 October 2026. The title has 12 words and the unstructured abstract 127 words. All Methods, mathematical derivations and reproducibility instructions are in the main manuscript; funding is in Acknowledgements; tables follow References. Artwork uses RGB 600 dpi PNGs, editable Arial SVGs, lowercase bold panels, at least 8 pt text and 1 pt strokes at the intended review-page size. Simulation figures also show the full retained repetition distributions. Exact guidance and literature comparisons are in `extensions/journal_audit_20261015/`.

The real-data contribution is the dissociation between reproducible ground collection scores and transfer of genome-wide correction loadings. Ground-trained sensor scores achieved AUC 1 in every core holdout, while original outside-sensor prediction worsened by 9.1% in cerebellum and 3.5% in RR-23 on average. The 362 splits overlap and involve 30 ground source animals. A qualified hippocampal Circadian flight-by-collection interaction is −0.405 log2 (wild p=0.000838, q=0.026816 in its separately declared 32-test exploratory family). Its negative direction agrees with count models and omission checks; metadata exclusions lose corrected support and collapse age-cell representation. Primary count-model gene interactions identify 109 cerebellar and seven hippocampal candidates, of which 90 and two retain global significance after exclusions. These are cohort-dependent candidates, not independently replicated causal microgravity genes.

Seven deposits contain 129 libraries, 126 collapsed animal–region observations and **83 animals**, with 24 core and ten RR-10 paired regional animals. They do not constitute seven independent missions. Eleven negative-binomial fits include seven regional robustness fits and four direct interaction fits. All 232 primary/sensitivity panel tests, 104 paired panel tests, 32 additional interaction tests, 54,989 cross-cohort gene screen entries, 72,023 evaluable interaction gene tests, and every simulation output are retained. Separate exploratory families and unavailable values are explicit.

The score operator is an oblique projection with real-data gain 8.13–15.31. A closest-loading constraint caps gain at two and retains sensor-score removal. It can leave nuisance behind and does not establish biological recovery. The known-truth tests contain 2,600 independent generated datasets; they are mechanism simulations, not additional real animals. No cross-cohort gene or paired panel survives its respective correction.

## Files to review

- `MANUSCRIPT_COMPLETE.pdf` and `.docx`: complete Article, with figure legends on the same page as each figure.
- `SUPPLEMENTARY_INFORMATION.pdf` and `.docx`: supplementary results, electronic result index and five supplementary figures. Derivations, configurations and reporting details are in the main Methods.
- `publication_figures/`: final numbered 600 dpi PNGs and editable-text SVGs; `GRAPHICAL_ABSTRACT` is optional artwork. All plots use actual saved results.
- `results/validation_report.json`: numerical and rendered-document checks. `validation_fault_injection_report.json` records rejection of deliberately altered numerical and text content.
- `extensions/literature_novelty/NOVELTY_ASSESSMENT.md`: verified prior-work comparison and claim boundaries.
- `JOURNAL_TARGET.md`: journal fit, verified SJR2025 Q1 categories and publisher-guidance links. No JCR quartile or acceptance is implied.
- `authors.json`: confirmed author identities, CRediT roles, funding and interests.
- `supplementary_data/INDEX.csv`: numbered electronic Supplementary Data, complete-table source paths, dimensions and SHA-256 hashes. Original analysis paths remain available for reproduction.
- `results/completion_environment.json`: the separate final artwork/document/verification environment; the original scientific-run versions remain in `results/run_summary.json` and audit records.
- `PACKAGE_FILE_MANIFEST.json`: SHA-256 hashes for the sealed portable package; the validator checks these when present.

The Word body and tables are editable. Display equations are rendered images with matching LaTeX in their accessibility descriptions and `figures/equations/manifest.json`; SVG/PNG equation sources are retained. Figures use Arial, lowercase panel labels and vector sources. Initial journal submission permits a combined Word/PDF; no mandatory graphical abstract was found in the checked guidance. No journal branding or acceptance information is fabricated.

## Reproduction

The scientific analyses used Python 3.12.13 with the versions recorded in `results/run_summary.json` and pinned in `requirements.txt`. Final artwork, document assembly and verification used a separate Python 3.12.14 environment recorded in `results/completion_environment.json`; existing count models and simulations were reused. Create a dedicated environment and install `requirements.txt` for scientific reconstruction. Plot reproduction requires Arial or an intentional documented font substitution; supplied artwork already uses Arial. Limit BLAS threads with:

```sh
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
python -m unittest discover -s tests -v
python scripts/validate_revision.py
```

The validator checks supplied outputs without refitting. Full scientific reconstruction is substantially longer:

```sh
python scripts/run_revision.py
python scripts/run_count_model_validation.py --force
python scripts/run_reproducibility_screen.py
python scripts/run_validation_extension.py
python extensions/paired_biology/run_paired_biology.py
python extensions/mathematics_review/independent_operator_tests.py
python extensions/mathematics_review/fixed_sensor_noise_inflation.py
python extensions/mathematics_review/high_gain_benchmark.py
python extensions/ground_holdout/run_ground_holdout.py
python extensions/ground_holdout/verify_ground_holdout.py
python extensions/real_data_structure/run_structure_qc.py
python extensions/flight_collection_interaction/run_flight_collection_interaction.py
python extensions/flight_collection_interaction/describe_biological_context.py
python extensions/flight_collection_interaction/validate_flight_collection_interaction.py
python extensions/independent_audit_20261015/run_independent_audit.py
python extensions/independent_audit_20261015/verify_exact_ground_null.py
python extensions/independent_audit_20261015/build_audit_figure.py
python scripts/prepare_manuscript_sources.py
python scripts/build_draft.py --figures-only
python scripts/build_extended_figures.py
python scripts/refine_panel_layouts.py
python scripts/build_graphical_abstract.py
python scripts/prepare_final_content.py
python scripts/finalize_for_journal.py
node scripts/refine_journal_artwork.cjs
python scripts/build_final_documents.py
python scripts/build_supplementary_data.py
```

Use `--force` for actual count-model refits; without it cached fits may be reused. The legacy `build_draft.py` must use **`--figures-only`**; its historical text/PDF builder is superseded by `prepare_final_content.py`, `finalize_for_journal.py` and `build_final_documents.py`. Render the resulting Word files using LibreOffice (provided PDFs were rendered from these Word files). The supplied equation images need no regeneration for normal reproduction. To regenerate equations or raster artwork, install Node, `mathjax-full` 3.2.2 and `sharp`, and set both `P1_MATH_MODULES` and `P1_RUNTIME_MODULES` to the absolute installed `node_modules` directory. The default Node module locations are local fallbacks, not portable requirements. Run `node scripts/render_equations.cjs` when equation sources change. The completion finalizer must run after `prepare_final_content.py`, which intentionally rebuilds the source blocks before journal formatting. Rerunning analyses invalidates a sealed package's file hashes; regenerate the manifest only after reviewing the new outputs.

Inputs are processed public NASA counts and metadata, with source URLs and hashes. This work does not repeat raw-read alignment, create new animal experiments, or resolve deposited location conflicts. All inferential assumptions and failed transfer configurations are retained.

## Release status

This is a complete computational **review package**, not a submitted article. All five authors must review and approve the final text and contributions. A permanent public data/code revision archive and its real identifier must be created before publication. Journal portal reporting-summary/checklist fields must be completed using factual main Table 4. The overlapping Paper 2 must be disclosed and its actual submission status checked. Originality, simultaneous-submission and author-approval declarations must not be asserted without confirmation. `--submission-ready` deliberately fails until these release items are handled. No external submission, upload or payment has been made.
