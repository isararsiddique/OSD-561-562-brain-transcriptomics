# Strengthened spaceflight brain transcriptomics review package

This is the current 15 October 2026 computational revision for **npj Microgravity**. The Article is **Gain constrained regression improves transcriptomic prediction across spaceflight and Earth brain datasets**. It contains 44 manuscript pages, 11 supplementary pages, 13 main figures, seven supplementary figures, an updated optional graphical overview and 50 complete numbered Supplementary Data tables. The earlier sealed review package is preserved separately.

The new constructive result is a gain-two convex multi-sensor regression method. With training-only preprocessing and nested tuning, mean outside-sensor prediction error falls by **1.2–5.2% in the three original settings and 1.4–7.2% in four external datasets**. All 16,084 fitted optimization problems meet a relative gap of 1e-7, and 96 independent full-coefficient CVXPY/Clarabel checks pass. A matched-grid comparator adds 852 rows across 426 holdouts and 28 fixed-alpha objective checks. The constrained estimator has lower equal-unit/split mean error than matched clipping in all seven settings, but OSD-611 pooled error and paired median favor clipping. Unbounded ridge can predict better at larger amplification. There is no universal-superiority claim.

The original resource retains seven deposits, 129 libraries and 83 animals. The external resource has 66 acquired columns and 64 primary expression units: **46 individual animals and 18 two-animal pools**. RR-3 adds a distinct flight mission; the other datasets are separate Earth irradiation/unloading experiments. External testing refits the algorithm within each dataset, rather than transporting frozen fitted weights. Pool donors, RNA integrity, lane structure, anatomy and exposure differences remain explicit limitations.

The additional real-count-calibrated benchmark contains 1,400 independent negative-binomial datasets and 14 pipelines, including official ComBat, published-algebra RUVg, known-batch regression and distinct PCA adjustment. All 19,600 method runs completed. Outcomes show tradeoffs and failures rather than universal biological recovery. These simulations are not additional animals. Earlier 2,600 mechanism/stress simulations remain fully available.

Biological analysis retains 1,488 effect/screen rows and 1,047 finite tests in one new BH/BY family. The hippocampal Circadian interaction remains negative across ages and changes only 1.81% after the fixed erythroid-marker adjustment, but adjusted statistical support is incomplete and no cross-context biological replication is established. All unfavorable sensitivities and unavailable probabilities are retained. There are 1,337 independent biological checks; the additional whole-extension validator passes 1,187 identity, provenance, summary and family checks. These are internal computational validations, not causal certification.

## Review files

- MANUSCRIPT_COMPLETE.pdf and .docx: complete Article with all Methods, 12 displayed equations, 38 references, 6 tables and 13 figure pages. Markdown and equation LaTeX sources are included.
- SUPPLEMENTARY_INFORMATION.pdf and .docx: supplementary results, original table index, extension output descriptions and 7 figure pages.
- publication_figures/: 600 dpi RGB PNGs and editable Arial SVGs. The graphical overview is optional; no mandatory graphical abstract was found in publisher guidance.
- supplementary_data/INDEX.csv: 50 complete-table exports, dimensions, source paths and SHA-256 hashes. No discovery-only filtering is used.
- extensions/strengthening_20261015/method_benchmark/METHOD_THEORY_AND_RESULTS.md: proofs, nested comparison, external results and matched-grid exceptions.
- extensions/strengthening_20261015/external_data/DATASET_PROTOCOL.md: acquisition, experimental units and design limits.
- extensions/strengthening_20261015/biological_validation/BIOLOGICAL_VALIDATION_REPORT.md: complete biological sensitivity interpretation.
- extensions/strengthening_20261015/literature/STRENGTHENING_NOVELTY_ASSESSMENT.md: prior work, contribution scope and journal evidence.
- results/validation_report.json: 1858 current numerical, provenance and rendered-content checks. Independent reports and two rejected fault probes accompany it.
- results/visual_review.json: every rendered page reviewed against the final PDF hashes.
- PACKAGE_FILE_MANIFEST.json and PACKAGE_VERIFICATION.json: portable-package seal and verification, added during packaging.

## Reproduction

Preserve the supplied outputs before refitting. Use a dedicated Python 3.12 environment; requirements_strengthening.txt and runtime_strengthening.txt describe this extension, while requirements.txt and runtime_versions.txt describe historical original analyses. Module versions can differ from stale bundled distribution metadata; the observed extension module versions are recorded. Scientific plotting requires Arial or a documented replacement. BLAS thread limits used were OPENBLAS_NUM_THREADS=1, OMP_NUM_THREADS=1 and VECLIB_MAXIMUM_THREADS=1.

To verify supplied results, run python scripts/validate_revision.py and python -m unittest discover -s tests -v. The validator recomputes full aggregate families and document consistency without rerunning all models. Source CSV/GZ matrices, metadata, protocols and complete tables accompany the package; historical acquisition URLs and original absolute paths are provenance, while analysis manifests resolve portable relative paths.

For extension reconstruction, run the scripts in extensions/strengthening_20261015 in this order: method_benchmark/run_method_benchmark.py; method_benchmark/run_strengthened_methods.py; method_benchmark/run_matched_grid_comparator.py; correction_benchmark/run_count_benchmark.py --repetitions 100; biological_validation/prepare_external_metadata.py; biological_validation/run_biological_validation.py. External inputs are already included; reacquisition is optional and its script documents the original source locations. Independent verifiers accompany each analysis. The biological design document has an exact at-analysis snapshot plus a dated correction of the preliminary OSD-611 group expectation; the fitted data always used the verified primary 3/5/2 design.

For document/artwork reconstruction, run biological_validation/build_biological_report.py, scripts/build_strengthening_figures.py, scripts/integrate_strengthening.py, node scripts/render_equations.cjs, node scripts/refine_journal_artwork.cjs, scripts/build_final_documents.py and scripts/build_supplementary_data.py. The integration script reads supplied pre-strengthening text snapshots and the current numerical tables. MathJax 3.2.2 and sharp are required only to regenerate equation/vector raster artwork; point P1_MATH_MODULES and P1_RUNTIME_MODULES to their installed node_modules paths. Render Word to PDF with LibreOffice, inspect every page, then validate. The original analysis reproduction order remains in the main Methods and its source scripts. Recalculation invalidates a sealed package's hashes; create a new reviewed seal afterward.

## Publication status

Publisher guidance was checked on 15 October 2026. The 12-word title and 127-word unstructured abstract meet the checked limits. Methods and derivations are in the main file; funding is in Acknowledgements; figure legends are on their figure pages. All final artwork has at least 8 pt text and 1 pt strokes at the intended review-page size. JOURNAL_TARGET.md distinguishes the verified SJR 2025 quartiles from unverified JCR claims.

This is a complete computational author-review package. It has not been submitted. Final all-author approval, a permanent public revision archive/identifier, journal portal reporting forms and disclosure of the related Paper 2 remain necessary. Authorship, funding and competing-interest statements follow the confirmed author record. No acceptance, unprecedented general algorithm or independently replicated biological mechanism is asserted.

The complete historical original-analysis command sequence is preserved in results/README_pre_strengthening.md. Its earlier title, figure counts and conclusions describe the pre-extension revision; use the current files and extension workflow above for this release. The biological figure builder is extensions/strengthening_20261015/biological_validation/build_biological_report.py when invoked from the package root.
