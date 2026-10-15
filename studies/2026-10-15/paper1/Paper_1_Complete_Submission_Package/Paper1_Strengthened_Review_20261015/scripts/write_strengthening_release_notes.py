#!/usr/bin/env python3
"""Generate current review-package notes only after numerical and visual checks."""
from pathlib import Path
import json
R=Path(__file__).resolve().parents[1]
v=json.loads((R/'results/validation_report.json').read_text());qa=json.loads((R/'results/visual_review.json').read_text())
assert v['computational_validation_passed'] and qa['all_pages_reviewed']
history=R/'results/README_pre_strengthening.md'
if not history.exists():history.write_text((R/'README.md').read_text())
main=v['pdf_pages']['MANUSCRIPT_COMPLETE'];supp=v['pdf_pages']['SUPPLEMENTARY_INFORMATION'];checks=v['computational_checks']
(R/'README.md').write_text(f'''# Strengthened spaceflight brain transcriptomics review package

This is the current 15 October 2026 computational revision for **npj Microgravity**. The Article is **Gain constrained regression improves transcriptomic prediction across spaceflight and Earth brain datasets**. It contains {main} manuscript pages, {supp} supplementary pages, 13 main figures, seven supplementary figures, an updated optional graphical overview and 50 complete numbered Supplementary Data tables. The earlier sealed review package is preserved separately.

The new constructive result is a gain-two convex multi-sensor regression method. With training-only preprocessing and nested tuning, mean outside-sensor prediction error falls by **1.2–5.2% in the three original settings and 1.4–7.2% in four external datasets**. All 16,084 fitted optimization problems meet a relative gap of 1e-7, and 96 independent full-coefficient CVXPY/Clarabel checks pass. A matched-grid comparator adds 852 rows across 426 holdouts and 28 fixed-alpha objective checks. The constrained estimator has lower equal-unit/split mean error than matched clipping in all seven settings, but OSD-611 pooled error and paired median favor clipping. Unbounded ridge can predict better at larger amplification. There is no universal-superiority claim.

The original resource retains seven deposits, 129 libraries and83 animals. The external resource has66 acquired columns and64 primary expression units: **46 individual animals and18 two-animal pools**. RR-3 adds a distinct flight mission; the other datasets are separate Earth irradiation/unloading experiments. External testing refits the algorithm within each dataset, rather than transporting frozen fitted weights. Pool donors, RNA integrity, lane structure, anatomy and exposure differences remain explicit limitations.

The additional real-count-calibrated benchmark contains1,400 independent negative-binomial datasets and14 pipelines, including official ComBat, published-algebra RUVg, known-batch regression and distinct PCA adjustment. All19,600 method runs completed. Outcomes show tradeoffs and failures rather than universal biological recovery. These simulations are not additional animals. Earlier2,600 mechanism/stress simulations remain fully available.

Biological analysis retains1,488 effect/screen rows and1,047 finite tests in one new BH/BY family. The hippocampal Circadian interaction remains negative across ages and changes only1.81% after the fixed erythroid-marker adjustment, but adjusted statistical support is incomplete and no cross-context biological replication is established. All unfavorable sensitivities and unavailable probabilities are retained. There are1,337 independent biological checks; the additional whole-extension validator passes1,187 identity, provenance, summary and family checks. These are internal computational validations, not causal certification.

## Review files

- MANUSCRIPT_COMPLETE.pdf and .docx: complete Article with all Methods,12 displayed equations,38 references,6 tables and13 figure pages. Markdown and equation LaTeX sources are included.
- SUPPLEMENTARY_INFORMATION.pdf and .docx: supplementary results, original table index, extension output descriptions and 7 figure pages.
- publication_figures/:600 dpi RGB PNGs and editable Arial SVGs. The graphical overview is optional; no mandatory graphical abstract was found in publisher guidance.
- supplementary_data/INDEX.csv:50 complete-table exports, dimensions, source paths and SHA-256 hashes. No discovery-only filtering is used.
- extensions/strengthening_20261015/method_benchmark/METHOD_THEORY_AND_RESULTS.md: proofs, nested comparison, external results and matched-grid exceptions.
- extensions/strengthening_20261015/external_data/DATASET_PROTOCOL.md: acquisition, experimental units and design limits.
- extensions/strengthening_20261015/biological_validation/BIOLOGICAL_VALIDATION_REPORT.md: complete biological sensitivity interpretation.
- extensions/strengthening_20261015/literature/STRENGTHENING_NOVELTY_ASSESSMENT.md: prior work, contribution scope and journal evidence.
- results/validation_report.json: {checks} current numerical, provenance and rendered-content checks. Independent reports and two rejected fault probes accompany it.
- results/visual_review.json: every rendered page reviewed against the final PDF hashes.
- PACKAGE_FILE_MANIFEST.json and PACKAGE_VERIFICATION.json: portable-package seal and verification, added during packaging.

## Reproduction

Preserve the supplied outputs before refitting. Use a dedicated Python3.12 environment; requirements_strengthening.txt and runtime_strengthening.txt describe this extension, while requirements.txt and runtime_versions.txt describe historical original analyses. Module versions can differ from stale bundled distribution metadata; the observed extension module versions are recorded. Scientific plotting requires Arial or a documented replacement. BLAS thread limits used were OPENBLAS_NUM_THREADS=1, OMP_NUM_THREADS=1 and VECLIB_MAXIMUM_THREADS=1.

To verify supplied results, run python scripts/validate_revision.py and python -m unittest discover -s tests -v. The validator recomputes full aggregate families and document consistency without rerunning all models. Source CSV/GZ matrices, metadata, protocols and complete tables accompany the package; historical acquisition URLs and original absolute paths are provenance, while analysis manifests resolve portable relative paths.

For extension reconstruction, run the scripts in extensions/strengthening_20261015 in this order: method_benchmark/run_method_benchmark.py; method_benchmark/run_strengthened_methods.py; method_benchmark/run_matched_grid_comparator.py; correction_benchmark/run_count_benchmark.py --repetitions100; biological_validation/prepare_external_metadata.py; biological_validation/run_biological_validation.py. External inputs are already included; reacquisition is optional and its script documents the original source locations. Independent verifiers accompany each analysis. The biological design document has an exact at-analysis snapshot plus a dated correction of the preliminary OSD-611 group expectation; the fitted data always used the verified primary3/5/2 design.

For document/artwork reconstruction, run biological_validation/build_biological_report.py, scripts/build_strengthening_figures.py, scripts/integrate_strengthening.py, node scripts/render_equations.cjs, node scripts/refine_journal_artwork.cjs, scripts/build_final_documents.py and scripts/build_supplementary_data.py. The integration script reads supplied pre-strengthening text snapshots and the current numerical tables. MathJax3.2.2 and sharp are required only to regenerate equation/vector raster artwork; point P1_MATH_MODULES and P1_RUNTIME_MODULES to their installed node_modules paths. Render Word to PDF with LibreOffice, inspect every page, then validate. The original analysis reproduction order remains in the main Methods and its source scripts. Recalculation invalidates a sealed package's hashes; create a new reviewed seal afterward.

## Publication status

Publisher guidance was checked on15 October2026. The12-word title and127-word unstructured abstract meet the checked limits. Methods and derivations are in the main file; funding is in Acknowledgements; figure legends are on their figure pages. All final artwork has at least 8 pt text and 1 pt strokes at the intended review-page size. JOURNAL_TARGET.md distinguishes the verified SJR2025 quartiles from unverified JCR claims.

This is a complete computational author-review package. It has not been submitted. Final all-author approval, a permanent public revision archive/identifier, journal portal reporting forms and disclosure of the related Paper2 remain necessary. Authorship, funding and competing-interest statements follow the confirmed author record. No acceptance, unprecedented general algorithm or independently replicated biological mechanism is asserted.
'''.replace('and83','and 83').replace('has66','has 66').replace('and64','and 64').replace('and18','and 18').replace('contains1','contains 1').replace('and14','and 14').replace('All19','All 19').replace('Earlier2','Earlier 2').replace('retains1','retains 1').replace('and1','and 1').replace('only1','only 1').replace('are1','are 1').replace('passes1','passes 1').replace('Methods,12','Methods, 12').replace('equations,38','equations, 38').replace('references,6','references, 6').replace('and13','and 13').replace('and 7','and 7').replace('/:600','/: 600').replace(':50',': 50').replace('Python3','Python 3').replace('--repetitions100','--repetitions 100').replace('and100','and 100').replace('primary3','primary 3').replace('MathJax3','MathJax 3').replace('on15','on 15').replace('October2026','October 2026').replace('The12','The 12').replace('and127','and 127').replace('least8','least 8').replace('and 1 pt','and 1 pt').replace('SJR2025','SJR 2025').replace('Paper2','Paper 2'))
(R/'REVISION_STATUS.md').write_text(f'''# Author review status — 15 October 2026

The computational strengthening and journal-format revision are complete. The final Article has {main} pages, 13 main figures, 12 main equations and 38 references; Supplementary Information has {supp} pages and 7 figures. Fifty numbered electronic tables retain complete results. Every final rendered page was visually reviewed. Numerical/provenance/rendered-content validation passed {checks} checks, plus independent 96-check convex verification, 1,337-check biological verification and 11 statistical unit tests. The two deliberate numerical/text corruptions were rejected and restored before clean revalidation.

New evidence comprises gain-constrained convex regression, four external datasets with 64 primary animal/pool units,852 matched-grid comparator rows,28 fixed-alpha comparisons,1,400 new count simulations with 14 pipelines and1,488 new biological effect/screen rows. The complete manuscript reports improved mean prediction, stability bounds, comparator advantages, OSD-611 weighting exceptions and absence of corrected cross-context biological replication.

The strongest defensible novelty is the applied, certified correction framework and its external evaluation. Ridge, spectral constraints, singular-value clipping and proximal optimization are established. Q1 suitability remains an editorial judgement; no publication or acceptance claim has been manufactured.

Release items: all-author final approval; permanent public revision data/code archive; journal portal reporting documents; related Paper2 disclosure with its actual status. No journal submission, external upload, outreach or payment has been made.
''')
(R/'COVER_LETTER_DRAFT.md').write_text('''# Cover letter for author review

15 October 2026

Dear Editors of npj Microgravity,

Please consider our Article, “Gain constrained regression improves transcriptomic prediction across spaceflight and Earth brain datasets”. The corresponding author is Isarar Siddique, Convolity AI Lab, Lucknow, India; isarar@convolity.com; ORCID0009-0005-2305-9579.

The study addresses reliable analysis of collection-associated variation in spaceflight brain transcriptomics. We connect a multi-sensor correction operator to a convex prediction objective with a fixed spectral gain bound and verifiable optimality certificates. Training-only nested validation improves mean outside-sensor prediction by1.2–5.2% across the original settings and1.4–7.2% across64held-out animal/pool units from four additional datasets, including a separate flight mission. A matched-grid comparison isolates tuning choices, and all unfavorable individual-unit and alternative-weighting results are retained.

Seven original and four external public NASA deposits, explicit animal/pool provenance, complete endpoint families and reproducible code accompany the manuscript. Fourteen pipelines in1,400real-count-calibrated simulations show where established methods perform better and why stable prediction does not establish biological recovery. A hippocampal Circadian candidate is characterized through age, marker, omission and external-context sensitivity without asserting independent replication. The contribution combines constructive methodology with practical evidence and limits for comparative spaceflight omics; the established mathematical and unwanted-variation literature is acknowledged.

The authors report no external funding and no competing interests. Substantial AI assistance is disclosed in Methods, with scientific responsibility assigned to the authors.

Before this draft is sent, the corresponding author must confirm every author's approval, originality and exclusive-submission status, insert the permanent revision archive identifier, and disclose the related Paper2 with its verified title and current status. These declarations are not assumed in this draft.

Sincerely,

Isarar Siddique
Convolity AI Lab, Lucknow, India
isarar@convolity.com
'''.replace('ORCID0009','ORCID 0009').replace('by1.2','by 1.2').replace('and1.4','and 1.4').replace('across64held','across 64 held').replace('in1,400real','in 1,400 real').replace('Paper2','Paper 2'))
p=R/'JOURNAL_TARGET.md';text=p.read_text().replace('Prepare an Article around cohort confounding and the limits of ground-anchor transfer.','Prepare an Article around gain-constrained prediction, external algorithm evaluation and the distinction between stability and biological recovery.')
text+='\nThe strengthened revision adds four external datasets, matched-grid regression comparisons, certified convex optimization and established-method count simulations. Current novelty and design reports are in extensions/strengthening_20261015/literature/STRENGTHENING_NOVELTY_ASSESSMENT.md and external_data/DATASET_PROTOCOL.md. These improvements address the earlier within-resource validation and missing-comparator weaknesses; they do not establish a universally superior algorithm or independently replicated biological mechanism.\n';p.write_text(text)
for name in ['README.md','REVISION_STATUS.md','COVER_LETTER_DRAFT.md']:
 p=R/name;t=p.read_text()
 for old,new in [('units,852','units, 852'),('rows,28','rows, 28'),('comparisons,1,400','comparisons, 1,400'),('and1,488','and 1,488'),('Paper2','Paper 2')]:t=t.replace(old,new)
 if name=='README.md':t+='\nThe complete historical original-analysis command sequence is preserved in results/README_pre_strengthening.md. Its earlier title, figure counts and conclusions describe the pre-extension revision; use the current files and extension workflow above for this release. The biological figure builder is extensions/strengthening_20261015/biological_validation/build_biological_report.py when invoked from the package root.\n'
 p.write_text(t)
print('Wrote current release notes after passed numerical and visual review.')
