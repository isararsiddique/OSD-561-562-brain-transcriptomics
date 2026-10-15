# Portable transcript-context audit

This is a full-universe reference-annotation and descriptive read-support audit, completed after the original statistical findings were known. Read `SCIENTIFIC_REPORT.md` and `COMPLETION_AND_VERIFICATION_PLAN_2026-10-15.md`. It does not add independent expression, protein, full-transcript or causal-flight evidence.

From this directory, using Python with NumPy, pandas, SciPy and Matplotlib:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
TC_PYTHON=/path/to/python
TC_PACKAGE=../..
"$TC_PYTHON" scripts/analyse_context.py --base "$TC_PACKAGE"
"$TC_PYTHON" scripts/verify_context.py --base "$TC_PACKAGE"
"$TC_PYTHON" scripts/build_context_outputs.py
```

All writes resolve to this script-relative extension. The NEW package supplies `inputs/all_animal_junction_counts.csv.gz`, animal metadata, sibling multiway definitions/compositions/local results, and gene omnibus event/candidate tables. The frozen local `inputs/required_gene_context.GRCm39.112.gtf.gz` is sufficient for offline reproduction; network access or the original full GTF is unnecessary. `prepare_context.py` is the separate original initialization step requiring the full source GTF and publisher directory; do not rerun it for routine offline reproduction. Its original provenance is preserved.

The verifier checks immutable local GTF/publisher hashes and current input semantics. Parent result CSV.gz bytes can change when regenerated due to gzip headers/timestamps. Their historical byte hashes are therefore recorded but not used to reject semantically identical tables.

Key outputs include the full 923-cluster and 1,983-gene audits; all 3,248 multiway junction match statuses (including 1,199 unmatched); all exact transcript memberships and CDS/UTR summaries; three retained blank-symbol gene IDs; all Tcaf1 animal counts, conditional fractions and companion intron support; independent verification; and scientific report.

`figures/Figure_Transcript_Context.pdf`, `.svg`, `.png` and `Figure_Transcript_Context_caption.txt` show all six matching transcript memberships, common coding structure and actual 26-animal read support. The figure is 183 mm wide, uses Arial/Helvetica at 8 pt and ≥1 pt lines, and the PNG is 600 dpi. The caption is ≤350 words. All deliverables except the manifest itself are hashed in `MANIFEST.csv`; software versions are recorded in `tables/software_versions.json`.
