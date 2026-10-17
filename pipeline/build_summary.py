"""Write a concise scientific report from the current run's tables."""
import json
import os
from pathlib import Path
import pandas as pd
from pipeline.paths import stage_root, results_root, publication_font

C, E, D = (stage_root(s) for s in ["core", "extended", "deep"])
root = results_root()

def markdown(frame):
    columns = frame.columns.tolist()
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    for row in frame.itertuples(index=False, name=None):
        lines.append("| " + " | ".join(str(x).replace("|", "/").replace("\n", " ") for x in row) + " |")
    return "\n".join(lines)

summary = []
families = [
    ("Original bulk / DESeq2", C / "tables/bulk_panel_complete.csv", "q_panel_primary_full_family", None),
    ("Original GeoMx", C / "tables/spatial_panel_complete.csv", "q_panel_primary_full_family", None),
    ("Bulk method comparison", E / "tables/bulk_method_comparison_panel.csv", "q_full_family", ["method"]),
    ("Independent bulk", E / "tables/independent_bulk_panel.csv", "q_full_family", ["method"]),
    ("RR8 enriched populations", E / "tables/sorted_cell_panel_complete.csv", "q_full_3624_family", ["QC", "method"]),
    ("Additional and regional bulk", D / "tables/additional_panel_all_regions.csv", "q_expanded_1812_family", ["method"]),
    ("Article 5% mitochondrial gate", D / "tables/article_5pct_panel_complete.csv", "q_full_3624_family", ["method"]),
]
for label, path, qcol, keys in families:
    data = pd.read_csv(path)
    groups = data.groupby(keys) if keys else [("", data)]
    for key, sub in groups:
        summary.append({"Analysis": label, "Specification": str(key), "Planned tests": len(sub),
                        "q < 0.05": int(sub[qcol].lt(.05).sum()), "Minimum q": f"{sub[qcol].min():.6g}"})
pd.DataFrame(summary).to_csv(root / "metadata/analysis_family_summary.csv", index=False)
consensus = pd.read_csv(D / "tables/forebrain_all_eight_QC_method_consensus.csv")
different_genes = consensus.loc[consensus.all_eight_corrected, "mouse_symbol_primary"].nunique()
bulk = pd.read_csv(D / "tables/additional_panel_all_regions.csv")
focus = bulk[bulk.analyte.isin(["AIM2", "PYCARD", "SLC1A3", "S100B"]) & bulk.contrast.isin(
    ["flight_vs_habitat_ground", "flight_vs_vivarium", "flight_WT", "flight_balanced_genotype", "flight_balanced_age", "flight_vs_ground"])].copy()
focus = focus[["OSD", "analyte", "method", "contrast", "log2FC", "q_expanded_1812_family"]]
focus["log2FC"] = focus.log2FC.map(lambda x: f"{x:+.3f}" if pd.notna(x) else "unavailable")
focus["q_expanded_1812_family"] = focus.q_expanded_1812_family.map(lambda x: f"{x:.5g}")
checks = [json.loads((C / "tables/verification_checks.json").read_text()),
          json.loads((E / "tables/verification_summary.json").read_text()),
          json.loads((E / "tables/figure_zero_value_audit.json").read_text()),
          json.loads((D / "tables/deep_verification.json").read_text())]
counts = [checks[0]["n_checks"], *[r["checks"] for r in checks[1:]]]
report = f"""# CNS transcriptomics: reproducible analysis report

This run analyzes the supplied 151-target CNS panel against deposited count matrices and verified biological identities. The complete planned testing families are retained, with unavailable effects left missing and their hypothesis p-values set to one. Every numerical summary below comes from this output folder.

{different_genes} distinct mapped RR8 forebrain genes pass q < 0.05 in all eight method/QC combinations. Those combinations reuse the same eight mice and do not add independent biological replication.

## Complete testing families

{markdown(pd.DataFrame(summary))}

## Selected biological contrasts

{markdown(focus)}

All target effects, including uncorrected and unavailable estimates, are in the full CSV tables. Mouse-deletion results are in [the expanded deletion summary](deep/tables/reported_expanded_deletion_summary.csv). Expression-state, technical-state, fixed-composition, depth, paired-region, pathway, spatial and promoter checks retain their original families and experimental units.

## Interpretation limits

- Original terminal bulk preservation differs by flight group. RR8 exposure and processing day are confounded. RR23 collection dates and habitat/vivarium comparator differences matter.
- Cells, nuclei, technical libraries, spatial spots, sections and shared modalities do not increase the number of exposure animals. Original GeoMx animal links remain unresolved.
- Expression states are numerical exploratory clusters learned and tested in the same mice. RNA fractions are not absolute cell fractions; adjusted associations are not causal mediation.
- CAMERA results depend on intergene correlation and residual degrees of freedom. An uncorrected test or zero discovery count is not evidence of equivalence.
- These data measure transcripts. They do not establish EV cargo, protein activation, inflammasome function, protection or a first-ever mechanism. Protein-source scale flags remain recorded.

## Verification and outputs

The independent core, extension, figure-value and deeper audits completed {sum(counts)} checks ({' + '.join(map(str, counts))}). Source integrity is recorded in [the input manifest](metadata/input_manifest.json), package versions in [the environment](metadata/environment.json), and file hashes in [the output manifest](metadata/output_manifest.json).

Figures are in [core](core/figures/), [integrated analyses](extended/figures/) and [deeper validation](deep/figures/). Publication figures use black headings and the installed {publication_font()} font, embedded in PDFs. Every stage has a separate log. The complete script order is available with `python run_pipeline.py list`.
"""
(root / "ANALYSIS_REPORT.md").write_text(report)
print(f"Built report from this run's tables; {sum(counts)} numerical checks.")
