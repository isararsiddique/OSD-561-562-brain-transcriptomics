#!/usr/bin/env python3
"""Collect every result table into one Excel workbook, plus per-contrast CSVs.

The CSVs written by the analyze_* scripts are the source of truth; this script
only assembles them.  The gene-level table is large, so it is summarised into
the workbook and written in full only as CSV.

Outputs
  results/tables/OSD-682_685_698_699_results.xlsx
  results/tables/csv/DE_<region>_<contrast>.csv   one per contrast
  results/tables/all_contrasts_DEG_summary.csv

Usage:  python scripts/export_results.py
"""
from __future__ import annotations

import sys

import pandas as pd

import common as C

# (sheet name, csv file, short description)
SHEETS = [
    ("README", None, "Contents and caveats"),
    ("cohort_summary", "cohort_summary.csv",
     "ROIs, slides and targets per region"),
    ("sample_metadata", "sample_metadata.csv",
     "One row per ROI: region, group, slide, ROI id, factor values"),
    ("qc_metrics", "qc_metrics.csv",
     "Sequencing, ROI geometry and NegProbe background per ROI"),
    ("design_confounds", "design_confounds.csv",
     "Structural limits of the design (slide confound, replication level)"),
    ("normalisation_check", "normalisation_check.csv",
     "Cross-file comparability of the four Q3-normalised matrices"),
    ("log2FC_orientation", "log2_orientation_check.csv",
     "Empirical resolution of the deposited log2 sign convention"),
    ("deposited_stats_audit", "deposited_stats_audit.csv",
     "Deposited adjusted p-values vs correct Benjamini-Hochberg FDR"),
    ("deposited_vs_derived", "deposited_vs_derived.csv",
     "Reproduction of the deposited fold changes by this pipeline"),
    ("PCA_variance", "PCA_variance.csv", "Variance explained per component"),
    ("PCA_coordinates", "PCA_coordinates.csv", "PC1-PC4 per ROI"),
    ("variance_partition", "variance_partition.csv",
     "Share of per-target variance by factor"),
    ("DEG_summary_all", "DEG_summary_all.csv",
     "Target counts per region x contrast under each significance definition"),
    ("sensitivity_MDE", "sensitivity_minimum_detectable_effect.csv",
     "Smallest detectable fold change at 80% power"),
    ("attenuation_summary", "attenuation_summary.csv",
     "BuOE attenuation test, including the variance-assumption sensitivity"),
    ("interaction_top", "interaction_top_targets.csv",
     "Top 25 interaction targets per region by nominal p"),
    ("target_mapping", "target_panel_mapping.csv",
     "CNS target panel: human ENSID -> mouse ortholog -> GeoMx target, with corrections"),
    ("target_coverage", "target_panel_coverage_by_category.csv",
     "Panel coverage by functional category"),
    ("target_stats", "target_stats_all.csv",
     "Per-target statistics with panel-restricted and genome-wide FDR"),
    ("target_sig_panelFDR", "target_significant_panelFDR.csv",
     "Targets significant at BH<0.05 within the panel"),
    ("target_setenrichment", "target_setenrichment.csv",
     "Panel vs transcriptome, rank and label-permutation tests"),
    ("target_categories", "target_category_summary.csv",
     "Functional-category aggregate effects per region x contrast"),
    ("target_attenuation", "target_attenuation.csv",
     "BuOE attenuation restricted to panel targets"),
    ("target_recurrence", "target_recurrence.csv",
     "Panel targets nominal in two or more regions"),
    ("ev_marker_status", "ev_marker_status.csv",
     "CD9/CD63/CD81 mapping and measurability - NOT part of the CNS panel"),
    ("ev_marker_results", "ev_marker_results.csv",
     "CD9/CD63/CD81 tissue mRNA response per region x contrast"),
    ("recurrent_targets", "recurrent_targets.csv",
     "Transcriptome-wide cross-region recurrence"),
    ("recurrence_null", "recurrence_null_test.csv",
     "Permutation null for cross-region recurrence"),
]

README_LINES = [
    ("OSD-682 / 685 / 698 / 699 -- GeoMx DSP re-analysis", ""),
    ("", ""),
    ("Source", "Four OSDR accessions = four brain regions of one GeoMx DSP study."),
    ("", "Processed layer from NCBI GEO GSE239336; OSDR hosts only raw FASTQ."),
    ("Design", "2x2 factorial: Spaceflight (Flight/Ground) x Treatment (BuOE/saline),"),
    ("", "n = 3 ROIs per cell per region, 48 ROIs, 15,782 targets."),
    ("", ""),
    ("READ THIS FIRST", ""),
    ("1", "The adjusted p-values deposited with GSE239336 are unusable: they are"),
    ("", "rounded to two decimals and are smaller than the corresponding raw"),
    ("", "p-value for about 77% of targets, which no correction can produce."),
    ("", "See sheet 'deposited_stats_audit'."),
    ("2", "After correct Benjamini-Hochberg control, NO target reaches FDR < 0.05"),
    ("", "in any region or contrast. Counts in 'DEG_summary_all' labelled"),
    ("", "'nominal' are exploratory and are not multiplicity-controlled."),
    ("3", "With n = 3 the smallest detectable effect at 80% power is a ~1.6-fold"),
    ("", "change nominally and a ~6-8-fold change genome-wide. See"),
    ("", "'sensitivity_MDE'. Absence of hits reflects power, not absence of biology."),
    ("4", "Each experimental group occupies its own slide, so batch and treatment"),
    ("", "cannot be separated. See 'design_confounds'."),
    ("5", "The BuOE attenuation conclusion depends on whether one pooled variance"),
    ("", "or each arm's own variance is used. Both are reported in"),
    ("", "'attenuation_summary'; the arm-specific version is the primary result."),
    ("8", "This is a CNS injury / neuroinflammation / neurovascular panel, NOT an"),
    ("", "extracellular-vesicle panel. CD9/CD63/CD81 are not in it; they are"),
    ("", "reported separately in 'ev_marker_status' and 'ev_marker_results'."),
    ("", "GeoMx measures tissue mRNA, which is not EV cargo."),
    ("6", "The CNS target panel is analysed separately. Because the panel is"),
    ("", "pre-specified, FDR is applied within its 111 measurable members rather"),
    ("", "than across 15,782 targets, and results DO reach significance there."),
    ("", "See 'target_sig_panelFDR'. The 'padj_genomewide' column in"),
    ("", "'target_stats' shows what the same tests give transcriptome-wide."),
    ("7", "Six ENSIDs in the source spreadsheet pointed to unrelated genes and were"),
    ("", "corrected against Ensembl. See 'target_mapping' columns"),
    ("", "ensid_supplied / ensid / correction_reason / correction_verified."),
]


def main() -> int:
    C.ensure_dirs()

    long_path = C.CSV_DIR / "DE_all_contrasts_long.csv"
    if not long_path.exists():
        print(f"ERROR: {long_path} missing -- run scripts/analyze_de.py first",
              file=sys.stderr)
        return 1
    long = pd.read_csv(long_path)

    # Per-contrast CSVs.
    n = 0
    for (region, contrast), g in long.groupby(["region_short", "contrast"], sort=False):
        C.write_csv(
            g.sort_values("pvalue"), f"DE_{region}_{contrast}.csv"
        )
        n += 1
    print(f"  wrote {n} per-contrast gene-level CSVs")

    summary = pd.read_csv(C.CSV_DIR / "DEG_summary_all.csv")
    summary.to_csv(C.TAB_DIR / "all_contrasts_DEG_summary.csv", index=False)
    print(f"  [csv]  {(C.TAB_DIR / 'all_contrasts_DEG_summary.csv').relative_to(C.ROOT)}")

    out = C.TAB_DIR / "OSD-682_685_698_699_results.xlsx"
    written, missing = [], []
    with pd.ExcelWriter(out, engine="openpyxl") as xl:
        pd.DataFrame(README_LINES, columns=["item", "detail"]).to_excel(
            xl, sheet_name="README", index=False
        )
        written.append("README")
        for sheet, fname, _desc in SHEETS:
            if fname is None:
                continue
            path = C.CSV_DIR / fname
            if not path.exists():
                missing.append(fname)
                continue
            df = pd.read_csv(path)
            df.to_excel(xl, sheet_name=sheet[:31], index=False)
            written.append(sheet)

        # A compact gene-level extract: everything at nominal p < 0.01.
        top = long[long["pvalue"] < 0.01].sort_values(["region_short", "contrast", "pvalue"])
        top.to_excel(xl, sheet_name="nominal_p001_targets", index=False)
        written.append("nominal_p001_targets")

    print(f"  [xlsx] {out.relative_to(C.ROOT)}  ({len(written)} sheets)")
    if missing:
        print(f"  WARNING: skipped missing tables: {', '.join(missing)}")
        return 1
    print(f"  sheets: {', '.join(written)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
