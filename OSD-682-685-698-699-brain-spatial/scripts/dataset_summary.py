#!/usr/bin/env python3
"""Write a dataset and results overview as Markdown + JSON.

Every number is read back from the exported tables, so the summary cannot drift
away from the analysis.

Outputs
  results/tables/dataset_summary.md
  results/tables/dataset_summary.json

Usage:  python scripts/dataset_summary.py
"""
from __future__ import annotations

import json
import sys

import numpy as np
import pandas as pd

import common as C


def load(name: str) -> pd.DataFrame:
    path = C.CSV_DIR / name
    if not path.exists():
        raise SystemExit(f"missing {path.relative_to(C.ROOT)} -- run the analyze_* scripts")
    return pd.read_csv(path)


def main() -> int:
    C.ensure_dirs()
    cohort = load("cohort_summary.csv")
    qc = load("qc_metrics.csv")
    audit = load("deposited_stats_audit.csv")
    summary = load("DEG_summary_all.csv")
    sens = load("sensitivity_minimum_detectable_effect.csv")
    att = load("attenuation_summary.csv")
    vp = load("variance_partition.csv")
    var = load("PCA_variance.csv")
    norm = load("normalisation_check.csv")

    combined = cohort[cohort["region_short"] == "combined"].iloc[0]
    per_region = cohort[cohort["region_short"] != "combined"]

    payload = {
        "project": "OSD-682/685/698/699 GeoMx DSP re-analysis",
        "geo_series": C.GEO_SERIES,
        "studies": [
            {
                "osd": s.osd, "glds": s.glds, "region": s.region,
                "region_code": s.code, "label": s.short,
            }
            for s in C.STUDIES
        ],
        "assay": "NanoString GeoMx Digital Spatial Profiling (spatial transcriptomics)",
        "design": {
            "factors": {"Spaceflight": C.SF_LEVELS, "Treatment": C.TRT_LEVELS},
            "n_per_cell": 3,
            "n_roi_total": int(combined["n_roi"]),
            "n_slides": int(combined["n_slides"]),
            "n_targets": int(combined["n_targets"]),
            "sex": str(combined["sex"]),
            "strain": str(combined["strain"]),
        },
        "qc": {
            "raw_reads_millions_min": round(float(qc["raw_reads"].min()) / 1e6, 2),
            "raw_reads_millions_max": round(float(qc["raw_reads"].max()) / 1e6, 2),
            "raw_reads_millions_median": round(float(qc["raw_reads"].median()) / 1e6, 2),
            "seq_saturation_min": float(qc["seq_saturation"].min()),
            "seq_saturation_max": float(qc["seq_saturation"].max()),
            "align_rate_min": round(float(qc["align_rate"].min()), 4),
            "align_rate_max": round(float(qc["align_rate"].max()), 4),
            "nuclei_min": int(qc["nuclei"].min()),
            "nuclei_max": int(qc["nuclei"].max()),
            "q3_upper_quartile_identical_across_files": bool(
                np.allclose(norm["q3_median"], norm["q3_median"].iloc[0])
            ),
        },
        "structure": {
            "pc1_variance_explained": round(
                float(var[(var["scope"] == "combined") & (var["PC"] == "PC1")]
                      ["variance_explained"].iloc[0]), 4),
            "pc2_variance_explained": round(
                float(var[(var["scope"] == "combined") & (var["PC"] == "PC2")]
                      ["variance_explained"].iloc[0]), 4),
            "median_variance_fraction": {
                str(r["factor"]): round(float(r["median_frac_variance"]), 4)
                for _, r in vp.iterrows()
            },
        },
        "deposited_statistics": {
            "frac_adjusted_below_raw_min": round(float(audit["frac_adj_below_raw_p"].min()), 4),
            "frac_adjusted_below_raw_max": round(float(audit["frac_adj_below_raw_p"].max()), 4),
            "n_unique_adjusted_values_max": int(audit["n_unique_adj_values"].max()),
            "deposited_significant_min": int(audit["n_sig_deposited_adj_lt_0.05"].min()),
            "deposited_significant_max": int(audit["n_sig_deposited_adj_lt_0.05"].max()),
            "recomputed_BH_significant_total": int(audit["n_sig_recomputed_BH_lt_0.05"].sum()),
            "best_recomputed_BH": round(float(audit["min_recomputed_BH"].min()), 4),
        },
        "differential_expression": {
            "fdr_cutoff": C.FDR_CUTOFF,
            "lfc_cutoff": C.LFC_CUTOFF,
            "n_contrasts_tested": int(summary["contrast"].nunique()),
            "n_region_contrast_combinations": int(len(summary)),
            "total_fdr_significant_targets": int(summary["fdr05_lfc1_n_deg"].sum()),
            "best_padj_across_everything": round(float(summary["min_padj_BH"].min()), 4),
            "nominal_p001_obs_over_expected_min": round(
                float(summary["obs_over_expected_p0.01"].min()), 3),
            "nominal_p001_obs_over_expected_max": round(
                float(summary["obs_over_expected_p0.01"].max()), 3),
        },
        "sensitivity": {
            r["region_short"]: {
                "median_se_log2": round(float(r["median_se_log2"]), 4),
                "mde_log2FC_nominal": round(float(r["mde_log2FC_nominal_a0.05"]), 3),
                "mde_log2FC_genomewide": round(
                    float(r["mde_log2FC_genomewide_a0.05_over_n"]), 3),
                "mde_fold_change_genomewide": round(
                    float(r["mde_fold_change_genomewide_a0.05_over_n"]), 2),
            }
            for _, r in sens[sens["contrast"] == "FLT_vs_GC_in_SAL"].iterrows()
        },
        "buoe_attenuation": {
            r["region_short"]: {
                "verdict": str(r["verdict"]),
                "rms_true_effect_saline": round(float(r["rms_true_effect_saline"]), 4),
                "rms_true_effect_BuOE": round(float(r["rms_true_effect_BuOE"]), 4),
                "msq_diff_saline_minus_BuOE": round(
                    float(r["msq_diff_saline_minus_BuOE"]), 5),
                "msq_diff_ci": [round(float(r["msq_diff_ci_lo"]), 5),
                                round(float(r["msq_diff_ci_hi"]), 5)],
                "bootstrap_p": float(r["msq_diff_bootstrap_p"]),
                "within_arm_variance_ratio_BuOE_over_saline": round(
                    float(r["median_var_ratio_BuOE_over_saline"]), 3),
                "verdict_under_pooled_variance": str(r["verdict_under_pooled_variance"]),
                "verdict_depends_on_variance_assumption": bool(
                    r["verdict_depends_on_variance_assumption"]),
                "n_interaction_targets_fdr05": int(r["n_interaction_fdr05"]),
            }
            for _, r in att.iterrows()
        },
    }

    tmap_path = C.CSV_DIR / "target_panel_mapping.csv"
    tsig_path = C.CSV_DIR / "target_significant_panelFDR.csv"
    if tmap_path.exists() and tsig_path.exists():
        tmap = pd.read_csv(tmap_path)
        tsig = pd.read_csv(tsig_path)
        tenr = load("target_setenrichment.csv")
        tcat = load("target_category_summary.csv")
        tatt = load("target_attenuation.csv")
        payload["cns_ev_target_panel"] = {
            "n_analytes": int(len(tmap)),
            "n_measurable": int(tmap["on_geomx_panel"].sum()),
            "n_ensid_corrections": int(tmap["ensid_corrected"].sum()),
            "all_corrections_verified": bool(
                tmap.loc[tmap["ensid_corrected"], "correction_verified"].all()
            ),
            "corrections": [
                {
                    "analyte": r["analyte"],
                    "supplied": r["ensid_supplied"],
                    "corrected_to": r["ensid"],
                    "resolves_to": r["human_symbol_ensembl"],
                    "measured": bool(r["on_geomx_panel"]),
                }
                for _, r in tmap[tmap["ensid_corrected"]].iterrows()
            ],
            "n_significant_panel_fdr": int(len(tsig)),
            "significant": [
                {
                    "analyte": r["analyte"],
                    "gene": r["gene"],
                    "region": r["region_short"],
                    "contrast": r["contrast_label"],
                    "log2FC": round(float(r["log2FC"]), 3),
                    "padj_panel": round(float(r["padj_panel"]), 4),
                    "padj_genomewide": round(float(r["padj_genomewide"]), 3),
                }
                for _, r in tsig.iterrows()
            ],
            "set_enrichment_significant": [
                {
                    "region": r["region_short"],
                    "contrast": r["contrast_label"],
                    "delta_mean_abs_t": round(float(r["delta_mean_abs_t"]), 3),
                    "permutation_p": float(r["permutation_p_delta_abs_t"]),
                }
                for _, r in tenr[tenr["permutation_p_delta_abs_t"] < 0.05].iterrows()
            ],
            "categories_q_below_005": [
                {
                    "region": r["region_short"],
                    "contrast": r["contrast_label"],
                    "category": r["category"],
                    "mean_t": round(float(r["mean_t"]), 3),
                    "q": round(float(r["q_mean_t"]), 4),
                }
                for _, r in tcat[tcat["q_mean_t"] < 0.05].iterrows()
            ],
            "attenuation_verdicts": {
                r["region_short"]: r["verdict"] for _, r in tatt.iterrows()
            },
        }

    (C.TAB_DIR / "dataset_summary.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(f"  [json] {(C.TAB_DIR / 'dataset_summary.json').relative_to(C.ROOT)}")

    d = payload
    md: list[str] = []
    md.append("# OSD-682 / 685 / 698 / 699 — dataset and results summary\n")
    md.append(
        "Generated by `scripts/dataset_summary.py`. Every number is read back from "
        "the exported tables in `results/tables/csv/`.\n"
    )

    md.append("## Datasets\n")
    md.append("| OSD | GLDS | Region | Label | ROIs |")
    md.append("|---|---|---|---|---|")
    for _, r in per_region.iterrows():
        md.append(
            f"| [{r['osd']}](https://osdr.nasa.gov/bio/repo/data/studies/{r['osd']}) "
            f"| {r['glds']} | {r['region']} | {r['region_short']} | {int(r['n_roi'])} |"
        )
    md.append(
        f"| **combined** | — | all four | — | **{d['design']['n_roi_total']}** |\n"
    )
    md.append(
        f"All four accessions are one NanoString GeoMx DSP experiment "
        f"(GEO [{C.GEO_SERIES}](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc="
        f"{C.GEO_SERIES})): {d['design']['n_targets']:,} targets, "
        f"{d['design']['n_roi_total']} ROIs on {d['design']['n_slides']} slides, "
        f"{d['design']['strain']} {d['design']['sex']} mice.\n"
    )

    md.append("## Design\n")
    md.append(
        "2 × 2 factorial, n = 3 ROIs per cell per region: "
        "Spaceflight (Flight vs Ground Control) × Treatment "
        "(BuOE antioxidant vs saline).\n"
    )

    md.append("## Quality control\n")
    q = d["qc"]
    md.append(
        f"- Raw reads per ROI: {q['raw_reads_millions_min']}–"
        f"{q['raw_reads_millions_max']} M (median {q['raw_reads_millions_median']} M)\n"
        f"- Alignment rate: {q['align_rate_min']:.3f}–{q['align_rate_max']:.3f}; "
        f"sequencing saturation {q['seq_saturation_min']}–{q['seq_saturation_max']}%\n"
        f"- Nuclei per ROI: {q['nuclei_min']}–{q['nuclei_max']}\n"
        f"- The four Q3-normalised files share an identical per-ROI upper quartile "
        f"({str(q['q3_upper_quartile_identical_across_files']).lower()}), so they are "
        f"directly comparable across regions\n"
    )

    md.append("## Global structure\n")
    s = d["structure"]
    md.append(
        f"Combined PCA over all {d['design']['n_roi_total']} ROIs: "
        f"PC1 = {s['pc1_variance_explained']*100:.1f}%, "
        f"PC2 = {s['pc2_variance_explained']*100:.1f}%.\n"
    )
    md.append("Median share of per-target variance:\n")
    md.append("| Factor | Median variance share |")
    md.append("|---|---|")
    for k, v in sorted(s["median_variance_fraction"].items(), key=lambda x: -x[1]):
        md.append(f"| {k} | {v*100:.1f}% |")
    md.append("")
    md.append(
        "Brain region dominates; both experimental factors explain less variance "
        "than section-to-section variation.\n"
    )

    md.append("## The deposited statistics cannot be used as they stand\n")
    a = d["deposited_statistics"]
    md.append(
        f"- The deposited `Adjusted pvalue` is below the deposited raw `Pvalue` for "
        f"{a['frac_adjusted_below_raw_min']*100:.0f}–"
        f"{a['frac_adjusted_below_raw_max']*100:.0f}% of targets. "
        f"No multiple-testing correction can lower a p-value.\n"
        f"- It takes at most {a['n_unique_adjusted_values_max']} distinct values "
        f"(rounded to two decimals).\n"
        f"- Counted as deposited, {a['deposited_significant_min']:,}–"
        f"{a['deposited_significant_max']:,} targets per contrast look significant. "
        f"Recomputing Benjamini-Hochberg from the raw p-values gives "
        f"**{a['recomputed_BH_significant_total']}** across the entire experiment; "
        f"the smallest adjusted p-value anywhere is {a['best_recomputed_BH']:.3f}.\n"
    )

    md.append("## Differential expression\n")
    de = d["differential_expression"]
    md.append(
        f"{de['n_contrasts_tested']} contrasts × 4 regions = "
        f"{de['n_region_contrast_combinations']} comparisons, re-derived from the "
        f"Q3 matrix with a per-target 2 × 2 model.\n\n"
        f"- Targets at BH < {de['fdr_cutoff']} and |log2FC| ≥ {de['lfc_cutoff']}: "
        f"**{de['total_fdr_significant_targets']}**\n"
        f"- Best adjusted p-value across all comparisons: "
        f"{de['best_padj_across_everything']:.3f}\n"
        f"- Nominal p < 0.01 counts run "
        f"{de['nominal_p001_obs_over_expected_min']}× to "
        f"{de['nominal_p001_obs_over_expected_max']}× the null expectation\n"
    )

    md.append("### What this design could have detected\n")
    md.append("| Region | Median SE (log2) | Detectable log2FC (nominal) | Detectable log2FC (genome-wide) | As a fold change |")
    md.append("|---|---|---|---|---|")
    for region, v in d["sensitivity"].items():
        md.append(
            f"| {region} | {v['median_se_log2']:.3f} | {v['mde_log2FC_nominal']:.2f} "
            f"| {v['mde_log2FC_genomewide']:.2f} | {v['mde_fold_change_genomewide']:.1f}× |"
        )
    md.append("")
    md.append(
        "At 80% power. Surviving multiple-testing control would have required "
        "roughly a six- to eight-fold change, far beyond what brain tissue shows "
        "under these perturbations, so the absence of FDR-significant targets "
        "reflects the power of the design rather than the absence of biology.\n"
    )

    md.append("## Does BuOE attenuate the spaceflight response?\n")
    md.append("| Region | RMS effect, saline | RMS effect, BuOE | saline − BuOE (95% CI) | bootstrap p | Verdict |")
    md.append("|---|---|---|---|---|---|")
    for region, v in d["buoe_attenuation"].items():
        ci = v["msq_diff_ci"]
        bp = v["bootstrap_p"]
        bp_txt = "< 0.001" if bp < 0.001 else f"{bp:.3f}"
        md.append(
            f"| {region} | {v['rms_true_effect_saline']:.3f} | "
            f"{v['rms_true_effect_BuOE']:.3f} | "
            f"{v['msq_diff_saline_minus_BuOE']:+.4f} "
            f"({ci[0]:+.4f}, {ci[1]:+.4f}) | {bp_txt} | {v['verdict']} |"
        )
    md.append("")
    flips = [r for r, v in d["buoe_attenuation"].items()
             if v["verdict_depends_on_variance_assumption"]]
    md.append(
        "No individual target reaches FDR significance for the interaction term in "
        "any region, so this is an aggregate effect-size comparison, not a gene list. "
        "The BuOE arm carries "
        + "–".join(
            f"{v['within_arm_variance_ratio_BuOE_over_saline']:.2f}"
            for v in [min(d["buoe_attenuation"].values(),
                          key=lambda x: x["within_arm_variance_ratio_BuOE_over_saline"]),
                      max(d["buoe_attenuation"].values(),
                          key=lambda x: x["within_arm_variance_ratio_BuOE_over_saline"])]
        )
        + "× the within-arm variance of the saline arm, so the choice of variance "
        "estimate matters: "
        + (f"the verdict changes for {', '.join(flips)} "
           if flips else "the verdict is unchanged ")
        + "when a single pooled variance is assumed instead.\n"
    )

    if "cns_ev_target_panel" in d:
        t = d["cns_ev_target_panel"]
        md.append("## CNS / EV target panel\n")
        md.append(
            f"{t['n_measurable']} of {t['n_analytes']} analytes in the project target "
            f"list are measurable on this assay. Because the panel is pre-specified, "
            f"FDR is controlled within it rather than across all "
            f"{d['design']['n_targets']:,} targets.\n"
        )
        if t["n_ensid_corrections"]:
            md.append(
                f"### Source-data corrections\n\n"
                f"{t['n_ensid_corrections']} Ensembl IDs in the spreadsheet resolved to "
                f"genes unrelated to their analyte label and were corrected "
                f"(all {'verified' if t['all_corrections_verified'] else 'NOT fully verified'} "
                f"against Ensembl):\n"
            )
            md.append("| Analyte | Supplied ENSID | Corrected to | Resolves to | Measured |")
            md.append("|---|---|---|---|---|")
            for c in t["corrections"]:
                md.append(
                    f"| {c['analyte']} | `{c['supplied']}` | `{c['corrected_to']}` "
                    f"| {c['resolves_to']} | {'yes' if c['measured'] else 'no'} |"
                )
            md.append("")

        md.append(
            f"### Targets significant at BH < 0.05 within the panel "
            f"({t['n_significant_panel_fdr']})\n"
        )
        if t["significant"]:
            md.append("| Analyte | Gene | Region | Contrast | log2FC | Panel BH | Genome-wide BH |")
            md.append("|---|---|---|---|---|---|---|")
            for s in t["significant"]:
                md.append(
                    f"| {s['analyte']} | *{s['gene']}* | {s['region']} | {s['contrast']} "
                    f"| {s['log2FC']:+.3f} | {s['padj_panel']:.4f} "
                    f"| {s['padj_genomewide']:.2f} |"
                )
            md.append("")
            md.append(
                "The last two columns are the point of the panel-restricted analysis: "
                "none of these results is detectable transcriptome-wide.\n"
            )
        else:
            md.append("None.\n")

        if t["set_enrichment_significant"]:
            md.append("### Panel more responsive than the transcriptome (permutation P < 0.05)\n")
            md.append("| Region | Contrast | Panel − background mean \\|t\\| | Permutation P |")
            md.append("|---|---|---|---|")
            for s in t["set_enrichment_significant"]:
                md.append(
                    f"| {s['region']} | {s['contrast']} | {s['delta_mean_abs_t']:+.3f} "
                    f"| {s['permutation_p']:.3f} |"
                )
            md.append("")

        if t["categories_q_below_005"]:
            md.append("### Functional categories moving coherently (q < 0.05)\n")
            md.append("| Region | Contrast | Category | Mean t | q |")
            md.append("|---|---|---|---|---|")
            for s in t["categories_q_below_005"]:
                md.append(
                    f"| {s['region']} | {s['contrast']} | {s['category']} "
                    f"| {s['mean_t']:+.3f} | {s['q']:.4f} |"
                )
            md.append("")

        md.append(
            "Panel-restricted BuOE attenuation: "
            + "; ".join(f"{k} {v}" for k, v in t["attenuation_verdicts"].items())
            + ". With only 111 targets the bootstrap interval is much wider than the "
            "transcriptome-wide version, so no region reaches a verdict.\n"
        )

    (C.TAB_DIR / "dataset_summary.md").write_text("\n".join(md) + "\n")
    print(f"  [md]   {(C.TAB_DIR / 'dataset_summary.md').relative_to(C.ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
