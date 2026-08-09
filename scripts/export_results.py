#!/usr/bin/env python3
"""Export all analysis results as CSV files + one multi-sheet Excel workbook.

Uses the full-data contrast logic from analyze.py (all_unique_contrasts) so the
exported numbers cover EVERY unique comparison in the DE tables — Spaceflight,
Age, Environment and Confounded — not just the five matched spaceflight ones.

Outputs
  results/tables/csv/*.csv                      one file per result table
                                                + one full DE CSV per contrast
  results/tables/OSD-561_562_results.xlsx       consolidated workbook

Run:
  ./.venv/bin/python scripts/export_results.py
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

import analyze as A  # same folder

ROOT = Path(__file__).resolve().parent.parent
TAB = ROOT / "results" / "tables"
CSV = TAB / "csv"
XLSX = TAB / "OSD-561_562_results.xlsx"


def gather() -> dict:
    cohort, meta_all, pca_coords, pca_var = [], [], [], []
    deg_summary, de_long, focus_long = [], [], []
    lfc_blocks = []
    vsts, metas = {}, {}
    full_de = {}                       # (study, tag) -> full DE frame
    focus_up = {g.upper() for g in A.read_focus_genes()}

    for st in A.STUDIES:
        info = A.STUDIES[st]
        vst = A.load_counts(st, "VST_Counts")
        dge = A.load_dge(st)
        meta = A.sample_metadata(vst.columns)
        vsts[st], metas[st] = vst, meta

        cohort.append({
            "study": st, "region": info["label"], "glds": info["glds"],
            "samples": int(vst.shape[1]), "genes_measured": int(vst.shape[0]),
            "space_flight": int((meta["treatment"] == "FLT").sum()),
            "ground_control": int((meta["treatment"] == "GC").sum()),
            "old_29wk": int((meta["age"] == "OLD").sum()),
            "young_12wk": int((meta["age"] == "YNG").sum()),
            "iss_t": int((meta["environment"] == "ISS-T").sum()),
            "lar": int((meta["environment"] == "LAR").sum()),
        })

        md = meta.copy()
        md.insert(0, "study", st)
        md.insert(1, "region_label", info["label"])
        meta_all.append(md.reset_index().rename(columns={"index": "sample"}))

        # per-study PCA
        frame, evr = A._pca_frame(vst, A.N_TOP_VARIABLE)
        frame.index.name = "sample"
        pc = frame.join(meta).reset_index()
        pc.insert(0, "pca_set", st)
        pca_coords.append(pc)
        pca_var.append({"pca_set": st, "region": info["label"],
                        "PC1_pct": round(float(evr[0]) * 100, 2),
                        "PC2_pct": round(float(evr[1]) * 100, 2)})

        # ---- ALL unique contrasts ----
        contrasts = A.all_unique_contrasts(dge)
        lfc_cols = {}
        for c in contrasts:
            s = A.contrast_stats_oriented(dge, c["num"], c["den"])
            counts = A.deg_counts(s)
            deg_summary.append({
                "study": st, "region": info["label"], "category": c["category"],
                "label": c["label"], "factors_differ": c["factors_differ"],
                "num_group": c["num"], "den_group": c["den"],
                "genes_tested": counts["genes_tested"],
                "up_in_num": counts["up"], "down_in_num": counts["down"],
                "total_DEG": counts["total"],
            })
            block = s.copy()
            block.insert(0, "label", c["label"])
            block.insert(0, "category", c["category"])
            block.insert(0, "region", info["label"])
            block.insert(0, "study", st)
            de_long.append(block)
            full_de[(st, c["tag"])] = block
            # focus-gene log2FC only for clean single-factor contrasts
            if c["category"] != "Confounded":
                lfc_cols[f"{info['label'][:2]}: {c['label']}"] = (
                    s.set_index("SYMBOL")["log2FC"].groupby(level=0).first())

        mat = pd.DataFrame(lfc_cols)
        mat = mat.loc[mat.index.to_series().str.upper().isin(focus_up)]
        # preserve focus-list ordering
        order = [g for g in A.read_focus_genes() if g in mat.index]
        lfc_blocks.append(mat.loc[order])

    # ---- combined PCA + cohort ----
    combined = pd.concat(vsts.values(), axis=1, join="inner")
    meta_comb = pd.concat(metas.values())
    frame, evr = A._pca_frame(combined, A.N_TOP_VARIABLE)
    frame.index.name = "sample"
    pc = frame.join(meta_comb).reset_index()
    pc.insert(0, "pca_set", "Combined")
    pca_coords.append(pc)
    pca_var.append({"pca_set": "Combined", "region": "Cb+HPC",
                    "PC1_pct": round(float(evr[0]) * 100, 2),
                    "PC2_pct": round(float(evr[1]) * 100, 2)})
    cohort.append({
        "study": "COMBINED", "region": "Cb+HPC", "glds": "-",
        "samples": int(combined.shape[1]), "genes_measured": int(combined.shape[0]),
        "space_flight": int((meta_comb["treatment"] == "FLT").sum()),
        "ground_control": int((meta_comb["treatment"] == "GC").sum()),
        "old_29wk": int((meta_comb["age"] == "OLD").sum()),
        "young_12wk": int((meta_comb["age"] == "YNG").sum()),
        "iss_t": int((meta_comb["environment"] == "ISS-T").sum()),
        "lar": int((meta_comb["environment"] == "LAR").sum()),
    })

    de_all = pd.concat(de_long, ignore_index=True)
    sig_degs = de_all[(de_all["padj"] < A.ADJP_THRESH)
                      & (de_all["log2FC"].abs() >= A.LFC_THRESH)].copy()
    sig_degs = sig_degs.sort_values(["study", "category", "label", "padj"])

    focus_stats = de_all[de_all["SYMBOL"].str.upper().isin(focus_up)].copy()
    focus_stats = focus_stats.sort_values(["study", "category", "label", "SYMBOL"])

    lfc_matrix = pd.concat(lfc_blocks, axis=1)
    lfc_matrix.index.name = "gene"

    deg_df = pd.DataFrame(deg_summary)
    deg_df["_o"] = deg_df["category"].map(
        {"Spaceflight": 0, "Age": 1, "Environment": 2, "Confounded": 3})
    deg_df = deg_df.sort_values(["study", "_o", "total_DEG"],
                                ascending=[True, True, False]).drop(columns="_o")

    return {
        "cohort_summary": pd.DataFrame(cohort),
        "sample_metadata": pd.concat(meta_all, ignore_index=True),
        "PCA_coordinates": pd.concat(pca_coords, ignore_index=True),
        "PCA_variance": pd.DataFrame(pca_var),
        "DEG_summary_all": deg_df,
        "significant_DEGs_all": sig_degs,
        "focus_gene_stats_all": focus_stats,
        "focus_log2FC_matrix": lfc_matrix.reset_index(),
        "_de_all": de_all,
        "_full_de": full_de,
    }


README_ROWS = [
    ("cohort_summary", "Sample & gene counts per study and combined (Flight/Ground, Old/Young, ISS-T/LAR)."),
    ("sample_metadata", "Every sample with region, treatment, age, environment, animal ID."),
    ("PCA_coordinates", "PC1/PC2 per sample for each PCA set (OSD-561, OSD-562, Combined)."),
    ("PCA_variance", "Percent variance explained by PC1/PC2 for each PCA set."),
    ("DEG_summary_all", "DEG counts for ALL 31 unique contrasts, tagged Spaceflight/Age/Environment/Confounded (FDR<0.05 & |log2FC|>=1)."),
    ("significant_DEGs_all", "Gene-level significant DEGs across every unique contrast."),
    ("focus_gene_stats_all", "log2FC/p/adj-p for the focus-gene panel across every unique contrast."),
    ("focus_log2FC_matrix", "Focus genes (rows) x clean single-factor contrasts (cols) log2FC matrix."),
    ("[CSV] DE_all_contrasts_long", "Full gene-level DE for every unique contrast (CSV only)."),
    ("[CSV] DE_<study>_<category>_<contrast>", "Full gene-level DE, one CSV per unique contrast."),
]


def main() -> int:
    CSV.mkdir(parents=True, exist_ok=True)
    d = gather()

    sheet_tables = {k: v for k, v in d.items() if not k.startswith("_")}
    for name, df in sheet_tables.items():
        df.to_csv(CSV / f"{name}.csv", index=False)
    d["_de_all"].to_csv(CSV / "DE_all_contrasts_long.csv", index=False)
    for (st, tag), df in d["_full_de"].items():
        df.to_csv(CSV / f"DE_{st}_{tag}.csv", index=False)

    readme = pd.DataFrame(README_ROWS, columns=["sheet / file", "description"])
    with pd.ExcelWriter(XLSX, engine="openpyxl") as xw:
        readme.to_excel(xw, sheet_name="README", index=False)
        for name, df in sheet_tables.items():
            df.to_excel(xw, sheet_name=name[:31], index=False)
        for ws in xw.book.worksheets:
            for col in ws.columns:
                width = max((len(str(c.value)) for c in col if c.value is not None),
                            default=10)
                ws.column_dimensions[col[0].column_letter].width = min(48, width + 2)

    n_csv = len(list(CSV.glob("*.csv")))
    print(f"[ok] {n_csv} CSV files -> {CSV.relative_to(ROOT)}")
    print(f"[ok] Excel workbook   -> {XLSX.relative_to(ROOT)}")
    print("\nWorkbook sheets:")
    for name, df in sheet_tables.items():
        print(f"  - {name:24s} {df.shape[0]:>7,} rows x {df.shape[1]} cols")
    print(f"\nFull gene-level DE: DE_all_contrasts_long.csv "
          f"({d['_de_all'].shape[0]:,} rows) + {len(d['_full_de'])} per-contrast files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
