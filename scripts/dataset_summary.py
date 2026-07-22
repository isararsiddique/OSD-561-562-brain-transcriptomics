#!/usr/bin/env python3
"""Compute the real dataset / analysis numbers used to design the graphical
abstract. Reuses the loaders and contrast logic from analyze.py so the numbers
match the figures exactly.

Writes:
  results/tables/dataset_summary.json   (machine-readable)
  results/tables/dataset_summary.md      (human-readable)
and prints a summary to stdout.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

import analyze as A  # same folder; reuse loaders + contrast detection

ROOT = Path(__file__).resolve().parent.parent
TAB = ROOT / "results" / "tables"


def counts_dict(series) -> dict:
    return {str(k): int(v) for k, v in series.value_counts().sort_index().items()}


def study_block(study: str) -> dict:
    info = A.STUDIES[study]
    vst = A.load_counts(study, "VST_Counts")
    dge = A.load_dge(study)
    meta = A.sample_metadata(vst.columns)

    frame, evr = A._pca_frame(vst, A.N_TOP_VARIABLE)

    # DEGs per spaceflight contrast (FDR<0.05 & |LFC|>=1), whole transcriptome
    deg = []
    for c in A.spaceflight_contrasts(dge):
        s = A.contrast_stats(dge, c["flight"], c["ground"]).dropna(
            subset=["log2FC", "padj"])
        sig = (s["padj"] < A.ADJP_THRESH) & (s["log2FC"].abs() >= A.LFC_THRESH)
        up = int((sig & (s["log2FC"] > 0)).sum())
        down = int((sig & (s["log2FC"] < 0)).sum())
        deg.append({"contrast": c["short"], "genes_tested": int(s.shape[0]),
                    "up_in_flight": up, "down_in_flight": down, "total_DEG": up + down})

    return {
        "label": info["label"],
        "glds": info["glds"],
        "n_samples": int(vst.shape[1]),
        "n_genes_measured": int(vst.shape[0]),
        "by_treatment": counts_dict(meta["treatment_full"]),
        "by_age": counts_dict(meta["age"]),
        "by_environment": counts_dict(meta["environment"]),
        "pca_pc1_pct": round(float(evr[0]) * 100, 1),
        "pca_pc2_pct": round(float(evr[1]) * 100, 1),
        "n_spaceflight_contrasts": len(deg),
        "degs": deg,
    }


def combined_block(studies: list[str]) -> dict:
    vsts, metas = [], []
    for st in studies:
        v = A.load_counts(st, "VST_Counts")
        vsts.append(v)
        metas.append(A.sample_metadata(v.columns))
    combined = pd.concat(vsts, axis=1, join="inner")
    meta = pd.concat(metas)
    frame, evr = A._pca_frame(combined, A.N_TOP_VARIABLE)
    return {
        "n_samples": int(combined.shape[1]),
        "n_shared_genes": int(combined.shape[0]),
        "by_region": counts_dict(meta["region"]),
        "pca_pc1_pct": round(float(evr[0]) * 100, 1),
        "pca_pc2_pct": round(float(evr[1]) * 100, 1),
    }


def focus_highlights() -> list[dict]:
    path = TAB / "focus_gene_stats_by_contrast.csv"
    if not path.exists():
        return []
    df = pd.read_csv(path)
    df["abs"] = df["log2FC"].abs()
    top = df.sort_values("abs", ascending=False).head(8)
    return [
        {"study": r.study, "contrast": r.contrast, "gene": r.SYMBOL,
         "log2FC": round(float(r.log2FC), 2),
         "pvalue": round(float(r.pvalue), 4),
         "padj": (None if pd.isna(r.padj) else round(float(r.padj), 4))}
        for r in top.itertuples()
    ]


def main() -> int:
    studies = list(A.STUDIES)
    summary = {
        "studies": {st: study_block(st) for st in studies},
        "combined": combined_block(studies),
        "focus_genes": A.read_focus_genes(),
        "focus_highlights": focus_highlights(),
        "thresholds": {"FDR": A.ADJP_THRESH, "abs_log2FC": A.LFC_THRESH,
                       "pca_top_variable_genes": A.N_TOP_VARIABLE},
    }

    TAB.mkdir(parents=True, exist_ok=True)
    (TAB / "dataset_summary.json").write_text(json.dumps(summary, indent=2))
    _write_markdown(summary)
    print(json.dumps(summary, indent=2))
    print(f"\n[written] {TAB/'dataset_summary.json'}\n[written] {TAB/'dataset_summary.md'}")
    return 0


def _write_markdown(s: dict) -> None:
    lines = ["# Dataset summary (auto-generated)\n"]
    c = s["combined"]
    lines.append(f"- Combined: **{c['n_samples']} samples**, "
                 f"**{c['n_shared_genes']:,} shared genes**, "
                 f"region split {c['by_region']}, "
                 f"PCA PC1 {c['pca_pc1_pct']}% / PC2 {c['pca_pc2_pct']}%\n")
    for st, b in s["studies"].items():
        lines.append(f"\n## {st} — {b['label']} ({b['glds']})\n")
        lines.append(f"- {b['n_samples']} samples, {b['n_genes_measured']:,} genes measured")
        lines.append(f"- Treatment: {b['by_treatment']}")
        lines.append(f"- Age: {b['by_age']}")
        lines.append(f"- Environment: {b['by_environment']}")
        lines.append(f"- PCA: PC1 {b['pca_pc1_pct']}% / PC2 {b['pca_pc2_pct']}%")
        lines.append(f"- Spaceflight DEGs (FDR<{s['thresholds']['FDR']}, "
                     f"|log2FC|>={s['thresholds']['abs_log2FC']}):")
        lines.append("\n  | Contrast | Genes tested | Up | Down | Total |")
        lines.append("  |---|---|---|---|---|")
        for d in b["degs"]:
            lines.append(f"  | {d['contrast']} | {d['genes_tested']:,} | "
                         f"{d['up_in_flight']} | {d['down_in_flight']} | {d['total_DEG']} |")
    (TAB / "dataset_summary.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    raise SystemExit(main())
