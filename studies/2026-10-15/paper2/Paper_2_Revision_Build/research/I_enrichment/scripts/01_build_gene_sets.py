#!/usr/bin/env python3
"""Build query gene sets and the measured-gene background for Paper 2 enrichment.

Reads ONLY existing gene-level omnibus results from the submission package.
No new statistics on the raw data; gene sets are definitions, not discoveries.

Outputs (research/I_enrichment/tables/):
  background_measured_genes.txt      all measured genes (universe)
  set_FLT_GC_BY.txt                  BY joint q<=0.05 in FLT-GC
  set_FLT_GC_top100.txt              strongest FLT-GC exact evidence, ties kept
  set_GC_VIV_top100.txt              strongest GC-VIV exact evidence, ties kept
  gene_sets_definition.csv           per-gene membership + ranking statistics
  gene_sets_summary.json             sizes, tie handling, provenance
"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path

import pandas as pd

PKG = Path(__file__).resolve().parents[4] / "Paper_2_Complete_Submission_Package"
OMNIBUS = PKG / "extensions/gene_omnibus/results/gene_omnibus_complete.csv.gz"
OUT = Path(__file__).resolve().parents[1] / "tables"
OUT.mkdir(parents=True, exist_ok=True)

TOP_N = 100


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_list(path: Path, genes) -> None:
    path.write_text("\n".join(genes) + "\n")


def top_set(df: pd.DataFrame, contrast: str, n: int) -> tuple[list[str], dict]:
    """Top-n genes by exact p, with ALL ties at the cut-off retained.

    p_exact is discrete (permutation-based), so the n-th p value is usually
    shared by several genes. Breaking those ties with a secondary effect-size
    statistic would mix two different kinds of evidence, so instead the whole
    tie block is kept and the realised set size is reported.
    """
    sub = df[df.contrast == contrast].dropna(subset=["gene_symbol"]).copy()
    sub = sub.sort_values(["p_exact", "gene_symbol"], kind="mergesort")
    if len(sub) <= n:
        cut = sub.p_exact.max()
    else:
        cut = sub.p_exact.iloc[n - 1]
    keep = sub[sub.p_exact <= cut]
    info = {
        "contrast": contrast,
        "requested_n": n,
        "p_cutoff_inclusive": float(cut),
        "realised_n_genes": int(keep.gene_symbol.nunique()),
        "realised_n_rows": int(len(keep)),
        "n_tied_at_cutoff": int((sub.p_exact == cut).sum()),
        "tie_rule": "all genes with p_exact <= p of the n-th ranked gene are kept",
    }
    return sorted(keep.gene_symbol.unique().tolist()), info


def main() -> None:
    om = pd.read_csv(OMNIBUS)
    prov = {
        "source_file": str(OMNIBUS),
        "source_sha256": sha256(OMNIBUS),
        "rows": int(len(om)),
        "contrasts": sorted(om.contrast.unique().tolist()),
    }

    # Background = every gene measured by the junction-level pipeline.
    # Genes without an Ensembl symbol cannot be queried, so they are dropped
    # from the background as well as from the query sets (same universe both
    # sides of the test).
    uni = om.drop_duplicates("gene_id")[["gene_id", "gene_symbol"]]
    n_no_symbol = int(uni.gene_symbol.isna().sum())
    background = sorted(uni.gene_symbol.dropna().unique().tolist())

    fl = om[om.contrast == "FLT_minus_GC"]
    by_genes = sorted(
        fl.loc[fl.q_BY_joint_3G <= 0.05, "gene_symbol"].dropna().unique().tolist()
    )
    by_rows = int((fl.q_BY_joint_3G <= 0.05).sum())

    flt_top, flt_info = top_set(om, "FLT_minus_GC", TOP_N)
    gcv_top, gcv_info = top_set(om, "GC_minus_VIV", TOP_N)

    write_list(OUT / "background_measured_genes.txt", background)
    write_list(OUT / "set_FLT_GC_BY.txt", by_genes)
    write_list(OUT / "set_FLT_GC_top100.txt", flt_top)
    write_list(OUT / "set_GC_VIV_top100.txt", gcv_top)

    # Per-gene definition table (wide, one row per gene).
    wide = om.pivot_table(
        index=["gene_id", "gene_symbol"],
        columns="contrast",
        values=["p_exact", "q_BY_joint_3G", "rms_delta_percentage_points", "event_count"],
        aggfunc="first",
    )
    wide.columns = [f"{a}__{b}" for a, b in wide.columns]
    wide = wide.reset_index()
    wide["in_FLT_GC_BY"] = wide.gene_symbol.isin(by_genes)
    wide["in_FLT_GC_top100"] = wide.gene_symbol.isin(flt_top)
    wide["in_GC_VIV_top100"] = wide.gene_symbol.isin(gcv_top)
    wide.to_csv(OUT / "gene_sets_definition.csv", index=False)

    summary = {
        "provenance": prov,
        "background": {
            "n_measured_gene_ids": int(len(uni)),
            "n_without_gene_symbol_dropped": n_no_symbol,
            "n_background_symbols": len(background),
            "definition": "all genes with >=1 analysed junction cluster in the "
            "Paper 2 gene-level omnibus (the MEASURED universe, not the genome)",
        },
        "sets": {
            "FLT_GC_BY": {
                "n": len(by_genes),
                "n_rows": by_rows,
                "definition": "q_BY_joint_3G <= 0.05 in FLT_minus_GC",
            },
            "FLT_GC_top100": flt_info | {"n": len(flt_top)},
            "GC_VIV_top100": gcv_info | {"n": len(gcv_top)},
        },
        "overlaps": {
            "BY_in_FLT_top": len(set(by_genes) & set(flt_top)),
            "FLT_top_in_GC_VIV_top": len(set(flt_top) & set(gcv_top)),
            "BY_in_GC_VIV_top": len(set(by_genes) & set(gcv_top)),
        },
    }
    (OUT / "gene_sets_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
