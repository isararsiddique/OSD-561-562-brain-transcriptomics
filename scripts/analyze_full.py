#!/usr/bin/env python3
"""Full-data analysis for OSD-561 & OSD-562.

Where analyze.py reports only the five matched spaceflight contrasts, this
script uses EVERY comparison in the GeneLab DE tables and every downloaded
count table:

  * all unique group comparisons (reverse-direction duplicates collapsed),
    categorised as Spaceflight / Age / Environment / Confounded, with DEG
    counts for each  ->  DEG-landscape figures (one per region)
  * volcano plots for every *clean* single-factor contrast (Spaceflight, Age,
    Environment) -- not just the spaceflight ones
  * a focus-gene log2FC matrix spanning all clean contrasts (per region)
  * QC panels built from the RSEM (raw) and Normalized count tables so that
    the two otherwise-unused files are put to work:
        - library size per sample            (RSEM)
        - genes detected per sample          (RSEM)
        - VST vs Normalized concordance      (VST + Normalized)

Individual figures land in results/figures/ ; make_panels.py composes them.

Run:
  ./.venv/bin/python scripts/analyze_full.py
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch

import analyze as A

ROOT = A.ROOT
FIG = A.FIG
TAB = A.TAB

CAT_ORDER = {"Spaceflight": 0, "Age": 1, "Environment": 2, "Confounded": 3}
# per-category volcano wording (up_label, down_label, xlabel)
VOLC_WORDS = {
    "Spaceflight": ("Up in flight", "Down in flight",
                    "log$_2$ fold change  (Flight / Ground)"),
    "Age": ("Up in old (29 wk)", "Down in old (29 wk)",
            "log$_2$ fold change  (Old / Young)"),
    "Environment": ("Up on ISS", "Down on ISS",
                    "log$_2$ fold change  (On ISS / On Earth)"),
}


# --------------------------------------------------------------------------- #
# 1. Contrast table (all unique comparisons) + DEG counts
# --------------------------------------------------------------------------- #
def build_contrast_table(study: str, dge: pd.DataFrame):
    label = A.STUDIES[study]["label"]
    contrasts = A.all_unique_contrasts(dge)
    rows = []
    for c in contrasts:
        stats = A.contrast_stats_oriented(dge, c["num"], c["den"])
        counts = A.deg_counts(stats)
        rows.append({
            "study": study, "region": label, "category": c["category"],
            "label": c["label"], "factors_differ": c["factors_differ"],
            "num_group": c["num"], "den_group": c["den"],
            **counts,
        })
    df = pd.DataFrame(rows)
    df["_o"] = df["category"].map(CAT_ORDER)
    df = df.sort_values(["_o", "total", "label"],
                        ascending=[True, False, True]).drop(columns="_o")
    return df.reset_index(drop=True), contrasts


# --------------------------------------------------------------------------- #
# 2. DEG landscape (diverging bars, every unique contrast)
# --------------------------------------------------------------------------- #
def plot_deg_landscape(df: pd.DataFrame, region_label: str, out_stem: str):
    n = len(df)
    y = np.arange(n)[::-1]
    fig_h = 1.5 + 0.34 * n
    fig, ax = plt.subplots(figsize=(7.8, fig_h))

    ax.barh(y, -df["down"].values, color=A.PAL_TRT["Ground Control"],
            height=0.72, edgecolor="white", linewidth=0.4, zorder=3)
    ax.barh(y, df["up"].values, color=A.PAL_TRT["Space Flight"],
            height=0.72, edgecolor="white", linewidth=0.4, zorder=3)
    ax.axvline(0, color="0.25", lw=0.9, zorder=4)

    # symlog (linear near 0, log beyond) keeps the tiny spaceflight bars
    # readable alongside the huge confounded/environment effects
    xmax = int(max(df["up"].max(), df["down"].max(), 1))
    ax.set_xscale("symlog", linthresh=5, linscale=0.9)
    ax.set_xlim(-(xmax * 2.2 + 2), xmax * 2.2 + 2)
    ax.xaxis.set_major_formatter(plt.FuncFormatter(
        lambda v, _pos: f"{abs(int(v))}"))

    ax.set_yticks(y)
    ax.set_yticklabels(df["label"], fontsize=7.5)
    for tick, cat in zip(ax.get_yticklabels(), df["category"]):
        tick.set_color(A.PAL_CAT[cat])
        tick.set_fontweight("bold" if cat != "Confounded" else "normal")

    for yi, u, d in zip(y, df["up"].values, df["down"].values):
        if u > 0:
            ax.annotate(str(int(u)), xy=(u, yi), xytext=(3, 0),
                        textcoords="offset points", va="center", ha="left",
                        fontsize=6.5, color=A.PAL_TRT["Space Flight"])
        if d > 0:
            ax.annotate(str(int(d)), xy=(-d, yi), xytext=(-3, 0),
                        textcoords="offset points", va="center", ha="right",
                        fontsize=6.5, color=A.PAL_TRT["Ground Control"])

    ax.set_xlabel("DEGs  (FDR < 0.05, |log$_2$FC| $\\geq$ 1)   "
                  "$\\blacktriangleleft$ down  |  up $\\blacktriangleright$  in first-named group"
                  "   (symlog scale)")
    ax.set_title(f"{region_label} — differential expression across all "
                 f"{n} unique contrasts", pad=8)

    dir_leg = [Patch(facecolor=A.PAL_TRT["Space Flight"], label="Up"),
               Patch(facecolor=A.PAL_TRT["Ground Control"], label="Down")]
    cat_leg = [Patch(facecolor=A.PAL_CAT[c], label=c)
               for c in CAT_ORDER if (df["category"] == c).any()]
    leg1 = ax.legend(handles=dir_leg, title="Direction", loc="lower right",
                     frameon=False, fontsize=7, title_fontsize=7.5)
    ax.add_artist(leg1)
    ax.legend(handles=cat_leg, title="Contrast type (label colour)",
              loc="upper right", frameon=False, fontsize=7, title_fontsize=7.5)

    sns.despine(ax=ax, left=True)
    ax.tick_params(axis="y", length=0)
    A._save(fig, out_stem)


# --------------------------------------------------------------------------- #
# 3. Volcanoes for every clean single-factor contrast
# --------------------------------------------------------------------------- #
def clean_volcanoes(study: str, dge: pd.DataFrame, contrasts: list[dict],
                    focus: list[str]):
    label = A.STUDIES[study]["label"]
    made = []
    for c in contrasts:
        if c["category"] == "Confounded":
            continue
        stats = A.contrast_stats_oriented(dge, c["num"], c["den"])
        up_l, down_l, xl = VOLC_WORDS[c["category"]]
        stem = f"ALL_volcano_{study}_{c['tag']}"
        A.plot_volcano(
            stats, set(focus),
            title=f"{study} {label}\n{c['label']}",
            out_stem=stem, up_label=up_l, down_label=down_l, xlabel=xl,
        )
        made.append({**c, "stem": stem})
    return made


# --------------------------------------------------------------------------- #
# 4. Focus-gene log2FC across all clean contrasts (per region)
# --------------------------------------------------------------------------- #
def focus_matrix_all(study: str, dge: pd.DataFrame, contrasts: list[dict],
                     focus: list[str]):
    label = A.STUDIES[study]["label"]
    focus_up = {g.upper() for g in focus}
    cols = {}
    for c in contrasts:
        if c["category"] == "Confounded":
            continue
        stats = A.contrast_stats_oriented(dge, c["num"], c["den"])
        fsub = stats[stats["SYMBOL"].str.upper().isin(focus_up)]
        cols[c["label"]] = fsub.set_index("SYMBOL")["log2FC"].groupby(level=0).first()
    if not cols:
        return
    mat = pd.DataFrame(cols)
    # keep the focus-list ordering for rows
    order = [g for g in focus if g in mat.index]
    mat = mat.loc[order]
    A.plot_lfc_matrix(
        mat,
        title=f"{study} {label} — focus-gene log$_2$FC across all clean contrasts",
        out_stem=f"ALL_focus_log2FC_{study}",
    )


# --------------------------------------------------------------------------- #
# 5. QC from RSEM (raw) + Normalized counts  (the two unused files)
# --------------------------------------------------------------------------- #
def _libsize_frame():
    frames = []
    for st in A.STUDIES:
        rsem = A.load_counts(st, "RSEM_Unnormalized_Counts")
        meta = A.sample_metadata(rsem.columns)
        f = pd.DataFrame({
            "sample": rsem.columns,
            "library_size_M": rsem.sum(axis=0).values / 1e6,
            "genes_detected": (rsem > 0).sum(axis=0).values,
        }).set_index("sample").join(meta)
        f["region"] = A.STUDIES[st]["label"]
        frames.append(f)
    return pd.concat(frames)


def plot_qc_metric(lib: pd.DataFrame, ycol: str, ylabel: str, title: str,
                   out_stem: str):
    fig, ax = plt.subplots(figsize=(4.6, 3.6))
    order = ["Cerebellum", "Hippocampus"]
    sns.boxplot(data=lib, x="region", y=ycol, order=order, hue="treatment_full",
                hue_order=["Ground Control", "Space Flight"],
                palette=A.PAL_TRT, width=0.6, fliersize=0, linewidth=0.8, ax=ax)
    sns.stripplot(data=lib, x="region", y=ycol, order=order,
                  hue="treatment_full",
                  hue_order=["Ground Control", "Space Flight"],
                  palette=A.PAL_TRT, dodge=True, size=3.2, linewidth=0.4,
                  edgecolor="white", ax=ax, legend=False)
    ax.set_xlabel("")
    ax.set_ylabel(ylabel)
    ax.set_title(title, pad=8)
    h, l = ax.get_legend_handles_labels()
    ax.legend(h[:2], l[:2], title="Treatment", frameon=False,
              loc="best", fontsize=7, title_fontsize=7.5)
    sns.despine(ax=ax)
    A._save(fig, out_stem)


def plot_qc_concordance(out_stem: str):
    fig, ax = plt.subplots(figsize=(4.6, 3.8))
    rng = np.random.default_rng(0)
    colors = {"Cerebellum": A.PAL_REG["Cb"], "Hippocampus": A.PAL_REG["HPC"]}
    for st in A.STUDIES:
        label = A.STUDIES[st]["label"]
        vst = A.load_counts(st, "VST_Counts")
        norm = A.load_counts(st, "Normalized_Counts")
        genes = vst.index.intersection(norm.index)
        samples = vst.columns.intersection(norm.columns)
        v = vst.loc[genes, samples].values.ravel()
        nlog = np.log2(norm.loc[genes, samples].values.ravel() + 1.0)
        # drop structural zeros so the correlation reflects expressed genes
        m = (v > 0) | (nlog > 0)
        v, nlog = v[m], nlog[m]
        r = np.corrcoef(v, nlog)[0, 1]
        idx = rng.choice(v.size, size=min(8000, v.size), replace=False)
        ax.scatter(nlog[idx], v[idx], s=3, alpha=0.25, linewidths=0,
                   color=colors[label], rasterized=True,
                   label=f"{label} (r = {r:.3f})")
    ax.set_xlabel("log$_2$(Normalized counts + 1)")
    ax.set_ylabel("VST counts")
    ax.set_title("VST vs Normalized concordance", pad=8)
    leg = ax.legend(frameon=False, fontsize=7, markerscale=3,
                    loc="upper left")
    for lh in leg.legend_handles:
        lh.set_alpha(1)
    sns.despine(ax=ax)
    A._save(fig, out_stem)


# --------------------------------------------------------------------------- #
def main() -> int:
    A.set_journal_style()
    FIG.mkdir(parents=True, exist_ok=True)
    focus = A.read_focus_genes()
    print(f"Focus genes ({len(focus)}): {', '.join(focus)}")

    all_tables = []
    for st in A.STUDIES:
        info = A.STUDIES[st]
        print(f"\n=== {st} ({info['label']}) : all contrasts ===")
        dge = A.load_dge(st)
        table, contrasts = build_contrast_table(st, dge)
        all_tables.append(table)
        n_clean = int((table["category"] != "Confounded").sum())
        print(f"  {len(table)} unique contrasts "
              f"({n_clean} clean single-factor, "
              f"{len(table) - n_clean} confounded)")
        print(table[["category", "label", "up", "down", "total"]]
              .to_string(index=False))

        plot_deg_landscape(table, info["label"], f"ALL_deg_landscape_{st}")
        made = clean_volcanoes(st, dge, contrasts, focus)
        print(f"  clean-contrast volcanoes: {len(made)}")
        focus_matrix_all(st, dge, contrasts, focus)

    # QC using RSEM (raw) + Normalized counts
    print("\n=== QC (RSEM raw + Normalized counts) ===")
    lib = _libsize_frame()
    print(lib.groupby("region")[["library_size_M", "genes_detected"]]
          .agg(["mean", "min", "max"]).round(1).to_string())
    plot_qc_metric(lib, "library_size_M", "Library size (million reads)",
                   "Sequencing depth per sample", "QC_library_size")
    plot_qc_metric(lib, "genes_detected", "Genes detected (count > 0)",
                   "Genes detected per sample", "QC_genes_detected")
    plot_qc_concordance("QC_concordance")

    # persist the all-contrast DEG table for export/inspection
    full = pd.concat(all_tables, ignore_index=True)
    TAB.mkdir(parents=True, exist_ok=True)
    out = TAB / "all_contrasts_DEG_summary.csv"
    full.to_csv(out, index=False)
    print(f"\n[table] {out.relative_to(ROOT)}  ({len(full)} contrasts)")
    print("\nAll full-data figures written to results/figures/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
