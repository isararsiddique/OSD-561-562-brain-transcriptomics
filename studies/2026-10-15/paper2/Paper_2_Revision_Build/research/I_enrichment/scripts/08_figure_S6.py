#!/usr/bin/env python3
"""Figure S6 -- background-corrected enrichment of the Paper 2 junction gene sets.

Three panels, exactly 183 mm wide, Arial throughout at >= 7 pt, Okabe-Ito fills:
 (a) observed vs expected overlap with cerebellar cell-type markers;
 (b) observed odds ratio against the minimum odds ratio the test could have
     detected, which is what makes the null interpretable;
 (c) strongest functional terms: fold enrichment vs evidence, with the BH and
     uncorrected significance frontiers drawn.

Typography notes: no mathtext is used anywhere. Matplotlib renders mathtext
sub/superscripts at 0.7x the surrounding size and, under the default fontset,
in DejaVu rather than Arial -- so a "$10^{-5}$" tick label would have been
4.9 pt DejaVu in an Arial 7 pt figure. Panel (c) therefore plots -log10(p) on a
linear axis with plain integer ticks, and every label uses literal Unicode.
The canvas is also saved WITHOUT a tight bounding box, which would otherwise
expand the width past 183 mm.
"""
from __future__ import annotations
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator

HERE = Path(__file__).resolve().parents[1]
TAB = HERE / "tables"
FIG = HERE / "figures"
FIG.mkdir(parents=True, exist_ok=True)

MM = 1 / 25.4
LE = "\u2264"   # less-than-or-equal
MI = "\u2212"   # true minus sign
MD = "\u00b7"   # middle dot

OKABE = {
    "FLT_GC_BY": "#0072B2",      # blue
    "FLT_GC_top100": "#56B4E9",  # sky blue
    "GC_VIV_top100": "#E69F00",  # orange
}
LABEL = {
    "FLT_GC_BY": f"FLT{MI}GC, BY q {LE} 0.05 (n = 47)",
    "FLT_GC_top100": f"FLT{MI}GC, strongest 106",
    "GC_VIV_top100": f"GC{MI}VIV, strongest 100",
}
GREY = "#666666"
GREEN = "#009E73"
VERM = "#D55E00"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica"],
    "font.size": 7,
    "axes.labelsize": 7.5,
    "axes.titlesize": 9,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.2,
    "ytick.major.size": 2.2,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

SHORT = {
    "negative regulation of cell communication": "neg. reg. of cell communication",
    "translational initiation": "translational initiation",
}


def main() -> None:
    ct = pd.read_csv(TAB / "celltype_enrichment_results.csv")
    ct = ct[ct.marker_definition == "all_markers"].copy()
    mde = pd.read_csv(TAB / "celltype_minimum_detectable_effect.csv")
    mde = mde[mde.marker_definition == "all_markers"]
    gp = pd.read_csv(TAB / "gprofiler_recomputed_with_BH.csv.gz")

    order = (
        ct.groupby("cell_class").n_markers_in_background.first()
        .sort_values(ascending=False).index.tolist()
    )
    ypos = {c: i for i, c in enumerate(order)}

    fig = plt.figure(figsize=(183 * MM, 104 * MM))
    gs = fig.add_gridspec(
        1, 3, width_ratios=[1.10, 1.0, 1.22], wspace=0.44,
        left=0.128, right=0.962, top=0.86, bottom=0.145,
    )
    ax_a, ax_b, ax_c = (fig.add_subplot(gs[0, i]) for i in range(3))

    # ---- (a) observed vs expected overlap -------------------------------
    for cls in order:
        ax_a.axhline(ypos[cls], color="#EEEEEE", lw=0.4, zorder=0)
    for i, (gsx, sub) in enumerate(ct.groupby("gene_set")):
        off = (i - 1) * 0.24
        sub = sub.set_index("cell_class").reindex(order)
        ys = np.array([ypos[c] for c in order]) + off
        ax_a.hlines(ys, sub.expected_overlap, sub.overlap,
                    color=OKABE[gsx], lw=0.8, alpha=0.55, zorder=2)
        ax_a.scatter(sub.expected_overlap, ys, s=7, facecolor="white",
                     edgecolor=OKABE[gsx], linewidth=0.6, zorder=3)
        ax_a.scatter(sub.overlap, ys, s=11, color=OKABE[gsx],
                     edgecolor="none", zorder=4)
    ax_a.set_yticks(range(len(order)))
    ax_a.set_yticklabels(order)
    ax_a.invert_yaxis()
    ax_a.set_xlabel("Genes shared with cell-type marker set")
    ax_a.set_xlim(-0.4, 10.6)
    ax_a.set_xticks([0, 2, 4, 6, 8, 10])
    ax_a.spines[["top", "right"]].set_visible(False)
    ax_a.set_title("a)  Cell-type markers", loc="left", fontweight="bold", pad=5)
    ax_a.text(0.985, 0.035,
              "open = expected under the\nmeasured background\nfilled = observed\n\n"
              f"no class passes BH q {LE} 0.05",
              transform=ax_a.transAxes, ha="right", va="bottom", fontsize=7,
              color=GREY, linespacing=1.3)

    # ---- (b) observed OR vs minimum detectable OR ------------------------
    m = mde.merge(
        ct[["gene_set", "cell_class", "odds_ratio_haldane", "q_BH_within_set"]],
        on=["gene_set", "cell_class"], how="left",
    )
    lo, hi = 0.3, 160.0
    ax_b.fill_between([lo, hi], [lo, hi], [hi, hi], color=GREEN, alpha=0.09, zorder=1)
    ax_b.plot([lo, hi], [lo, hi], color=GREY, lw=0.7, ls="--", zorder=2)
    for gsx, sub in m.groupby("gene_set"):
        ax_b.scatter(sub.min_detectable_odds_ratio_p05, sub.odds_ratio_haldane,
                     s=12, color=OKABE[gsx], edgecolor="none", alpha=0.85, zorder=3)
    ax_b.set_xscale("log"); ax_b.set_yscale("log")
    ax_b.set_xlim(lo, hi); ax_b.set_ylim(lo, hi)
    for ax in (ax_b.xaxis, ax_b.yaxis):
        ax.set_major_locator(FixedLocator([0.5, 1, 2, 5, 10, 20, 50, 100]))
        ax.set_major_formatter(FixedFormatter(["0.5", "1", "2", "5", "10", "20", "50", "100"]))
        ax.set_minor_locator(NullLocator())
    ax_b.set_xlabel("Minimum detectable odds ratio")
    ax_b.set_ylabel("Observed odds ratio")
    ax_b.spines[["top", "right"]].set_visible(False)
    ax_b.set_title("b)  Detectability", loc="left", fontweight="bold", pad=5)
    ax_b.text(0.045, 0.55, f"detectable at p {LE} 0.05", transform=ax_b.transAxes,
              ha="left", va="center", fontsize=7, color="#00705A")
    ax_b.text(0.975, 0.045, "below the floor\n(all 51 class tests)",
              transform=ax_b.transAxes, ha="right", va="bottom", fontsize=7,
              color=GREY, linespacing=1.3)
    fig.legend(
        handles=[Line2D([], [], marker="o", ls="none", ms=3.2, color=OKABE[k],
                        label=LABEL[k]) for k in OKABE],
        loc="upper center", bbox_to_anchor=(0.5, 1.0), ncol=3, frameon=False,
        handletextpad=0.4, borderpad=0.1, columnspacing=1.6,
    )

    # ---- (c) functional terms -------------------------------------------
    gp = gp.assign(evidence=-np.log10(gp.p_hypergeom_uncorrected))
    fam = int(gp.groupby("gene_set").n_terms_in_BH_family.first().max())
    bh_frontier = -np.log10(0.05 / fam)
    for gsx, sub in gp.groupby("gene_set"):
        s = sub.nlargest(300, "evidence")
        ax_c.scatter(s.fold_enrichment.clip(upper=55), s.evidence,
                     s=8, color=OKABE[gsx], alpha=0.55, edgecolor="none", zorder=3)
    ax_c.axhline(bh_frontier, color=VERM, lw=0.9, zorder=4)
    ax_c.axhline(-np.log10(0.05), color=GREY, lw=0.6, ls=":", zorder=4)
    ax_c.set_xscale("log")
    ax_c.set_xlim(0.85, 80)
    ax_c.xaxis.set_major_locator(FixedLocator([1, 2, 5, 10, 20, 50]))
    ax_c.xaxis.set_major_formatter(FixedFormatter(["1", "2", "5", "10", "20", "50"]))
    ax_c.xaxis.set_minor_locator(NullLocator())
    ax_c.set_ylim(0.6, bh_frontier + 2.35)
    ax_c.set_yticks([1, 2, 3, 4, 5])
    ax_c.set_xlabel("Fold enrichment over background")
    ax_c.set_ylabel(f"Evidence, {MI}log10 uncorrected p")
    ax_c.spines[["top", "right"]].set_visible(False)
    ax_c.set_title("c)  Functional terms", loc="left", fontweight="bold", pad=5)
    ax_c.text(1.0, bh_frontier - 0.16,
              f"BH q {LE} 0.05 frontier (p {LE} 0.000019)",
              ha="left", va="top", fontsize=7, color=VERM)
    ax_c.text(0.025, 0.737, "dotted line: uncorrected p = 0.05",
              transform=ax_c.transAxes, ha="left", va="top", fontsize=7, color=GREY)

    # Strongest term per set, named above the data. Effect sizes and p values
    # stay in the caption and the CSV: at 7 pt a name plus OR plus p is wider
    # than the panel, and shrinking the type below 7 pt is not an option.
    best = gp.loc[gp.groupby("gene_set").evidence.idxmax()]
    for j, k in enumerate(OKABE):
        r = best[best.gene_set == k].iloc[0]
        nm = SHORT.get(r["name"], r["name"])
        if j == 0:
            ax_c.text(0.025, 0.985, "strongest term per set:",
                      transform=ax_c.transAxes, ha="left", va="top",
                      fontsize=7, color=GREY)
        ax_c.text(0.025, 0.923 - j * 0.062, nm, transform=ax_c.transAxes,
                  ha="left", va="top", fontsize=7, color=OKABE[k])
    # Term count and source list live in the caption: at 7 pt there is no
    # room for them inside the panel without covering data.

    fig.savefig(FIG / "Figure_S6_enrichment.pdf", bbox_inches=None)
    fig.savefig(FIG / "Figure_S6_enrichment.png", dpi=450, bbox_inches=None)
    print("wrote", FIG / "Figure_S6_enrichment.pdf")
    print(f"class tests={len(ct)}  term tests={len(gp)}  "
          f"BH frontier -log10 p={bh_frontier:.2f}")


if __name__ == "__main__":
    main()
