#!/usr/bin/env python3
"""Journal figure style and reusable plot components.

Matches the sibling OSD-561/562 project: Arial, embedded TrueType fonts in the
PDFs, no title clutter, and every figure written as both PDF (print) and PNG
(preview).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import common as C

RC = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.labelsize": 8.5,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "legend.fontsize": 7.5,
    "legend.frameon": False,
    "figure.dpi": 120,
    "savefig.dpi": 400,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "lines.linewidth": 1.0,
    "grid.linewidth": 0.5,
    "grid.alpha": 0.3,
}

FACTOR_COLORS = {
    "Spaceflight": "#C44E52",
    "Treatment": "#4C72B0",
    "Interaction": "#8172B3",
}
REGION_COLORS = {
    "CA1": "#4C72B0",
    "DG": "#55A868",
    "FCtx": "#C44E52",
    "Ctx": "#DD8452",
}
NS_COLOR = "#BFC4CB"
UP_COLOR = "#C0392B"
DOWN_COLOR = "#2471A3"


def use_style() -> None:
    plt.rcParams.update(RC)


def save(fig, name: str, directory: Path | None = None) -> None:
    """Write a figure as PDF + PNG under results/figures (or `directory`)."""
    d = directory or C.FIG_DIR
    d.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        path = d / f"{name}.{ext}"
        fig.savefig(path)
    print(f"  [fig]  {(d / name).relative_to(C.ROOT)}.{{pdf,png}}")
    plt.close(fig)


def panel_label(ax, letter: str, dx: float = -0.13, dy: float = 1.06) -> None:
    ax.text(
        dx,
        dy,
        letter,
        transform=ax.transAxes,
        fontsize=11,
        fontweight="bold",
        va="top",
        ha="left",
    )


def volcano(
    ax,
    lfc: np.ndarray,
    pval: np.ndarray,
    *,
    genes=None,
    padj: np.ndarray | None = None,
    lfc_cut: float = C.LFC_CUTOFF,
    fdr_cut: float = C.FDR_CUTOFF,
    p_cut: float = 0.01,
    title: str = "",
    label_n: int = 8,
    annotate: bool = True,
):
    """Volcano on nominal p (y) with significance shading.

    Points are coloured by the *nominal* exploratory criterion (p < `p_cut` and
    |log2FC| >= `lfc_cut`) because, as documented in the README, no target in
    this experiment reaches FDR significance.  Targets that do pass FDR (if any)
    are ringed in black so the distinction is never lost.
    """
    lfc = np.asarray(lfc, dtype=float)
    pval = np.asarray(pval, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        y = -np.log10(pval)
    finite = np.isfinite(lfc) & np.isfinite(y)

    sig = finite & (pval < p_cut) & (np.abs(lfc) >= lfc_cut)
    up = sig & (lfc > 0)
    down = sig & (lfc < 0)
    rest = finite & ~sig

    ax.scatter(lfc[rest], y[rest], s=3, c=NS_COLOR, lw=0, alpha=0.55, rasterized=True)
    ax.scatter(lfc[down], y[down], s=6, c=DOWN_COLOR, lw=0, alpha=0.85, rasterized=True)
    ax.scatter(lfc[up], y[up], s=6, c=UP_COLOR, lw=0, alpha=0.85, rasterized=True)

    if padj is not None:
        padj = np.asarray(padj, dtype=float)
        fdr_ok = finite & np.isfinite(padj) & (padj < fdr_cut)
        if fdr_ok.any():
            ax.scatter(
                lfc[fdr_ok],
                y[fdr_ok],
                s=26,
                facecolors="none",
                edgecolors="black",
                linewidths=0.7,
                zorder=5,
            )

    ax.axvline(0, color="#909090", lw=0.5, ls=":")
    for x in (-lfc_cut, lfc_cut):
        ax.axvline(x, color="#909090", lw=0.6, ls="--")
    ax.axhline(-np.log10(p_cut), color="#909090", lw=0.6, ls="--")

    ax.set_xlabel("log$_2$ fold change")
    ax.set_ylabel("$-$log$_{10}$ nominal $P$")
    if title:
        ax.set_title(title, pad=4)

    # Full range plus a margin, so no target is clipped out of the panel.
    lim = float(np.nanmax(np.abs(lfc[finite]))) if finite.any() else 1.0
    lim = max(lim, lfc_cut * 1.4) * 1.06
    ax.set_xlim(-lim, lim)

    if annotate and genes is not None and label_n > 0 and sig.any():
        score = np.where(sig, y * np.abs(lfc), -np.inf)
        order = np.argsort(score)[::-1][:label_n]
        texts = []
        for i in order:
            if not np.isfinite(score[i]) or score[i] == -np.inf:
                continue
            texts.append(
                ax.text(
                    lfc[i],
                    y[i],
                    str(genes[i]),
                    fontsize=6,
                    style="italic",
                    ha="center",
                    va="bottom",
                )
            )
        if texts:
            try:
                from adjustText import adjust_text

                adjust_text(
                    texts,
                    ax=ax,
                    arrowprops=dict(arrowstyle="-", color="#666666", lw=0.4),
                    expand=(1.15, 1.3),
                )
            except Exception:
                pass

    # Bottom centre sits inside the empty notch of the volcano, so the count
    # never lands on top of a point or a gene label.
    n_txt = f"{int(sig.sum()):,} at $P$<{p_cut:g}, |log$_2$FC|$\\geq${lfc_cut:g}"
    ax.text(
        0.5,
        0.015,
        n_txt,
        transform=ax.transAxes,
        fontsize=6.5,
        va="bottom",
        ha="center",
        color="#444444",
    )
    return ax
