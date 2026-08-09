#!/usr/bin/env python3
"""Compose the individual result figures into labelled multi-panel figures.

Full-data panel set (every downloaded file contributes):
  Figure1_PCA                  global structure (combined + per region)
  Figure2_DEG_landscape        DEGs across ALL unique contrasts, per region
  Figure3_Volcanoes_Spaceflight   Flight vs Ground (matched), all 5
  Figure4_Volcanoes_Age           Old vs Young (matched), all 5
  Figure5_Volcanoes_Environment   On ISS vs On Earth (matched), all 4
  Figure6_FocusGene_log2FC     focus genes across all clean contrasts
  Figure7_FocusGene_heatmaps   focus-gene VST z-scores (combined + regions)
  Figure8_QC                   library size / genes detected / VST-Norm concordance

Reads the 400-dpi PNGs from results/figures/ and writes PDF (for print) + PNG
to results/figures/panels/.

Run:
  ./.venv/bin/python scripts/make_panels.py
"""
from __future__ import annotations

import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
FIG = ROOT / "results" / "figures"
OUT = FIG / "panels"
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({"font.family": "sans-serif",
                     "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"]})


def place(ax, name: str, letter: str | None):
    """Draw a PNG into an axes with an optional bold panel letter."""
    img = FIG / f"{name}.png"
    if not img.exists():
        ax.text(0.5, 0.5, f"[missing]\n{name}", ha="center", va="center",
                fontsize=8, color="crimson")
        ax.axis("off")
        return
    ax.imshow(mpimg.imread(img))
    ax.axis("off")
    if letter:
        ax.annotate(letter, xy=(0, 1), xytext=(2, 2),
                    xycoords="axes fraction", textcoords="offset points",
                    fontsize=17, fontweight="bold", va="bottom", ha="left")


def save(fig, stem: str):
    fig.savefig(OUT / f"{stem}.pdf", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  [panel] {stem}.pdf / .png")


def grid_panel(stems: list[str], ncols: int, out_stem: str,
               cell=(4.6, 4.3)):
    """Arrange a list of figure stems into an auto-sized grid with A,B,C… ."""
    n = len(stems)
    if n == 0:
        print(f"  [skip] {out_stem}: no source figures found")
        return
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(ncols * cell[0], nrows * cell[1]),
                             constrained_layout=True)
    axes = [axes] if n == 1 and not hasattr(axes, "flat") else list(axes.flat)
    for i, ax in enumerate(axes):
        if i < n:
            place(ax, stems[i], chr(65 + i))
        else:
            ax.axis("off")
    save(fig, out_stem)


def globstems(pattern: str) -> list[str]:
    """Sorted figure stems matching a glob (Cb OSD-561 sorts before HPC)."""
    return [p.stem for p in sorted(FIG.glob(pattern))]


# --------------------------------------------------------------------------- #
# Figure 1 — PCA (2x2)
# --------------------------------------------------------------------------- #
fig, ax = plt.subplots(2, 2, figsize=(9.0, 7.4), constrained_layout=True)
place(ax[0, 0], "COMBINED_PCA_by_region", "A")
place(ax[0, 1], "COMBINED_PCA_by_treatment", "B")
place(ax[1, 0], "OSD-561_PCA", "C")
place(ax[1, 1], "OSD-562_PCA", "D")
save(fig, "Figure1_PCA")

# --------------------------------------------------------------------------- #
# Figure 2 — DEG landscape across ALL unique contrasts (Cb | HPC)
# --------------------------------------------------------------------------- #
fig, ax = plt.subplots(1, 2, figsize=(15.5, 8.6), constrained_layout=True)
place(ax[0], "ALL_deg_landscape_OSD-561", "A")
place(ax[1], "ALL_deg_landscape_OSD-562", "B")
save(fig, "Figure2_DEG_landscape")

# --------------------------------------------------------------------------- #
# Figures 3-5 — volcanoes for every clean single-factor contrast
# --------------------------------------------------------------------------- #
grid_panel(globstems("ALL_volcano_*_Spaceflight_*.png"), ncols=3,
           out_stem="Figure3_Volcanoes_Spaceflight")
grid_panel(globstems("ALL_volcano_*_Age_*.png"), ncols=3,
           out_stem="Figure4_Volcanoes_Age")
grid_panel(globstems("ALL_volcano_*_Environment_*.png"), ncols=2,
           out_stem="Figure5_Volcanoes_Environment")

# --------------------------------------------------------------------------- #
# Figure 6 — focus-gene log2FC across all clean contrasts (Cb | HPC)
# --------------------------------------------------------------------------- #
fig, ax = plt.subplots(1, 2, figsize=(13.5, 5.6), constrained_layout=True)
place(ax[0], "ALL_focus_log2FC_OSD-561", "A")
place(ax[1], "ALL_focus_log2FC_OSD-562", "B")
save(fig, "Figure6_FocusGene_log2FC")

# --------------------------------------------------------------------------- #
# Figure 7 — focus-gene heat maps (combined full width + Cb/HPC)
# --------------------------------------------------------------------------- #
fig = plt.figure(figsize=(11.0, 9.2), constrained_layout=True)
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.1])
place(fig.add_subplot(gs[0, :]), "COMBINED_focus_heatmap", "A")
place(fig.add_subplot(gs[1, 0]), "OSD-561_focus_heatmap", "B")
place(fig.add_subplot(gs[1, 1]), "OSD-562_focus_heatmap", "C")
save(fig, "Figure7_FocusGene_heatmaps")

# --------------------------------------------------------------------------- #
# Figure 8 — QC (library size / genes detected / VST-Norm concordance)
# --------------------------------------------------------------------------- #
grid_panel(["QC_library_size", "QC_genes_detected", "QC_concordance"],
           ncols=3, out_stem="Figure8_QC", cell=(4.7, 3.9))

print("Done -> results/figures/panels/")
