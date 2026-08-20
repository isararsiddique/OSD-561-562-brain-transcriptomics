#!/usr/bin/env python3
"""Compose the numbered main figures from the exported result tables.

Reads only from results/tables/csv/, so it can be rerun without repeating any
analysis, and every panel is redrawn as vector art rather than stitched from
bitmaps.

Outputs results/figures/panels/Figure1..Figure7 as PDF + PNG.

Usage:  python scripts/make_panels.py   (after the analyze_* scripts)
"""
from __future__ import annotations

import sys

import numpy as np
import pandas as pd

import common as C
import plotting as P

CSV = C.CSV_DIR

LANDSCAPE = [
    ("FLT_vs_GC_in_SAL", "Flight vs Ground (saline)"),
    ("FLT_vs_GC_in_BuOE", "Flight vs Ground (BuOE)"),
    ("BuOE_vs_SAL_in_GC", "BuOE vs saline (ground)"),
    ("BuOE_vs_SAL_in_FLT", "BuOE vs saline (flight)"),
    ("spaceflight", "Spaceflight main effect"),
    ("treatment", "Treatment main effect"),
    ("interaction", "Flight $\\times$ BuOE interaction"),
]


def need(name: str) -> pd.DataFrame:
    path = CSV / name
    if not path.exists():
        raise SystemExit(
            f"missing {path.relative_to(C.ROOT)} -- run the analyze_* scripts first"
        )
    return pd.read_csv(path)


# --------------------------------------------------------------------------- #
def figure1(qc: pd.DataFrame, cohort: pd.DataFrame) -> None:
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(11.0, 6.4))
    gs = fig.add_gridspec(2, 3, hspace=0.42, wspace=0.32)

    # (a) design schematic
    ax = fig.add_subplot(gs[0, 0])
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")
    cells = [
        (1.1, 5.6, "GC_SAL"), (5.4, 5.6, "GC_BuOE"),
        (1.1, 1.7, "FLT_SAL"), (5.4, 1.7, "FLT_BuOE"),
    ]
    for x, y, grp in cells:
        ax.add_patch(plt.Rectangle((x, y), 3.5, 3.0, facecolor=C.GROUP_COLORS[grp],
                                   alpha=0.22, edgecolor=C.GROUP_COLORS[grp], lw=1.1))
        ax.text(x + 1.75, y + 2.05, C.GROUP_LABELS[grp].replace(" / ", "\n"),
                ha="center", va="center", fontsize=7.2, color="#222222")
        ax.text(x + 1.75, y + 0.72, "n = 3 ROIs\n$\\times$ 4 regions", ha="center",
                va="center", fontsize=6.2, color="#555555")
    ax.text(4.85, 9.35, "Treatment", ha="center", fontsize=7.8, fontweight="bold")
    ax.text(0.35, 5.2, "Spaceflight", va="center", rotation=90, fontsize=7.8,
            fontweight="bold")
    ax.text(4.85, 0.55, "48 ROIs total  |  15,782 targets", ha="center", fontsize=6.6,
            color="#444444")
    P.panel_label(ax, "a", dx=-0.02, dy=1.04)
    ax.set_title("2 $\\times$ 2 factorial design", pad=2)

    # (b) regions / accessions
    ax = fig.add_subplot(gs[0, 1])
    ax.axis("off")
    tab = cohort[cohort["region_short"] != "combined"]
    lines = [f"{r.osd}   {r.glds}   {r.region}" for r in tab.itertuples()]
    ax.text(0.0, 0.93, "One GeoMx DSP study, four OSDR accessions",
            fontsize=7.4, fontweight="bold", va="top", transform=ax.transAxes)
    for i, line in enumerate(lines):
        ax.text(0.0, 0.78 - i * 0.115, line, fontsize=6.9, va="top",
                family="monospace", transform=ax.transAxes)
    ax.text(0.0, 0.30,
            "Processed layer from GEO GSE239336.\n"
            "OSDR hosts raw FASTQ only for these\n"
            "four studies (no GeneLab count or DE\n"
            "tables), so counts come from GEO.",
            fontsize=6.4, va="top", color="#444444", transform=ax.transAxes)
    P.panel_label(ax, "b", dx=-0.02, dy=1.04)
    ax.set_title("Datasets", pad=2)

    # (c) sequencing depth
    ax = fig.add_subplot(gs[0, 2])
    for gi, g in enumerate(C.GROUPS):
        xs, ys = [], []
        for ri, region in enumerate(C.REGION_ORDER):
            sel = qc[(qc["region_short"] == region) & (qc["group"] == g)]
            xs.extend(np.full(len(sel), ri + (gi - 1.5) * 0.17))
            ys.extend(sel["raw_reads"] / 1e6)
        ax.scatter(xs, ys, s=13, color=C.GROUP_COLORS[g], lw=0.3, edgecolor="white",
                   label=C.GROUP_LABELS[g])
    ax.set_xticks(range(len(C.REGION_ORDER)))
    ax.set_xticklabels(C.REGION_ORDER)
    ax.set_ylabel("Raw reads (millions)")
    ax.legend(fontsize=5.6, ncols=2, columnspacing=0.7, handletextpad=0.25)
    ax.grid(axis="y", ls=":", lw=0.5, alpha=0.4)
    P.panel_label(ax, "c")
    ax.set_title("Sequencing depth per ROI", pad=2)

    # (d) saturation
    ax = fig.add_subplot(gs[1, 0])
    for gi, g in enumerate(C.GROUPS):
        xs, ys = [], []
        for ri, region in enumerate(C.REGION_ORDER):
            sel = qc[(qc["region_short"] == region) & (qc["group"] == g)]
            xs.extend(np.full(len(sel), ri + (gi - 1.5) * 0.17))
            ys.extend(sel["seq_saturation"])
        ax.scatter(xs, ys, s=13, color=C.GROUP_COLORS[g], lw=0.3, edgecolor="white")
    ax.set_xticks(range(len(C.REGION_ORDER)))
    ax.set_xticklabels(C.REGION_ORDER)
    ax.set_ylabel("Sequencing saturation (%)")
    ax.grid(axis="y", ls=":", lw=0.5, alpha=0.4)
    P.panel_label(ax, "d")
    ax.set_title("Library saturation", pad=2)

    # (e) background separation
    ax = fig.add_subplot(gs[1, 1])
    if "snr_median_over_negprobe_log2" in qc:
        for gi, g in enumerate(C.GROUPS):
            xs, ys = [], []
            for ri, region in enumerate(C.REGION_ORDER):
                sel = qc[(qc["region_short"] == region) & (qc["group"] == g)]
                xs.extend(np.full(len(sel), ri + (gi - 1.5) * 0.17))
                ys.extend(sel["snr_median_over_negprobe_log2"])
            ax.scatter(xs, ys, s=13, color=C.GROUP_COLORS[g], lw=0.3, edgecolor="white")
        ax.set_xticks(range(len(C.REGION_ORDER)))
        ax.set_xticklabels(C.REGION_ORDER)
        ax.set_ylabel("Median target $-$ NegProbe\n(log$_2$)")
        ax.grid(axis="y", ls=":", lw=0.5, alpha=0.4)
    P.panel_label(ax, "e")
    ax.set_title("Signal above background", pad=2)

    # (f) ROI layout for one slide
    ax = fig.add_subplot(gs[1, 2])
    for region in C.REGION_ORDER:
        s = qc[(qc["slide"] == "FLT_SAL_2") & (qc["region_short"] == region)]
        ax.scatter(s["roi_x"] / 1e3, s["roi_y"] / 1e3, s=40,
                   color=P.REGION_COLORS[region], lw=0.4, edgecolor="white",
                   label=region)
    ax.invert_yaxis()
    ax.set_xlabel("ROI x (10$^3$ px)")
    ax.set_ylabel("ROI y (10$^3$ px)")
    ax.legend(fontsize=6, title="Region", title_fontsize=6)
    P.panel_label(ax, "f")
    ax.set_title("ROI placement (slide FLT_SAL_2):\n3 sections $\\times$ 4 regions", pad=2)

    P.save(fig, "Figure1_design_and_QC", C.PANEL_DIR)


def figure2(coords: pd.DataFrame, var: pd.DataFrame, vp: pd.DataFrame) -> None:
    import matplotlib.pyplot as plt

    def pv(scope, pc):
        v = var[(var["scope"] == scope) & (var["PC"] == pc)]["variance_explained"]
        return float(v.iloc[0]) * 100 if len(v) else np.nan

    markers = {"GC_SAL": "o", "GC_BuOE": "s", "FLT_SAL": "^", "FLT_BuOE": "D"}
    fig = plt.figure(figsize=(11.0, 6.2))
    gs = fig.add_gridspec(2, 3, hspace=0.40, wspace=0.30)

    comb = coords[coords["scope"] == "combined"]
    ax = fig.add_subplot(gs[0, 0])
    for region in C.REGION_ORDER:
        s = comb[comb["region_short"] == region]
        ax.scatter(s["PC1"], s["PC2"], s=30, color=P.REGION_COLORS[region], lw=0.4,
                   edgecolor="white", label=region)
    ax.set_xlabel(f"PC1 ({pv('combined','PC1'):.1f}%)")
    ax.set_ylabel(f"PC2 ({pv('combined','PC2'):.1f}%)")
    ax.legend(title="Region", fontsize=6.3, title_fontsize=6.3)
    P.panel_label(ax, "a")
    ax.set_title("All 48 ROIs by region", pad=2)

    ax = fig.add_subplot(gs[0, 1])
    for g in C.GROUPS:
        s = comb[comb["group"] == g]
        ax.scatter(s["PC1"], s["PC2"], s=30, color=C.GROUP_COLORS[g], marker=markers[g],
                   lw=0.4, edgecolor="white", label=C.GROUP_LABELS[g])
    ax.set_xlabel(f"PC1 ({pv('combined','PC1'):.1f}%)")
    ax.set_ylabel(f"PC2 ({pv('combined','PC2'):.1f}%)")
    ax.legend(fontsize=5.8)
    P.panel_label(ax, "b")
    ax.set_title("Same ROIs by treatment group", pad=2)

    ax = fig.add_subplot(gs[0, 2])
    vp2 = vp.sort_values("median_frac_variance")
    colors = {"Region": "#7F8C8D", "Spaceflight": P.FACTOR_COLORS["Spaceflight"],
              "Treatment": P.FACTOR_COLORS["Treatment"], "Section (rep)": "#95A5A6"}
    ax.barh(vp2["factor"], vp2["median_frac_variance"] * 100,
            color=[colors.get(f, "#888") for f in vp2["factor"]], height=0.6)
    for y, v in enumerate(vp2["median_frac_variance"]):
        ax.text(v * 100 + 0.4, y, f"{v*100:.1f}%", va="center", fontsize=6.5)
    ax.set_xlabel("Median % of per-target variance")
    ax.tick_params(axis="y", labelsize=6.5)
    ax.grid(axis="x", ls=":", lw=0.5, alpha=0.4)
    P.panel_label(ax, "c", dx=-0.30)
    ax.set_title("Variance attribution", pad=2)

    for k, region in enumerate(C.REGION_ORDER):
        ax = fig.add_subplot(gs[1, k]) if k < 3 else None
        if ax is None:
            break
        s = coords[coords["scope"] == region]
        for g in C.GROUPS:
            ss = s[s["group"] == g]
            ax.scatter(ss["PC1"], ss["PC2"], s=34, color=C.GROUP_COLORS[g],
                       marker=markers[g], lw=0.4, edgecolor="white")
        ax.set_xlabel(f"PC1 ({pv(region,'PC1'):.1f}%)")
        ax.set_ylabel(f"PC2 ({pv(region,'PC2'):.1f}%)")
        P.panel_label(ax, "def"[k])
        ax.set_title(region, pad=2)
    P.save(fig, "Figure2_global_structure", C.PANEL_DIR)


def figure3(audit: pd.DataFrame) -> None:
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.5))

    ax = axes[0]
    de = C.read_provided_de("CA", "GCvsFLT-SAL")
    ax.scatter(de["pvalue"], de["padj_deposited"], s=3, c="#7F8C8D", lw=0, alpha=0.4,
               rasterized=True)
    ax.plot([0, 1], [0, 1], color="black", lw=0.9, ls="--")
    ax.set_xlabel("Deposited raw $P$")
    ax.set_ylabel("Deposited 'adjusted $P$'")
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.02, 1)
    frac = float(audit.loc[audit["region_short"] == "CA1", "frac_adj_below_raw_p"].iloc[0])
    ax.text(0.04, 0.96,
            f"{frac*100:.0f}% of targets sit below the\nidentity line: no multiple-testing\n"
            "correction can lower a $P$-value",
            transform=ax.transAxes, fontsize=6.3, va="top")
    P.panel_label(ax, "a")
    ax.set_title("CA1, Flight vs Ground (saline)", pad=3)

    ax = axes[1]
    vals = de["padj_deposited"].dropna()
    ax.hist(vals, bins=np.arange(0, 1.02, 0.01), color="#4C72B0", lw=0)
    ax.set_xlabel("Deposited 'adjusted $P$'")
    ax.set_ylabel("Targets")
    n_uniq = int(audit.loc[audit["region_short"] == "CA1", "n_unique_adj_values"].iloc[0])
    n_zero = int(audit.loc[audit["region_short"] == "CA1", "n_adj_exactly_zero"].iloc[0])
    ax.text(0.96, 0.96,
            f"only {n_uniq} distinct values\n(rounded to 2 dp);\n{n_zero:,} are exactly 0",
            transform=ax.transAxes, fontsize=6.3, va="top", ha="right")
    P.panel_label(ax, "b")
    ax.set_title("Resolution of the deposited FDR", pad=3)

    ax = axes[2]
    idx = np.arange(len(audit))
    ax.bar(idx - 0.2, audit["n_sig_deposited_adj_lt_0.05"], width=0.4, color="#C44E52",
           label="deposited 'adjusted $P$' < 0.05", lw=0)
    ax.bar(idx + 0.2, audit["n_sig_recomputed_BH_lt_0.05"], width=0.4, color="#4C72B0",
           label="Benjamini-Hochberg < 0.05", lw=0)
    ax.set_xticks(idx)
    ax.set_xticklabels(
        [f"{r}\n{c.replace('SALvsBuOE-', 'BuOE ').replace('GCvsFLT-', 'Flight ')}"
         for r, c in zip(audit["region_short"], audit["contrast"])],
        fontsize=5.0,
    )
    ax.set_ylabel("Significant targets")
    ax.legend(fontsize=6.2, loc="upper right")
    ax.grid(axis="y", ls=":", lw=0.5, alpha=0.4)
    P.panel_label(ax, "c")
    ax.set_title("Every blue bar is zero", pad=3)
    fig.tight_layout()
    P.save(fig, "Figure3_deposited_statistics_audit", C.PANEL_DIR)


def figure4(summary: pd.DataFrame, sens: pd.DataFrame) -> None:
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(11.4, 3.6),
                             gridspec_kw={"width_ratios": [1.4, 1, 1]})
    ax = axes[0]
    x = np.arange(len(LANDSCAPE))
    for ri, region in enumerate(C.REGION_ORDER):
        vals = [
            int(summary[(summary["region_short"] == region) &
                        (summary["contrast"] == n)]["nominal_p0.01_n"].iloc[0])
            for n, _ in LANDSCAPE
        ]
        ax.bar(x + (ri - 1.5) * 0.2, vals, width=0.2, color=P.REGION_COLORS[region],
               label=region, lw=0)
    exp_null = float(summary["expected_null_p0.01"].iloc[0])
    ax.axhline(exp_null, color="black", ls="--", lw=0.9)
    ax.text(len(LANDSCAPE) - 0.45, exp_null * 1.03,
            f"null expectation ({exp_null:.0f})", fontsize=6.2, ha="right", va="bottom")
    ax.set_xticks(x)
    ax.set_xticklabels([l for _, l in LANDSCAPE], fontsize=6.0, rotation=30, ha="right")
    ax.set_ylabel("Targets at nominal $P$ < 0.01")
    ax.legend(title="Region", fontsize=6.2, title_fontsize=6.2, ncols=2)
    ax.grid(axis="y", ls=":", lw=0.5, alpha=0.4)
    P.panel_label(ax, "a")
    ax.set_title("Nominal signal vs the global null", pad=3)

    ax = axes[1]
    piv = summary.pivot_table(index="contrast", columns="region_short",
                              values="fdr05_lfc1_n_deg", aggfunc="sum")
    piv = piv.reindex(index=[n for n, _ in LANDSCAPE], columns=C.REGION_ORDER)
    ax.imshow(piv.to_numpy(), cmap="Reds", vmin=0, vmax=1)
    ax.set_xticks(range(len(C.REGION_ORDER)))
    ax.set_xticklabels(C.REGION_ORDER)
    ax.set_yticks(range(len(piv)))
    ax.set_yticklabels([l for _, l in LANDSCAPE], fontsize=6.0)
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            ax.text(j, i, f"{int(piv.to_numpy()[i, j])}", ha="center", va="center",
                    fontsize=7.5)
    P.panel_label(ax, "b", dx=-0.85)
    ax.set_title("FDR-significant targets\n(BH<0.05, |log$_2$FC|$\\geq$1)", pad=3)

    ax = axes[2]
    s = sens[sens["contrast"] == "FLT_vs_GC_in_SAL"].set_index("region_short")
    s = s.reindex(C.REGION_ORDER)
    x = np.arange(len(C.REGION_ORDER))
    ax.bar(x - 0.2, s["mde_log2FC_nominal_a0.05"], width=0.4, color="#95A5A6",
           label="nominal $\\alpha$=0.05", lw=0)
    ax.bar(x + 0.2, s["mde_log2FC_genomewide_a0.05_over_n"], width=0.4, color="#2C3E50",
           label="genome-wide $\\alpha$/n", lw=0)
    ax.axhline(C.LFC_CUTOFF, color="#C44E52", ls="--", lw=0.9)
    ax.text(len(C.REGION_ORDER) - 0.5, C.LFC_CUTOFF * 1.04, "|log$_2$FC| = 1",
            fontsize=6.2, ha="right", va="bottom", color="#C44E52")
    ax.set_xticks(x)
    ax.set_xticklabels(C.REGION_ORDER)
    ax.set_ylabel("Smallest detectable |log$_2$FC|\nat 80% power")
    ax.legend(fontsize=6.2)
    ax.grid(axis="y", ls=":", lw=0.5, alpha=0.4)
    P.panel_label(ax, "c")
    ax.set_title("What this design could detect", pad=3)
    fig.tight_layout()
    P.save(fig, "Figure4_DE_landscape_and_sensitivity", C.PANEL_DIR)


def figure5(long: pd.DataFrame) -> None:
    import matplotlib.pyplot as plt

    rows = [
        ("FLT_vs_GC_in_SAL", "Flight vs Ground (saline)"),
        ("FLT_vs_GC_in_BuOE", "Flight vs Ground (BuOE)"),
        ("interaction", "Flight $\\times$ BuOE interaction"),
    ]
    fig, axes = plt.subplots(len(rows), 4, figsize=(12.2, 8.4))
    for i, (contrast, nice) in enumerate(rows):
        for j, region in enumerate(C.REGION_ORDER):
            ax = axes[i][j]
            g = long[(long["region_short"] == region) & (long["contrast"] == contrast)]
            P.volcano(
                ax,
                g["log2FC"].to_numpy(),
                g["pvalue"].to_numpy(),
                padj=g["padj_BH"].to_numpy(),
                genes=g["gene"].to_numpy(),
                lfc_cut=C.LFC_CUTOFF_RELAXED,
                title=region if i == 0 else "",
                label_n=5,
            )
            if j == 0:
                ax.text(-0.34, 0.5, nice, transform=ax.transAxes, rotation=90,
                        va="center", ha="center", fontsize=8, fontweight="bold")
        P.panel_label(axes[i][0], "abc"[i], dx=-0.24)
    fig.tight_layout()
    P.save(fig, "Figure5_volcanoes", C.PANEL_DIR)


def figure6(att: pd.DataFrame, long: pd.DataFrame) -> None:
    import matplotlib.pyplot as plt

    r = att.set_index("region_short").reindex(C.REGION_ORDER)
    x = np.arange(len(C.REGION_ORDER))
    fig = plt.figure(figsize=(11.0, 6.4))
    gs = fig.add_gridspec(2, 3, hspace=0.42, wspace=0.32)

    # (a) saline vs BuOE fold-change scatter for the strongest region
    ax = fig.add_subplot(gs[0, 0])
    region = "DG"
    a = long[(long["region_short"] == region) & (long["contrast"] == "FLT_vs_GC_in_SAL")]
    b = long[(long["region_short"] == region) & (long["contrast"] == "FLT_vs_GC_in_BuOE")]
    m = a[["gene", "log2FC"]].merge(b[["gene", "log2FC"]], on="gene",
                                    suffixes=("_sal", "_buoe"))
    ax.scatter(m["log2FC_sal"], m["log2FC_buoe"], s=3, c="#7F8C8D", lw=0, alpha=0.35,
               rasterized=True)
    lim = float(np.nanpercentile(np.abs(np.r_[m["log2FC_sal"], m["log2FC_buoe"]]), 99.5))
    ax.plot([-lim, lim], [-lim, lim], color="black", lw=0.9, ls="--",
            label="equal response")
    # The noise-corrected RMS ratio is used rather than a fitted slope: with
    # per-target noise this large, regression slopes (OLS or Deming) are
    # dominated by the assumed error-variance ratio and are not stable.
    slope = float(r.loc[region, "rms_ratio_BuOE_over_saline"])
    ax.plot([-lim, lim], [-lim * slope, lim * slope], color="#C44E52", lw=1.1,
            label=f"noise-corrected ratio {slope:.2f}")
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("log$_2$FC Flight vs Ground (saline)")
    ax.set_ylabel("log$_2$FC Flight vs Ground (BuOE)")
    ax.legend(fontsize=6.0, loc="upper left")
    P.panel_label(ax, "a")
    ax.set_title(f"{region}: per-target comparison", pad=3)

    # (b) noise-corrected magnitudes
    ax = fig.add_subplot(gs[0, 1])
    for off, key, color, lab in [(-0.16, "saline", "#C44E52", "saline"),
                                 (0.16, "BuOE", "#4C72B0", "BuOE")]:
        est = r[f"msq_true_{key}"].to_numpy()
        lo = est - r[f"msq_true_{key}_ci_lo"].to_numpy()
        hi = r[f"msq_true_{key}_ci_hi"].to_numpy() - est
        ax.errorbar(x + off, est, yerr=[lo, hi], fmt="o", color=color, ms=4.5,
                    capsize=2.5, lw=1.0, label=lab)
    ax.axhline(0, color="black", ls="--", lw=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels(C.REGION_ORDER)
    ax.set_ylabel("Noise-corrected mean squared\nspaceflight effect")
    ax.legend(fontsize=6.3)
    ax.grid(axis="y", ls=":", lw=0.5, alpha=0.4)
    P.panel_label(ax, "b")
    ax.set_title("Effect magnitude above noise", pad=3)

    # (c) attenuation test
    ax = fig.add_subplot(gs[0, 2])
    est = r["msq_diff_saline_minus_BuOE"].to_numpy()
    lo = est - r["msq_diff_ci_lo"].to_numpy()
    hi = r["msq_diff_ci_hi"].to_numpy() - est
    ax.errorbar(est, x, xerr=[lo, hi], fmt="o", color="#2C3E50", ms=5, capsize=2.5, lw=1.0)
    ax.axvline(0, color="black", ls="--", lw=0.9)
    ax.set_yticks(x)
    ax.set_yticklabels(C.REGION_ORDER)
    ax.set_xlabel("saline $-$ BuOE mean squared effect\n(>0 = attenuation)")
    ax.grid(axis="x", ls=":", lw=0.5, alpha=0.4)
    span = float(np.nanmax(r["msq_diff_ci_hi"]) - np.nanmin(r["msq_diff_ci_lo"]))
    ax.set_xlim(float(np.nanmin(r["msq_diff_ci_lo"])) - 0.08 * span,
                float(np.nanmax(r["msq_diff_ci_hi"])) + 0.34 * span)
    ax.set_ylim(-0.6, len(C.REGION_ORDER) - 0.4)
    for i, p in enumerate(r["msq_diff_bootstrap_p"]):
        txt = "$p$<0.001" if p < 0.001 else f"$p$={p:.3f}"
        ax.text(0.995, i + 0.30, txt, transform=ax.get_yaxis_transform(),
                ha="right", va="center", fontsize=6.0)
    P.panel_label(ax, "c")
    ax.set_title("Attenuation test (bootstrap 95% CI)", pad=3)

    # (d) variance ratio
    ax = fig.add_subplot(gs[1, 0])
    ax.bar(x, r["median_var_ratio_BuOE_over_saline"], width=0.55, color="#8172B3", lw=0)
    ax.axhline(1.0, color="black", ls="--", lw=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels(C.REGION_ORDER)
    ax.set_ylabel("Within-arm variance ratio\nBuOE / saline")
    ax.set_ylim(0, max(1.6, float(r["median_var_ratio_BuOE_over_saline"].max()) * 1.15))
    ax.grid(axis="y", ls=":", lw=0.5, alpha=0.4)
    P.panel_label(ax, "d")
    ax.set_title("The BuOE arm is noisier", pad=3)

    # (e) verdict depends on the variance assumption
    ax = fig.add_subplot(gs[1, 1])
    ax.bar(x - 0.19, r["msq_diff_pooled_variance"], width=0.36, color="#BDC3C7",
           label="single pooled variance", lw=0)
    ax.bar(x + 0.19, r["msq_diff_saline_minus_BuOE"], width=0.36, color="#2C3E50",
           label="each arm's own variance", lw=0)
    ax.axhline(0, color="black", ls="--", lw=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels(C.REGION_ORDER)
    ax.set_ylabel("saline $-$ BuOE\nmean squared effect")
    ax.legend(fontsize=6.0)
    ax.grid(axis="y", ls=":", lw=0.5, alpha=0.4)
    P.panel_label(ax, "e")
    ax.set_title("The sign depends on the\nvariance assumption", pad=3)

    # (f) selection bias
    ax = fig.add_subplot(gs[1, 2])
    ax.bar(x - 0.31, r["median_abs_lfc_top500_saline_SELECTION_BIASED"], width=0.19,
           color="#C44E52", label="saline, top 500 by saline $P$", lw=0)
    ax.bar(x - 0.11, r["median_abs_lfc_top500_BuOE_SELECTION_BIASED"], width=0.19,
           color="#E8A0A0", label="BuOE, same 500", lw=0)
    ax.bar(x + 0.11, r["median_abs_lfc_random500_saline"], width=0.19, color="#4C72B0",
           label="saline, 500 random", lw=0)
    ax.bar(x + 0.31, r["median_abs_lfc_random500_BuOE"], width=0.19, color="#9DB8D8",
           label="BuOE, same random", lw=0)
    ax.set_xticks(x)
    ax.set_xticklabels(C.REGION_ORDER)
    ax.set_ylabel("Median |log$_2$FC|")
    ax.legend(fontsize=5.3, loc="upper right")
    ax.grid(axis="y", ls=":", lw=0.5, alpha=0.4)
    P.panel_label(ax, "f")
    ax.set_title("Selecting on the saline response\nmanufactures attenuation", pad=3)

    P.save(fig, "Figure6_BuOE_attenuation", C.PANEL_DIR)


def figure7(tstats: pd.DataFrame, cats: pd.DataFrame, enr: pd.DataFrame,
            mapping: pd.DataFrame) -> None:
    """CNS target panel: coverage, panel-restricted hits, category pattern."""
    import matplotlib.pyplot as plt

    CATEGORY_ORDER = [
        "Neuroinflammation / cytokine", "Inflammasome / pyroptosis",
        "Blood-brain barrier / endothelial", "Coagulation / vascular",
        "Neuronal / synaptic", "Glial / myelin", "Injury biomarker",
        "Oxidative stress / metabolism", "Transcription / chromatin / other",
    ]
    CONTRASTS = {
        "FLT_vs_GC_in_SAL": "Flight vs Ground (saline)",
        "FLT_vs_GC_in_BuOE": "Flight vs Ground (BuOE)",
        "BuOE_vs_SAL_in_GC": "BuOE vs saline (ground)",
        "BuOE_vs_SAL_in_FLT": "BuOE vs saline (flight)",
        "spaceflight": "Spaceflight main effect",
        "treatment": "Treatment main effect",
        "interaction": "Flight $\\times$ BuOE interaction",
    }

    fig = plt.figure(figsize=(11.2, 7.4))
    gs = fig.add_gridspec(2, 3, hspace=0.55, wspace=0.42,
                          height_ratios=[1, 1.15])

    # (a) coverage
    ax = fig.add_subplot(gs[0, 0])
    cov = (mapping.groupby("category")
           .agg(n=("analyte", "size"), on_panel=("on_geomx_panel", "sum"))
           .reindex(CATEGORY_ORDER).fillna(0))
    y = np.arange(len(cov))
    ax.barh(y, cov["n"], color="#DDDDDD", height=0.68, label="in target list")
    ax.barh(y, cov["on_panel"], color="#4C72B0", height=0.68, label="measured")
    ax.set_yticks(y)
    ax.set_yticklabels([c.replace(" / ", "/\n") for c in cov.index], fontsize=5.2)
    ax.invert_yaxis()
    ax.set_xlabel("Analytes")
    ax.legend(fontsize=5.8, loc="lower right")
    ax.grid(axis="x", ls=":", lw=0.5, alpha=0.4)
    P.panel_label(ax, "a", dx=-0.52)
    ax.set_title("Panel coverage", pad=3)

    # (b) panel-restricted significant counts
    ax = fig.add_subplot(gs[0, 1])
    piv = tstats.pivot_table(index="contrast", columns="region_short",
                             values="padj_panel",
                             aggfunc=lambda s: int((s < C.FDR_CUTOFF).sum()))
    piv = piv.reindex(index=list(CONTRASTS), columns=C.REGION_ORDER).fillna(0)
    im = ax.imshow(piv.to_numpy(), cmap="YlOrRd", vmin=0,
                   vmax=max(1, np.nanmax(piv.to_numpy())))
    ax.set_xticks(range(len(C.REGION_ORDER)))
    ax.set_xticklabels(C.REGION_ORDER)
    ax.set_yticks(range(len(piv)))
    ax.set_yticklabels([CONTRASTS[c] for c in piv.index], fontsize=5.6)
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            ax.text(j, i, f"{int(piv.to_numpy()[i, j])}", ha="center", va="center",
                    fontsize=7)
    P.panel_label(ax, "b", dx=-0.95)
    ax.set_title("Targets at BH<0.05\nwithin the panel", pad=3)

    # (c) panel volcano for the flight contrast
    ax = fig.add_subplot(gs[0, 2])
    g = tstats[tstats["contrast"] == "FLT_vs_GC_in_SAL"]
    for region in C.REGION_ORDER:
        s = g[g["region_short"] == region]
        ax.scatter(s["log2FC"], -np.log10(s["pvalue"]), s=13,
                   color=P.REGION_COLORS[region], lw=0.3, edgecolor="white",
                   label=region)
    sig = g[g["padj_panel"] < C.FDR_CUTOFF]
    if len(sig):
        ax.scatter(sig["log2FC"], -np.log10(sig["pvalue"]), s=46, facecolors="none",
                   edgecolors="black", linewidths=0.8, zorder=5)
        for r in sig.itertuples():
            ax.annotate(r.analyte, (r.log2FC, -np.log10(r.pvalue)),
                        textcoords="offset points", xytext=(-4, 5), fontsize=5.6,
                        style="italic", ha="right")
    ax.axvline(0, color="#999", lw=0.5, ls=":")
    ax.set_xlabel("log$_2$ fold change")
    ax.set_ylabel("$-$log$_{10}$ nominal $P$")
    ax.legend(fontsize=5.6, title="Region", title_fontsize=5.6)
    P.panel_label(ax, "c")
    ax.set_title("Flight vs Ground (saline),\npanel targets only", pad=3)

    # (d, e) category mean-t heat maps for two contrasts
    vmax = float(np.nanpercentile(np.abs(cats["mean_t"]), 98)) or 1.0
    for k, contrast in enumerate(["FLT_vs_GC_in_SAL", "BuOE_vs_SAL_in_FLT"]):
        ax = fig.add_subplot(gs[1, k])
        sub = cats[cats["contrast"] == contrast]
        piv = sub.pivot_table(index="category", columns="region_short", values="mean_t")
        piv = piv.reindex(index=CATEGORY_ORDER, columns=C.REGION_ORDER)
        pq = sub.pivot_table(index="category", columns="region_short",
                             values="p_mean_t_vs_0").reindex(
            index=CATEGORY_ORDER, columns=C.REGION_ORDER)
        im2 = ax.imshow(piv.to_numpy(), cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                        aspect="auto", interpolation="nearest")
        for i in range(piv.shape[0]):
            for j in range(piv.shape[1]):
                p = pq.to_numpy()[i, j]
                if np.isfinite(p) and p < 0.05:
                    ax.text(j, i, "*", ha="center", va="center", fontsize=8)
        ax.set_xticks(range(len(C.REGION_ORDER)))
        ax.set_xticklabels(C.REGION_ORDER)
        if k == 0:
            ax.set_yticks(range(len(CATEGORY_ORDER)))
            ax.set_yticklabels([c.replace(" / ", "/\n") for c in CATEGORY_ORDER],
                               fontsize=5.2)
            P.panel_label(ax, "d", dx=-0.52)
        else:
            ax.set_yticks([])
            P.panel_label(ax, "e", dx=-0.08)
        ax.set_title(CONTRASTS[contrast] + "\n(* mean $t\\neq$0 at $P$<0.05)",
                     fontsize=7.2, pad=3)
        if k == 1:
            cb = fig.colorbar(im2, ax=ax, fraction=0.05, pad=0.03)
            cb.set_label("Mean $t$", fontsize=6.5)
            cb.ax.tick_params(labelsize=6)

    # (f) set enrichment permutation p
    ax = fig.add_subplot(gs[1, 2])
    piv = enr.pivot_table(index="contrast", columns="region_short",
                          values="permutation_p_delta_abs_t")
    piv = piv.reindex(index=list(CONTRASTS), columns=C.REGION_ORDER)
    im3 = ax.imshow(piv.to_numpy(), cmap="viridis_r", vmin=0, vmax=1)
    ax.set_xticks(range(len(C.REGION_ORDER)))
    ax.set_xticklabels(C.REGION_ORDER)
    ax.set_yticks(range(len(piv)))
    ax.set_yticklabels([CONTRASTS[c] for c in piv.index], fontsize=5.6)
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            v = piv.to_numpy()[i, j]
            if np.isfinite(v):
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=6,
                        color="white" if v < 0.5 else "black")
    P.panel_label(ax, "f", dx=-0.95)
    ax.set_title("Panel vs transcriptome,\npermutation $P$", pad=3)

    P.save(fig, "Figure7_CNS_target_panel", C.PANEL_DIR)


def main() -> int:
    C.ensure_dirs()
    C.PANEL_DIR.mkdir(parents=True, exist_ok=True)
    P.use_style()

    qc = need("qc_metrics.csv")
    cohort = need("cohort_summary.csv")
    coords = need("PCA_coordinates.csv")
    var = need("PCA_variance.csv")
    vp = need("variance_partition.csv")
    audit = need("deposited_stats_audit.csv")
    summary = need("DEG_summary_all.csv")
    sens = need("sensitivity_minimum_detectable_effect.csv")
    long = need("DE_all_contrasts_long.csv")
    att = need("attenuation_summary.csv")

    figure1(qc, cohort)
    figure2(coords, var, vp)
    figure3(audit)
    figure4(summary, sens)
    figure5(long)
    figure6(att, long)
    figure7(
        need("target_stats_all.csv"),
        need("target_category_summary.csv"),
        need("target_setenrichment.csv"),
        need("target_panel_mapping.csv"),
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
