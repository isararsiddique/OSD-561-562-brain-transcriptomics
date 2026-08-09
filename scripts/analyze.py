#!/usr/bin/env python3
"""CNS Injury EV manuscript -- bulk RNA-seq analysis for OSD-561 & OSD-562.

Uses the GeneLab-processed DESeq2 tables (no raw-read processing needed) to
produce, for each study and for the combined (Cb + HPC) dataset:

  1. PCA of global transcriptomic profiles (VST counts)
  2. Volcano plots for the spaceflight effect (Space Flight vs Ground Control,
     matched on age and environment)
  3. Heat maps of a user-supplied focus gene list (VST z-scores across samples)
  4. A focus-gene log2FC matrix across the spaceflight contrasts (+ heat map)
  5. A tidy table of focus-gene statistics for every spaceflight contrast

Datasets
  OSD-561 / GLDS-556 : Cerebellum (Cb)
  OSD-562 / GLDS-557 : Hippocampus (HPC)
Both are RRRM2 (Rodent Research Reference Mission 2) spaceflight studies.

Run:
  ./.venv/bin/python scripts/analyze.py
"""
from __future__ import annotations

import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #
ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
FIG = ROOT / "results" / "figures"
TAB = ROOT / "results" / "tables"
GENE_FILE = ROOT / "config" / "focus_genes.txt"

STUDIES = {
    "OSD-561": {"glds": "GLDS-556", "region": "Cb", "label": "Cerebellum"},
    "OSD-562": {"glds": "GLDS-557", "region": "HPC", "label": "Hippocampus"},
}

N_TOP_VARIABLE = 2000     # genes used for PCA
LFC_THRESH = 1.0          # |log2FC| line on volcano
ADJP_THRESH = 0.05        # adj p-value significance line/threshold

# Shared, colour-blind-friendly categorical palettes (consistent across figs)
PAL_TRT = {"Space Flight": "#c0392b", "Ground Control": "#2c6fb0"}
PAL_AGE = {"OLD": "#6a51a3", "YNG": "#2e8b7f"}
PAL_REG = {"Cb": "#e08214", "HPC": "#1f6fb2"}
NICE = {"treatment_full": "Treatment", "region": "Region", "age": "Age"}

# Contrast-category palette (used by the all-contrast / DEG-landscape figures)
PAL_CAT = {
    "Spaceflight": "#c0392b",   # Space Flight vs Ground Control (matched)
    "Age": "#6a51a3",           # 29 week vs 12 week (matched)
    "Environment": "#e08214",   # On ISS vs On Earth (matched)
    "Confounded": "#7f8c8d",    # >1 factor differs (not cleanly interpretable)
}
# Compact labels for group descriptors
ABBR = {
    "12 week": "12wk", "29 week": "29wk",
    "Space Flight": "SF", "Ground Control": "GC",
    "On Earth": "Earth", "On ISS": "ISS",
}


def _ab(text: str) -> str:
    for k, v in ABBR.items():
        text = text.replace(k, v)
    return text


def set_journal_style() -> None:
    """Consistent publication styling for every figure."""
    sns.set_theme(style="ticks", context="paper")
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8,
        "axes.titlesize": 10,
        "axes.titleweight": "bold",
        "axes.labelsize": 9,
        "axes.linewidth": 0.8,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "legend.fontsize": 7.5,
        "legend.title_fontsize": 8,
        "figure.dpi": 150,
        "savefig.dpi": 400,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.06,
        "pdf.fonttype": 42,   # keep text editable in Illustrator / Inkscape
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })


# --------------------------------------------------------------------------- #
# Loading helpers
# --------------------------------------------------------------------------- #
def load_counts(study: str, kind: str) -> pd.DataFrame:
    """kind in {'VST_Counts', 'Normalized_Counts'}; returns genes x samples."""
    glds = STUDIES[study]["glds"]
    path = DATA / study / f"{glds}_rna_seq_{kind}_GLbulkRNAseq.csv"
    df = pd.read_csv(path)
    df = df.rename(columns={df.columns[0]: "ENSEMBL"}).set_index("ENSEMBL")
    return df


def load_dge(study: str) -> pd.DataFrame:
    glds = STUDIES[study]["glds"]
    path = DATA / study / f"{glds}_rna_seq_differential_expression_GLbulkRNAseq.csv"
    return pd.read_csv(path, low_memory=False)


def sample_metadata(samples) -> pd.DataFrame:
    """Parse RRRM2_<region>_<FLT|GC>_<ISS-T|LAR>_<OLD|YNG>_<animal>."""
    rows = []
    for s in samples:
        t = s.split("_")
        rows.append(
            {
                "sample": s,
                "region": t[1] if len(t) > 1 else "NA",
                "treatment": t[2] if len(t) > 2 else "NA",   # FLT / GC
                "environment": t[3] if len(t) > 3 else "NA",  # ISS-T / LAR
                "age": t[4] if len(t) > 4 else "NA",          # OLD / YNG
                "animal": t[5] if len(t) > 5 else "NA",
            }
        )
    meta = pd.DataFrame(rows).set_index("sample")
    meta["treatment_full"] = meta["treatment"].map(
        {"FLT": "Space Flight", "GC": "Ground Control"}
    ).fillna(meta["treatment"])
    return meta


def read_focus_genes() -> list[str]:
    genes = []
    for line in GENE_FILE.read_text().splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            genes.append(line)
    return genes


# --------------------------------------------------------------------------- #
# Contrast parsing (spaceflight effect = Space Flight vs Ground Control,
# matched on age and environment)
# --------------------------------------------------------------------------- #
def dge_groups(dge: pd.DataFrame) -> set[str]:
    groups = set()
    for col in dge.columns:
        if col.startswith("Log2fc_"):
            a, b = split_contrast(col[len("Log2fc_"):])
            groups.add(a)
            groups.add(b)
    return groups


def split_contrast(name: str) -> tuple[str, str]:
    """'(A)v(B)' -> ('A', 'B')."""
    a, b = name.split(")v(")
    return a.lstrip("("), b.rstrip(")")


def parse_group(g: str) -> dict:
    parts = [p.strip() for p in g.split("&")]
    # e.g. ['12 week', 'Ground Control', 'On Earth']
    return {"age": parts[0], "treatment": parts[1], "environment": parts[2]}


def spaceflight_contrasts(dge: pd.DataFrame) -> list[dict]:
    """Return matched Flight-vs-Ground contrasts oriented Flight - Ground."""
    groups = dge_groups(dge)
    parsed = {g: parse_group(g) for g in groups}
    out = []
    seen = set()
    for g, p in parsed.items():
        if p["treatment"] != "Space Flight":
            continue
        # find the Ground Control group with same age & environment
        for g2, p2 in parsed.items():
            if (
                p2["treatment"] == "Ground Control"
                and p2["age"] == p["age"]
                and p2["environment"] == p["environment"]
            ):
                key = (p["age"], p["environment"])
                if key in seen:
                    continue
                seen.add(key)
                out.append(
                    {
                        "flight": g,
                        "ground": g2,
                        "age": p["age"],
                        "environment": p["environment"],
                        "short": f"{p['age']} | {p['environment']}".replace(" week", "wk"),
                    }
                )
    return sorted(out, key=lambda d: (d["age"], d["environment"]))


def contrast_stats_oriented(dge: pd.DataFrame, g_num: str, g_den: str) -> pd.DataFrame:
    """Return DE stats oriented as (g_num - g_den): positive log2FC = up in g_num.

    Resolves whichever direction the GeneLab table stored the contrast in and
    flips the sign if necessary, so callers get a consistent orientation.
    """
    lfc_fwd = f"Log2fc_({g_num})v({g_den})"
    lfc_rev = f"Log2fc_({g_den})v({g_num})"
    if lfc_fwd in dge.columns:
        base, sign = f"({g_num})v({g_den})", 1.0
    elif lfc_rev in dge.columns:
        base, sign = f"({g_den})v({g_num})", -1.0
    else:
        raise KeyError(f"No contrast column for {g_num} vs {g_den}")
    return pd.DataFrame(
        {
            "ENSEMBL": dge["ENSEMBL"],
            "SYMBOL": dge["SYMBOL"],
            "log2FC": sign * dge[f"Log2fc_{base}"],
            "pvalue": dge[f"P.value_{base}"],
            "padj": dge[f"Adj.p.value_{base}"],
        }
    )


def contrast_stats(dge: pd.DataFrame, flight: str, ground: str) -> pd.DataFrame:
    """log2FC oriented Flight - Ground (kept for the spaceflight-only pipeline)."""
    return contrast_stats_oriented(dge, flight, ground)


def all_unique_contrasts(dge: pd.DataFrame) -> list[dict]:
    """Every unique group comparison in the DE table (reverse-direction
    duplicates collapsed), each categorised and consistently oriented.

    Categories (by how many of age / treatment / environment differ):
      Spaceflight  - only treatment differs  -> oriented Space Flight - Ground
      Age          - only age differs         -> oriented 29 week - 12 week
      Environment  - only environment differs -> oriented On ISS - On Earth
      Confounded   - >1 factor differs         -> native orientation
    """
    bases = [c[len("Log2fc_"):] for c in dge.columns if c.startswith("Log2fc_")]
    seen: set[frozenset] = set()
    out: list[dict] = []
    for base in bases:
        a, b = split_contrast(base)
        key = frozenset((a, b))
        if key in seen:
            continue                      # reverse-direction duplicate
        seen.add(key)
        pa, pb = parse_group(a), parse_group(b)
        diffs = [f for f in ("age", "treatment", "environment") if pa[f] != pb[f]]

        if len(diffs) == 1:
            f = diffs[0]
            if f == "treatment":
                num, den = (a, b) if pa["treatment"] == "Space Flight" else (b, a)
                cat, unit = "Spaceflight", "SF / GC"
            elif f == "age":
                num, den = (a, b) if pa["age"] == "29 week" else (b, a)
                cat, unit = "Age", "Old / Young"
            else:
                num, den = (a, b) if pa["environment"] == "On ISS" else (b, a)
                cat, unit = "Environment", "ISS / Earth"
        else:
            cat, unit, num, den = "Confounded", "num / den", a, b

        pn, pd_ = parse_group(num), parse_group(den)
        if cat == "Spaceflight":
            label = f"SF vs GC | {_ab(pn['age'])} · {_ab(pn['environment'])}"
        elif cat == "Age":
            label = f"Old vs Young | {_ab(pn['treatment'])} · {_ab(pn['environment'])}"
        elif cat == "Environment":
            label = f"ISS vs Earth | {_ab(pn['age'])} · {_ab(pn['treatment'])}"
        else:
            label = f"{_ab(num).replace(' & ', '/')}  vs  {_ab(den).replace(' & ', '/')}"

        out.append({
            "category": cat, "num": num, "den": den, "unit": unit,
            "label": label, "n_factors": len(diffs),
            "factors_differ": ",".join(diffs),
            "tag": f"{cat}_{_ab(num).replace(' & ', '-').replace(' ', '')}"
                   f"__{_ab(den).replace(' & ', '-').replace(' ', '')}",
        })

    cat_order = {"Spaceflight": 0, "Age": 1, "Environment": 2, "Confounded": 3}
    return sorted(out, key=lambda d: (cat_order[d["category"]], d["label"]))


def deg_counts(stats: pd.DataFrame) -> dict:
    """Up / down / total DEGs at FDR < ADJP_THRESH and |log2FC| >= LFC_THRESH."""
    v = stats.dropna(subset=["log2FC", "padj"])
    sig = (v["padj"] < ADJP_THRESH) & (v["log2FC"].abs() >= LFC_THRESH)
    return {
        "genes_tested": int(v.shape[0]),
        "up": int((sig & (v["log2FC"] > 0)).sum()),
        "down": int((sig & (v["log2FC"] < 0)).sum()),
        "total": int(sig.sum()),
    }


# --------------------------------------------------------------------------- #
# Plots
# --------------------------------------------------------------------------- #
def _pca_frame(vst: pd.DataFrame, n_top: int) -> tuple[pd.DataFrame, np.ndarray]:
    variances = vst.var(axis=1).sort_values(ascending=False)
    top = variances.head(n_top).index
    X = vst.loc[top].T.values
    X = StandardScaler().fit_transform(X)
    pca = PCA(n_components=2, random_state=0)
    coords = pca.fit_transform(X)
    frame = pd.DataFrame(coords, columns=["PC1", "PC2"], index=vst.columns)
    return frame, pca.explained_variance_ratio_


def plot_pca(vst, meta, title, out_stem, hue="treatment_full", style="age",
             extra_hue=None):
    frame, evr = _pca_frame(vst, N_TOP_VARIABLE)
    frame = frame.join(meta).rename(columns=NICE)
    hue_col = NICE.get(extra_hue or hue, extra_hue or hue)
    style_col = NICE.get(style, style)
    palette = {"Treatment": PAL_TRT, "Region": PAL_REG, "Age": PAL_AGE}.get(hue_col)

    fig, ax = plt.subplots(figsize=(4.4, 3.5))
    ax.axhline(0, color="0.88", lw=0.7, zorder=0)
    ax.axvline(0, color="0.88", lw=0.7, zorder=0)
    sns.scatterplot(
        data=frame, x="PC1", y="PC2", hue=hue_col, style=style_col,
        palette=palette, s=55, ax=ax, edgecolor="black", linewidth=0.4,
        alpha=0.95, zorder=3,
    )
    ax.set_xlabel(f"PC1 ({evr[0] * 100:.1f}%)")
    ax.set_ylabel(f"PC2 ({evr[1] * 100:.1f}%)")
    ax.set_title(title, pad=8)
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False,
              handletextpad=0.4, labelspacing=0.3, borderaxespad=0.0)
    sns.despine(ax=ax)
    _save(fig, out_stem)


def plot_volcano(stats: pd.DataFrame, focus: set[str], title: str, out_stem: str,
                 up_label: str = "Up in flight", down_label: str = "Down in flight",
                 xlabel: str = "log$_2$ fold change  (Flight / Ground)"):
    d = stats.dropna(subset=["log2FC", "pvalue"]).copy()
    d["neglog10p"] = -np.log10(d["pvalue"].clip(lower=1e-300))
    sig = (d["padj"] < ADJP_THRESH) & (d["log2FC"].abs() >= LFC_THRESH)
    up = sig & (d["log2FC"] > 0)
    down = sig & (d["log2FC"] < 0)

    fig, ax = plt.subplots(figsize=(4.6, 4.2))
    ax.scatter(d.loc[~sig, "log2FC"], d.loc[~sig, "neglog10p"],
               s=6, c="#dcdcdc", alpha=0.7, linewidths=0, rasterized=True)
    ax.scatter(d.loc[down, "log2FC"], d.loc[down, "neglog10p"], s=8,
               c=PAL_TRT["Ground Control"], alpha=0.8, linewidths=0,
               rasterized=True, label=f"{down_label} (n={int(down.sum())})")
    ax.scatter(d.loc[up, "log2FC"], d.loc[up, "neglog10p"], s=8,
               c=PAL_TRT["Space Flight"], alpha=0.8, linewidths=0,
               rasterized=True, label=f"{up_label} (n={int(up.sum())})")

    # focus genes: on top, with repelled italic labels
    fg = d.loc[d["SYMBOL"].str.upper().isin({g.upper() for g in focus})]
    ax.scatter(fg["log2FC"], fg["neglog10p"], s=40, c="#f2c744",
               edgecolor="black", linewidth=0.6, zorder=6,
               label="Focus genes")
    texts = [ax.text(r["log2FC"], r["neglog10p"], r["SYMBOL"], fontsize=7,
                     fontstyle="italic", zorder=7) for _, r in fg.iterrows()]
    if texts:
        try:
            from adjustText import adjust_text
            adjust_text(texts, ax=ax,
                        arrowprops=dict(arrowstyle="-", color="0.45", lw=0.5),
                        expand=(1.3, 1.7), min_arrow_len=3,
                        force_text=(0.4, 0.6))
        except Exception:
            pass

    ax.axhline(-np.log10(ADJP_THRESH), ls=(0, (4, 3)), c="0.6", lw=0.7)
    ax.axvline(LFC_THRESH, ls=(0, (4, 3)), c="0.6", lw=0.7)
    ax.axvline(-LFC_THRESH, ls=(0, (4, 3)), c="0.6", lw=0.7)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("$-$log$_{10}$ $P$")
    ax.set_title(title, pad=8)
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False,
              title="FDR < 0.05, |LFC| ≥ 1", alignment="left",
              handletextpad=0.3, labelspacing=0.4)
    sns.despine(ax=ax)
    _save(fig, out_stem)


def _focus_zscores(vst, dge, focus):
    """Return a genes(symbol) x samples z-score matrix for the focus list."""
    ens2sym = dge.set_index("ENSEMBL")["SYMBOL"]
    sym2ens = {}
    for ens, sym in ens2sym.items():
        if isinstance(sym, str):
            sym2ens.setdefault(sym.upper(), ens)
    rows, labels = [], []
    for gsym in focus:                       # preserve the user's ordering
        ens = sym2ens.get(gsym.upper())
        if ens is not None and ens in vst.index:
            rows.append(ens)
            labels.append(gsym)
    if not rows:
        return None
    sub = vst.loc[rows]
    sub.index = labels
    z = sub.sub(sub.mean(axis=1), axis=0).div(
        sub.std(axis=1).replace(0, np.nan), axis=0)
    return z.dropna(how="all")


def plot_focus_heatmap(vst, meta, dge, focus, title, out_stem):
    from matplotlib.colors import to_rgb
    from matplotlib.patches import Patch
    from scipy.cluster.hierarchy import leaves_list, linkage
    from scipy.spatial.distance import pdist

    z = _focus_zscores(vst, dge, focus)
    if z is None or z.empty:
        print(f"  [heatmap] no focus genes found for {out_stem}; skipping")
        return

    # order samples into biologically meaningful, contiguous groups
    m = meta.loc[z.columns]
    sample_order = m.sort_values(["region", "age", "treatment_full"]).index.tolist()
    z = z[sample_order]
    m = m.loc[sample_order]

    # cluster genes (reorder rows) without drawing a dendrogram
    if z.shape[0] > 2:
        z = z.iloc[leaves_list(linkage(pdist(z.values), method="average"))]

    n_genes, n_samp = z.shape
    tracks = [("Region", "region", PAL_REG),
              ("Age", "age", PAL_AGE),
              ("Treatment", "treatment_full", PAL_TRT)]

    fig_w = 5.6 + 0.055 * n_samp
    fig_h = 1.9 + 0.30 * n_genes
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = fig.add_gridspec(
        nrows=2, ncols=2,
        height_ratios=[0.32 * len(tracks), max(1.2, 0.30 * n_genes)],
        width_ratios=[max(3.0, 0.055 * n_samp), 1.85],
        hspace=0.04, wspace=0.04,
    )
    ax_ann = fig.add_subplot(gs[0, 0])
    ax_hm = fig.add_subplot(gs[1, 0])
    gs_r = gs[:, 1].subgridspec(nrows=4, ncols=1,
                                height_ratios=[1, 1, 1, 1.5], hspace=0.55)
    leg_axes = [fig.add_subplot(gs_r[i]) for i in range(3)]
    ax_cb = fig.add_subplot(gs_r[3])
    for a in (*leg_axes, ax_cb):
        a.axis("off")

    # ---- annotation strips (Region / Age / Treatment) ----
    ann_rgb = np.zeros((len(tracks), n_samp, 3))
    for i, (_, col, pal) in enumerate(tracks):
        for j, s in enumerate(z.columns):
            ann_rgb[i, j] = to_rgb(pal.get(m.loc[s, col], "#dddddd"))
    ax_ann.imshow(ann_rgb, aspect="auto", interpolation="none")
    ax_ann.set_yticks(range(len(tracks)))
    ax_ann.set_yticklabels([t[0] for t in tracks], fontsize=7.5)
    ax_ann.set_xticks([])
    ax_ann.tick_params(length=0)
    for sp in ax_ann.spines.values():
        sp.set_visible(False)

    # ---- heat map ----
    im = ax_hm.imshow(z.values, aspect="auto", cmap="RdBu_r",
                      vmin=-2.5, vmax=2.5, interpolation="none")
    ax_hm.set_yticks(range(n_genes))
    ax_hm.set_yticklabels(z.index, fontstyle="italic", fontsize=8)
    ax_hm.set_xticks([])
    ax_hm.tick_params(length=0)
    for sp in ax_hm.spines.values():
        sp.set_visible(False)

    # white separators between sample groups (both strips and heatmap)
    grp = m[["region", "age", "treatment_full"]].astype(str).agg("|".join, axis=1).values
    for j in range(1, n_samp):
        if grp[j] != grp[j - 1]:
            ax_hm.axvline(j - 0.5, color="white", lw=1.1)
            ax_ann.axvline(j - 0.5, color="white", lw=1.1)

    # ---- category legends (one per annotation track) ----
    for a, (name, col, pal) in zip(leg_axes, tracks):
        present = [v for v in pal if (m[col] == v).any()]
        handles = [Patch(facecolor=pal[v], edgecolor="black", linewidth=0.3,
                         label=str(v)) for v in present]
        a.legend(handles=handles, title=name, loc="center left",
                 bbox_to_anchor=(0.0, 0.5), frameon=False, handlelength=1.1,
                 handleheight=1.1, labelspacing=0.3, borderaxespad=0.0)

    # ---- colorbar ----
    cax = ax_cb.inset_axes([0.0, 0.12, 0.32, 0.72])
    cb = fig.colorbar(im, cax=cax, ticks=[-2, -1, 0, 1, 2])
    cb.set_label("VST z-score", fontsize=8)
    cb.ax.tick_params(labelsize=7, length=2)
    cb.outline.set_linewidth(0.6)

    fig.suptitle(title, x=0.5, y=1.03, fontsize=10, fontweight="bold")
    _save(fig, out_stem)


def plot_lfc_matrix(lfc_df: pd.DataFrame, title: str, out_stem: str):
    """lfc_df: focus genes (rows) x contrasts (cols) of log2FC."""
    lfc_df = lfc_df.dropna(how="all")
    if lfc_df.empty:
        print(f"  [lfc-matrix] empty for {out_stem}; skipping")
        return
    n_g, n_c = lfc_df.shape
    fig, ax = plt.subplots(figsize=(2.8 + 1.05 * n_c, 1.7 + 0.34 * n_g))
    vmax = float(np.nanmax(np.abs(lfc_df.values))) or 1.0
    vmax = max(1.0, np.ceil(vmax * 10) / 10)
    sns.heatmap(lfc_df, cmap="RdBu_r", center=0, vmin=-vmax, vmax=vmax,
                annot=True, fmt=".2f", annot_kws={"size": 7},
                linewidths=0.6, linecolor="white", square=False,
                cbar_kws={"label": "log$_2$FC (Flight / Ground)", "shrink": 0.7,
                          "aspect": 12}, ax=ax)
    ax.set_title(title, pad=8)
    ax.set_xlabel("")
    ax.set_ylabel("")
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right",
             rotation_mode="anchor", fontsize=8)
    plt.setp(ax.get_yticklabels(), rotation=0, fontstyle="italic", fontsize=8)
    ax.tick_params(length=0)
    _save(fig, out_stem)


def _save(fig, stem: str):
    """Save a 400-dpi PNG (review) and a vector PDF (submission)."""
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / f"{stem}.png", dpi=400, bbox_inches="tight")
    fig.savefig(FIG / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  [fig] {stem}.png / .pdf")


# --------------------------------------------------------------------------- #
# Per-study and combined drivers
# --------------------------------------------------------------------------- #
def analyse_study(study: str, focus: list[str]) -> dict:
    info = STUDIES[study]
    print(f"\n=== {study} ({info['glds']}, {info['label']}) ===")
    vst = load_counts(study, "VST_Counts")
    dge = load_dge(study)
    meta = sample_metadata(vst.columns)

    # 1. PCA
    plot_pca(
        vst, meta,
        title=f"{study} {info['label']} — PCA\n(top {N_TOP_VARIABLE} variable genes)",
        out_stem=f"{study}_PCA",
    )

    # 2. spaceflight contrasts -> volcano + collect focus-gene log2FC
    contrasts = spaceflight_contrasts(dge)
    print(f"  spaceflight (Flight vs Ground) contrasts found: {len(contrasts)}")
    focus_up = {g.upper() for g in focus}
    lfc_cols = {}
    stat_rows = []
    for c in contrasts:
        stats = contrast_stats(dge, c["flight"], c["ground"])
        tag = c["short"].replace(" ", "").replace("|", "_")
        plot_volcano(
            stats, set(focus),
            title=f"{study} {info['label']}\nSpace Flight vs Ground Control ({c['short']})",
            out_stem=f"{study}_volcano_{tag}",
        )
        fsub = stats[stats["SYMBOL"].str.upper().isin(focus_up)]
        lfc_cols[c["short"]] = (
            fsub.set_index("SYMBOL")["log2FC"].groupby(level=0).first()
        )
        fs = fsub.copy()
        fs.insert(0, "contrast", c["short"])
        fs.insert(0, "study", study)
        stat_rows.append(fs)

    # 3. focus-gene heatmap (VST z-scores)
    plot_focus_heatmap(
        vst, meta, dge, focus,
        title=f"{study} {info['label']} — focus genes (VST z-score)",
        out_stem=f"{study}_focus_heatmap",
    )

    # 4. focus-gene log2FC matrix across spaceflight contrasts
    if lfc_cols:
        lfc_df = pd.DataFrame(lfc_cols)
        plot_lfc_matrix(
            lfc_df,
            title=f"{study} {info['label']} — focus-gene log$_2$FC",
            out_stem=f"{study}_focus_log2FC_matrix",
        )

    stats_all = pd.concat(stat_rows, ignore_index=True) if stat_rows else pd.DataFrame()
    return {"vst": vst, "meta": meta, "dge": dge, "focus_stats": stats_all}


def analyse_combined(results: dict, focus: list[str]):
    print("\n=== COMBINED (Cb + HPC) ===")
    vsts, metas = [], []
    for study, r in results.items():
        vsts.append(r["vst"])
        metas.append(r["meta"])
    # inner join on shared ENSEMBL ids
    combined = pd.concat(vsts, axis=1, join="inner")
    meta = pd.concat(metas)
    print(f"  combined matrix: {combined.shape[0]} shared genes x {combined.shape[1]} samples")

    plot_pca(
        combined, meta,
        title=f"Combined Cb + HPC — PCA\n(top {N_TOP_VARIABLE} variable genes)",
        out_stem="COMBINED_PCA_by_region", extra_hue="region",
    )
    plot_pca(
        combined, meta,
        title="Combined Cb + HPC — PCA\n(coloured by treatment)",
        out_stem="COMBINED_PCA_by_treatment", hue="treatment_full", style="region",
    )
    # combined focus-gene heatmap (use OSD-561 dge for ENSEMBL->SYMBOL mapping)
    any_dge = next(iter(results.values()))["dge"]
    plot_focus_heatmap(
        combined, meta, any_dge, focus,
        title="Combined Cb + HPC — focus genes (VST z-score)",
        out_stem="COMBINED_focus_heatmap",
    )


def main() -> int:
    set_journal_style()
    for d in (FIG, TAB):
        d.mkdir(parents=True, exist_ok=True)
    focus = read_focus_genes()
    print(f"Focus genes ({len(focus)}): {', '.join(focus)}")

    results = {}
    all_stats = []
    for study in STUDIES:
        r = analyse_study(study, focus)
        results[study] = r
        if not r["focus_stats"].empty:
            all_stats.append(r["focus_stats"])

    analyse_combined(results, focus)

    if all_stats:
        stats = pd.concat(all_stats, ignore_index=True)
        stats = stats[["study", "contrast", "SYMBOL", "ENSEMBL", "log2FC", "pvalue", "padj"]]
        out = TAB / "focus_gene_stats_by_contrast.csv"
        stats.to_csv(out, index=False)
        print(f"\n[table] {out.relative_to(ROOT)}  ({len(stats)} rows)")

    print("\nAll figures written to results/figures/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
