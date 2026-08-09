#!/usr/bin/env python3
"""Global expression structure: normalisation checks, PCA, clustering, variance.

Outputs
  results/tables/csv/normalisation_check.csv    cross-file Q3 comparability
  results/tables/csv/PCA_coordinates.csv        combined + per-region PCA
  results/tables/csv/PCA_variance.csv           variance explained
  results/tables/csv/variance_partition.csv     region vs flight vs drug
  results/tables/csv/sample_correlation.csv     ROI x ROI Pearson r
  results/figures/PCA_*.{pdf,png}
  results/figures/Clustering_*.{pdf,png}

Usage:  python scripts/analyze_expression.py
"""
from __future__ import annotations

import sys

import numpy as np
import pandas as pd
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform
from sklearn.decomposition import PCA

import common as C
import plotting as P

NEG_PROBE = "NegProbe-WTX"
N_TOP = 2000


def load_all(samples: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Concatenate the four per-region matrices into one 48-ROI matrix."""
    mats, metas = [], []
    for study in C.STUDIES:
        log2, sub = C.load_expression(study.code, samples)
        mats.append(log2.drop(index=[NEG_PROBE], errors="ignore"))
        metas.append(sub)
    shared = mats[0].index
    for m in mats[1:]:
        shared = shared.intersection(m.index)
    combined = pd.concat([m.loc[shared] for m in mats], axis=1)
    meta = (
        pd.concat(metas, ignore_index=True)
        .set_index("sample")
        .loc[list(combined.columns)]
        .rename_axis("sample")
        .reset_index()
    )
    return combined, meta


def normalisation_check(samples: pd.DataFrame) -> pd.DataFrame:
    """Q3-normalised files are only comparable across regions if the per-ROI
    upper quartiles agree; verify rather than assume."""
    rows = []
    for study in C.STUDIES:
        mat = C.read_q3(study.code).drop(index=[NEG_PROBE], errors="ignore")
        q3 = mat.quantile(0.75)
        rows.append(
            {
                "region_short": study.short,
                "osd": study.osd,
                "n_targets": int(mat.shape[0]),
                "q3_min": float(q3.min()),
                "q3_median": float(q3.median()),
                "q3_max": float(q3.max()),
                "q3_cv": float(q3.std() / q3.mean()),
                "library_median_min": float(mat.median().min()),
                "library_median_max": float(mat.median().max()),
            }
        )
    return pd.DataFrame(rows)


def top_variable(mat: pd.DataFrame, n: int = N_TOP) -> pd.DataFrame:
    v = mat.var(axis=1).sort_values(ascending=False)
    return mat.loc[v.index[: min(n, len(v))]]


def run_pca(mat: pd.DataFrame, meta: pd.DataFrame, scope: str, n_comp: int = 4):
    sub = top_variable(mat)
    X = sub.T.to_numpy()
    X = X - X.mean(axis=0, keepdims=True)
    n_comp = min(n_comp, X.shape[0] - 1, X.shape[1])
    pca = PCA(n_components=n_comp, random_state=0)
    scores = pca.fit_transform(X)
    coords = pd.DataFrame(
        scores, columns=[f"PC{i+1}" for i in range(n_comp)], index=sub.columns
    ).reset_index(names="sample")
    coords = coords.merge(
        meta[["sample", "region_short", "group", "spaceflight", "treatment", "rep", "slide"]],
        on="sample",
    )
    coords.insert(0, "scope", scope)
    var = pd.DataFrame(
        {
            "scope": scope,
            "PC": [f"PC{i+1}" for i in range(n_comp)],
            "variance_explained": pca.explained_variance_ratio_,
            "n_targets_used": sub.shape[0],
            "n_samples": sub.shape[1],
        }
    )
    return coords, var


def variance_partition(mat: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    """Fraction of each target's variance attributable to region / flight / drug.

    One-way sums of squares for each factor taken separately (they are mutually
    balanced by design, so the decomposition is orthogonal up to the interaction).
    """
    m = meta.set_index("sample").loc[mat.columns]
    X = mat.to_numpy()
    grand = X.mean(axis=1, keepdims=True)
    ss_total = ((X - grand) ** 2).sum(axis=1)

    def ss_factor(labels: pd.Series) -> np.ndarray:
        out = np.zeros(X.shape[0])
        for level in labels.unique():
            idx = np.where((labels == level).to_numpy())[0]
            grp = X[:, idx].mean(axis=1, keepdims=True)
            out += idx.size * ((grp - grand) ** 2).ravel()
        return out

    rows = []
    for name, labels in [
        ("Region", m["region_short"]),
        ("Spaceflight", m["spaceflight"]),
        ("Treatment", m["treatment"]),
        ("Section (rep)", m["rep"].astype(str)),
    ]:
        frac = np.divide(
            ss_factor(labels), ss_total, out=np.full(X.shape[0], np.nan),
            where=ss_total > 0,
        )
        rows.append(
            {
                "factor": name,
                "n_levels": int(labels.nunique()),
                "median_frac_variance": float(np.nanmedian(frac)),
                "mean_frac_variance": float(np.nanmean(frac)),
                "q75_frac_variance": float(np.nanpercentile(frac, 75)),
                "q95_frac_variance": float(np.nanpercentile(frac, 95)),
            }
        )
    return pd.DataFrame(rows)


def figures(coords: pd.DataFrame, var: pd.DataFrame, mat: pd.DataFrame,
            meta: pd.DataFrame, vp: pd.DataFrame) -> None:
    P.use_style()
    import matplotlib.pyplot as plt

    markers = {"GC_SAL": "o", "GC_BuOE": "s", "FLT_SAL": "^", "FLT_BuOE": "D"}

    def pv(scope, pc):
        v = var[(var["scope"] == scope) & (var["PC"] == pc)]["variance_explained"]
        return float(v.iloc[0]) * 100 if len(v) else np.nan

    # --- combined PCA: region + group ---------------------------------------- #
    comb = coords[coords["scope"] == "combined"]
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.9))
    ax = axes[0]
    for region in C.REGION_ORDER:
        s = comb[comb["region_short"] == region]
        ax.scatter(s["PC1"], s["PC2"], s=32, color=P.REGION_COLORS[region],
                   lw=0.4, edgecolor="white", label=region)
    ax.set_xlabel(f"PC1 ({pv('combined','PC1'):.1f}%)")
    ax.set_ylabel(f"PC2 ({pv('combined','PC2'):.1f}%)")
    ax.legend(title="Region", loc="best", fontsize=7, title_fontsize=7)
    ax.set_title("All 48 ROIs, coloured by region", pad=4)

    ax = axes[1]
    for g in C.GROUPS:
        s = comb[comb["group"] == g]
        ax.scatter(s["PC1"], s["PC2"], s=32, color=C.GROUP_COLORS[g],
                   marker=markers[g], lw=0.4, edgecolor="white", label=C.GROUP_LABELS[g])
    ax.set_xlabel(f"PC1 ({pv('combined','PC1'):.1f}%)")
    ax.set_ylabel(f"PC2 ({pv('combined','PC2'):.1f}%)")
    ax.legend(loc="best", fontsize=6.8)
    ax.set_title("Same ROIs, coloured by treatment group", pad=4)
    fig.tight_layout()
    P.save(fig, "PCA_combined")

    # --- per-region PCA ------------------------------------------------------- #
    fig, axes = plt.subplots(1, 4, figsize=(12.0, 3.2))
    for ax, region in zip(axes, C.REGION_ORDER):
        s = coords[coords["scope"] == region]
        for g in C.GROUPS:
            ss = s[s["group"] == g]
            ax.scatter(ss["PC1"], ss["PC2"], s=38, color=C.GROUP_COLORS[g],
                       marker=markers[g], lw=0.4, edgecolor="white",
                       label=C.GROUP_LABELS[g] if region == C.REGION_ORDER[0] else None)
        ax.set_xlabel(f"PC1 ({pv(region,'PC1'):.1f}%)")
        ax.set_ylabel(f"PC2 ({pv(region,'PC2'):.1f}%)")
        ax.set_title(region, pad=4)
    axes[0].legend(loc="best", fontsize=6.2)
    fig.tight_layout()
    P.save(fig, "PCA_per_region")

    # --- correlation heat map + dendrogram ------------------------------------ #
    sub = top_variable(mat)
    corr = np.corrcoef(sub.T.to_numpy())
    m = meta.set_index("sample").loc[sub.columns]
    dist = squareform(np.clip(1 - corr, 0, None), checks=False)
    link = hierarchy.linkage(dist, method="average")

    fig, axes = plt.subplots(
        1, 2, figsize=(11.4, 4.6), gridspec_kw={"width_ratios": [1.25, 1]}
    )
    ax = axes[0]
    dn = hierarchy.dendrogram(link, ax=ax, labels=list(m["title"]),
                              color_threshold=0, above_threshold_color="#555555",
                              leaf_rotation=90)
    ax.set_ylabel("1 $-$ Pearson r")
    for lbl in ax.get_xmajorticklabels():
        region = lbl.get_text().split("_")[0]
        key = {"CA": "CA1", "DG": "DG", "FCT": "FCtx", "CT": "Ctx"}.get(region, region)
        lbl.set_color(P.REGION_COLORS.get(key, "black"))
        lbl.set_fontsize(4.6)
    ax.set_title("Hierarchical clustering (labels coloured by region)", pad=4)

    ax = axes[1]
    order = np.argsort(
        [C.REGION_ORDER.index(r) * 100 + C.GROUPS.index(g)
         for r, g in zip(m["region_short"], m["group"])]
    )
    im = ax.imshow(corr[np.ix_(order, order)], cmap="RdYlBu_r", vmin=np.nanpercentile(corr, 2),
                   vmax=1.0, interpolation="nearest")
    bounds = np.cumsum([12, 12, 12, 12])
    for b in bounds[:-1]:
        ax.axhline(b - 0.5, color="black", lw=0.7)
        ax.axvline(b - 0.5, color="black", lw=0.7)
    ax.set_xticks(bounds - 6.5)
    ax.set_xticklabels(C.REGION_ORDER)
    ax.set_yticks(bounds - 6.5)
    ax.set_yticklabels(C.REGION_ORDER)
    ax.set_title("ROI $\\times$ ROI correlation", pad=4)
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label("Pearson r", fontsize=7.5)
    cb.ax.tick_params(labelsize=7)
    fig.tight_layout()
    P.save(fig, "Clustering_and_correlation")

    # --- variance partition --------------------------------------------------- #
    fig, ax = plt.subplots(figsize=(4.3, 3.0))
    vp2 = vp.sort_values("median_frac_variance", ascending=True)
    colors = {"Region": "#7F8C8D", "Spaceflight": P.FACTOR_COLORS["Spaceflight"],
              "Treatment": P.FACTOR_COLORS["Treatment"], "Section (rep)": "#95A5A6"}
    ax.barh(vp2["factor"], vp2["median_frac_variance"] * 100,
            color=[colors.get(f, "#888888") for f in vp2["factor"]], height=0.6)
    for y, (v, q) in enumerate(zip(vp2["median_frac_variance"], vp2["q95_frac_variance"])):
        ax.text(v * 100 + 0.6, y, f"{v*100:.1f}%", va="center", fontsize=7)
    ax.set_xlabel("Median % of per-target variance explained")
    ax.grid(axis="x", ls=":", lw=0.5, alpha=0.4)
    fig.tight_layout()
    P.save(fig, "Variance_partition")


def main() -> int:
    C.ensure_dirs()
    samples = C.load_sample_table()

    norm = normalisation_check(samples)
    C.write_csv(norm, "normalisation_check.csv")
    print("\nQ3 comparability across the four deposited files:")
    print(norm.to_string(index=False))

    mat, meta = load_all(samples)
    print(f"\nCombined matrix: {mat.shape[0]:,} shared targets x {mat.shape[1]} ROIs")

    all_coords, all_var = [], []
    coords, var = run_pca(mat, meta, "combined")
    all_coords.append(coords)
    all_var.append(var)
    for study in C.STUDIES:
        cols = meta.loc[meta["region_code"] == study.code, "sample"]
        c, v = run_pca(mat[cols.tolist()], meta, study.short)
        all_coords.append(c)
        all_var.append(v)
    coords = pd.concat(all_coords, ignore_index=True)
    var = pd.concat(all_var, ignore_index=True)

    vp = variance_partition(mat, meta)

    C.write_csv(coords, "PCA_coordinates.csv")
    C.write_csv(var, "PCA_variance.csv")
    C.write_csv(vp, "variance_partition.csv")

    sub = top_variable(mat)
    corr = pd.DataFrame(np.corrcoef(sub.T.to_numpy()), index=sub.columns, columns=sub.columns)
    C.write_csv(corr.reset_index(names="sample"), "sample_correlation.csv")

    figures(coords, var, mat, meta, vp)

    print("\nVariance partition (median fraction of per-target variance):")
    print(vp.to_string(index=False))
    print("\nCombined PCA variance explained:")
    print(var[var["scope"] == "combined"].to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
