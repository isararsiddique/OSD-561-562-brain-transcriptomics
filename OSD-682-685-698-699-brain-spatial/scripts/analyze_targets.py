#!/usr/bin/env python3
"""Deep analysis of the CNS / EV target panel across the four brain regions.

The panel is a pre-specified hypothesis, which changes the statistics in a way
that matters here.  The transcriptome-wide analysis in analyze_de.py has to
correct across 15,782 targets and finds nothing; restricting to the 111
measurable panel members reduces the multiple-testing burden by roughly 140-fold,
so the same data can support conclusions that a genome-wide scan cannot.  That is
only legitimate because the panel was fixed in advance, from the project's target
spreadsheet, and not chosen from these results.

Four analyses are run.

  Per-target      Panel-restricted Benjamini-Hochberg FDR, computed within the
                  panel for each region x contrast.
  Competitive     Is the panel as a whole shifted relative to the rest of the
                  transcriptome?  Reported both as a rank test and as a
                  label-permutation test, the latter respecting the correlation
                  between targets and the shared-section design.
  Category        The same question per functional category (neuroinflammation,
                  blood-brain barrier, inflammasome, and so on).
  Attenuation     Whether BuOE damps the panel's spaceflight response, using the
                  noise-corrected estimator from analyze_de.py restricted to the
                  panel, with each arm's own variance.

Outputs
  results/tables/csv/target_stats_all.csv
  results/tables/csv/target_significant_panelFDR.csv
  results/tables/csv/target_setenrichment.csv
  results/tables/csv/target_category_summary.csv
  results/tables/csv/target_attenuation.csv
  results/tables/csv/target_recurrence.csv
  results/figures/Targets_*.{pdf,png}

Usage:  python scripts/analyze_targets.py   (after map_targets.py and analyze_de.py)
"""
from __future__ import annotations

import sys

import numpy as np
import pandas as pd
from scipy import stats

import common as C
import plotting as P

NEG_PROBE = "NegProbe-WTX"
N_PERM = 500
PERM_SEED = 4242

KEY_CONTRASTS = {
    "FLT_vs_GC_in_SAL": "Flight vs Ground (saline)",
    "FLT_vs_GC_in_BuOE": "Flight vs Ground (BuOE)",
    "BuOE_vs_SAL_in_GC": "BuOE vs saline (ground)",
    "BuOE_vs_SAL_in_FLT": "BuOE vs saline (flight)",
    "spaceflight": "Spaceflight main effect",
    "treatment": "Treatment main effect",
    "interaction": "Flight x BuOE interaction",
}
CATEGORY_ORDER = [
    "Neuroinflammation / cytokine",
    "Inflammasome / pyroptosis",
    "Blood-brain barrier / endothelial",
    "Coagulation / vascular",
    "Neuronal / synaptic",
    "Glial / myelin",
    "Injury biomarker",
    "Oxidative stress / metabolism",
    "Transcription / chromatin / other",
]
CATEGORY_COLORS = {
    "Neuroinflammation / cytokine": "#C44E52",
    "Inflammasome / pyroptosis": "#E17C5B",
    "Blood-brain barrier / endothelial": "#4C72B0",
    "Coagulation / vascular": "#6BA3D6",
    "Neuronal / synaptic": "#55A868",
    "Glial / myelin": "#8172B3",
    "Injury biomarker": "#937860",
    "Oxidative stress / metabolism": "#DA8BC3",
    "Transcription / chromatin / other": "#8C8C8C",
}


def load_mapping() -> pd.DataFrame:
    path = C.CONFIG / "cns_ev_targets_mapped.csv"
    if not path.exists():
        raise SystemExit(
            f"missing {path.relative_to(C.ROOT)} -- run scripts/map_targets.py first"
        )
    m = pd.read_csv(path)
    m = m[m["on_geomx_panel"]].copy()
    m = m.drop_duplicates("geomx_target")
    return m


# --------------------------------------------------------------------------- #
# 1. Per-target statistics with panel-restricted FDR
# --------------------------------------------------------------------------- #
def per_target(long: pd.DataFrame, mapping: pd.DataFrame) -> pd.DataFrame:
    genes = mapping["geomx_target"].tolist()
    sub = long[long["gene"].isin(genes)].copy()
    meta = mapping.set_index("geomx_target")
    sub["analyte"] = sub["gene"].map(meta["analyte"])
    sub["human_symbol"] = sub["gene"].map(meta["human_symbol_ensembl"])
    sub["category"] = sub["gene"].map(meta["category"])
    sub["contrast_label"] = sub["contrast"].map(KEY_CONTRASTS)

    # BH within the panel, separately for each region x contrast.
    sub["padj_panel"] = np.nan
    for (region, contrast), idx in sub.groupby(["region_short", "contrast"]).groups.items():
        sub.loc[idx, "padj_panel"] = C.bh_fdr(sub.loc[idx, "pvalue"].to_numpy())
    sub["padj_genomewide"] = sub["padj_BH"]
    sub["t"] = sub["log2FC"] / sub["se"]
    cols = [
        "region_short", "osd", "region", "contrast", "contrast_label", "factor",
        "stratum", "analyte", "human_symbol", "gene", "category", "mean_log2_q3",
        "log2FC", "se", "t", "pvalue", "padj_panel", "padj_genomewide",
    ]
    return sub[cols].sort_values(["contrast", "region_short", "pvalue"])


# --------------------------------------------------------------------------- #
# 2. Competitive set enrichment
# --------------------------------------------------------------------------- #
def _fit(mat: np.ndarray, sf: np.ndarray, trt: np.ndarray, contrast: str):
    res = C.two_way_anova(mat, sf, trt)
    return res[f"log2FC_{contrast}"] / res[f"se_{contrast}"]


def set_enrichment(samples: pd.DataFrame, mapping: pd.DataFrame) -> pd.DataFrame:
    genes = set(mapping["geomx_target"])
    rng = np.random.default_rng(PERM_SEED)
    rows = []

    for study in C.STUDIES:
        log2, sub = C.load_expression(study.code, samples)
        log2 = log2.drop(index=[NEG_PROBE], errors="ignore")
        meta = sub.set_index("sample").loc[log2.columns]
        mat = log2.to_numpy()
        sf = meta["spaceflight"].to_numpy()
        trt = meta["treatment"].to_numpy()
        in_panel = np.array([g in genes for g in log2.index])

        # Permutations: relabel flight/ground within each arm, using one section
        # permutation shared across regions so the design structure is kept.
        perms = []
        for _ in range(N_PERM):
            new_sf = sf.copy()
            for arm in C.TRT_LEVELS:
                idx = np.where(trt == arm)[0]
                new_sf[idx] = rng.permutation(sf[idx])
            perms.append(new_sf)

        for contrast, label in KEY_CONTRASTS.items():
            t = _fit(mat, sf, trt, contrast)
            t_panel, t_bg = t[in_panel], t[~in_panel]

            # Rank test on absolute t: is the panel shifted in magnitude?
            u = stats.mannwhitneyu(np.abs(t_panel), np.abs(t_bg),
                                   alternative="two-sided")
            # Directional: mean signed t.
            obs_mean_t = float(np.nanmean(t_panel))
            obs_mean_abs = float(np.nanmean(np.abs(t_panel)))
            obs_delta = obs_mean_abs - float(np.nanmean(np.abs(t_bg)))

            # Permutation null for the panel's mean |t| minus background mean |t|.
            null = np.empty(N_PERM)
            for k, new_sf in enumerate(perms):
                tp = _fit(mat, new_sf, trt, contrast)
                null[k] = float(np.nanmean(np.abs(tp[in_panel]))
                                - np.nanmean(np.abs(tp[~in_panel])))
            perm_p = float((np.abs(null) >= abs(obs_delta)).mean())

            n_sig = int(
                (C.bh_fdr(2 * stats.t.sf(np.abs(t_panel), 8)) < C.FDR_CUTOFF).sum()
            )
            rows.append(
                {
                    "region_short": study.short,
                    "osd": study.osd,
                    "contrast": contrast,
                    "contrast_label": label,
                    "n_panel_targets": int(in_panel.sum()),
                    "mean_t_panel": obs_mean_t,
                    "mean_abs_t_panel": obs_mean_abs,
                    "mean_abs_t_background": float(np.nanmean(np.abs(t_bg))),
                    "delta_mean_abs_t": obs_delta,
                    "mannwhitney_p_abs_t": float(u.pvalue),
                    "permutation_p_delta_abs_t": perm_p,
                    "perm_null_mean": float(null.mean()),
                    "perm_null_sd": float(null.std(ddof=1)),
                    "n_permutations": N_PERM,
                    "n_panel_sig_panelFDR": n_sig,
                }
            )
        print(f"  {study.short}: enrichment done ({N_PERM} permutations x "
              f"{len(KEY_CONTRASTS)} contrasts)")
    out = pd.DataFrame(rows)
    for contrast in KEY_CONTRASTS:
        m = out["contrast"] == contrast
        out.loc[m, "mannwhitney_q_across_regions"] = C.bh_fdr(
            out.loc[m, "mannwhitney_p_abs_t"].to_numpy()
        )
    return out


# --------------------------------------------------------------------------- #
# 3. Category-level aggregation
# --------------------------------------------------------------------------- #
def category_summary(tstats: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (region, contrast, cat), g in tstats.groupby(
        ["region_short", "contrast", "category"], sort=False
    ):
        t = g["t"].to_numpy()
        t = t[np.isfinite(t)]
        if t.size == 0:
            continue
        # One-sample t-test of the category's mean t against 0: a directional
        # test of coordinated movement within the category.
        tt = stats.ttest_1samp(t, 0.0) if t.size > 1 else None
        rows.append(
            {
                "region_short": region,
                "contrast": contrast,
                "contrast_label": KEY_CONTRASTS.get(contrast, contrast),
                "category": cat,
                "n_targets": int(t.size),
                "mean_t": float(t.mean()),
                "mean_log2FC": float(g["log2FC"].mean()),
                "median_log2FC": float(g["log2FC"].median()),
                "frac_up": float((g["log2FC"] > 0).mean()),
                "p_mean_t_vs_0": float(tt.pvalue) if tt is not None else np.nan,
                "n_nominal_p05": int((g["pvalue"] < 0.05).sum()),
                "n_panelFDR_sig": int((g["padj_panel"] < C.FDR_CUTOFF).sum()),
            }
        )
    out = pd.DataFrame(rows)
    for contrast in out["contrast"].unique():
        m = out["contrast"] == contrast
        out.loc[m, "q_mean_t"] = C.bh_fdr(out.loc[m, "p_mean_t_vs_0"].to_numpy())
    return out


# --------------------------------------------------------------------------- #
# 4. Panel-restricted attenuation
# --------------------------------------------------------------------------- #
def panel_attenuation(samples: pd.DataFrame, mapping: pd.DataFrame) -> pd.DataFrame:
    genes = set(mapping["geomx_target"])
    rows = []
    for study in C.STUDIES:
        log2, sub = C.load_expression(study.code, samples)
        log2 = log2.drop(index=[NEG_PROBE], errors="ignore")
        keep = [g for g in log2.index if g in genes]
        meta = sub.set_index("sample").loc[log2.columns]
        res = C.two_way_anova(
            log2.loc[keep].to_numpy(), meta["spaceflight"].to_numpy(),
            meta["treatment"].to_numpy()
        )
        # Arm-specific noise, as in analyze_de.py.
        se = {}
        for arm in C.TRT_LEVELS:
            ss = np.zeros(len(keep))
            df = 0
            for sf in C.SF_LEVELS:
                cols = meta.index[(meta["treatment"] == arm) & (meta["spaceflight"] == sf)]
                block = log2.loc[keep, list(cols)].to_numpy()
                ss += ((block - block.mean(axis=1, keepdims=True)) ** 2).sum(axis=1)
                df += block.shape[1] - 1
            se[arm] = np.sqrt((ss / df) * (2 / 3))

        lfc_s = res["log2FC_FLT_vs_GC_in_SAL"]
        lfc_b = res["log2FC_FLT_vs_GC_in_BuOE"]
        msq_s = float(np.nanmean(lfc_s ** 2 - se["SAL"] ** 2))
        msq_b = float(np.nanmean(lfc_b ** 2 - se["BuOE"] ** 2))

        rng = np.random.default_rng(PERM_SEED + 1)
        n = len(keep)
        boot = np.empty(4000)
        for i in range(4000):
            idx = rng.integers(0, n, n)
            boot[i] = (np.nanmean(lfc_s[idx] ** 2 - se["SAL"][idx] ** 2)
                       - np.nanmean(lfc_b[idx] ** 2 - se["BuOE"][idx] ** 2))
        lo, hi = np.nanpercentile(boot, [2.5, 97.5])
        frac = float(np.mean(boot <= 0))
        boot_p = float(min(1.0, 2 * min(frac, 1 - frac)))
        rows.append(
            {
                "region_short": study.short,
                "n_panel_targets": n,
                "rms_true_saline": float(np.sqrt(msq_s)) if msq_s > 0 else 0.0,
                "rms_true_BuOE": float(np.sqrt(msq_b)) if msq_b > 0 else 0.0,
                "msq_diff_saline_minus_BuOE": msq_s - msq_b,
                "msq_diff_ci_lo": float(lo),
                "msq_diff_ci_hi": float(hi),
                "bootstrap_p": boot_p,
                "verdict": (
                    "attenuation" if lo > 0 else
                    "amplification" if hi < 0 else "inconclusive"
                ),
            }
        )
    return pd.DataFrame(rows)


def recurrence(tstats: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for contrast in KEY_CONTRASTS:
        g = tstats[(tstats["contrast"] == contrast) & (tstats["pvalue"] < 0.05)]
        if g.empty:
            continue
        agg = g.groupby(["analyte", "gene", "category"]).agg(
            n_regions=("region_short", "nunique"),
            regions=("region_short", lambda s: ",".join(
                sorted(set(s), key=C.REGION_ORDER.index))),
            mean_log2FC=("log2FC", "mean"),
            min_pvalue=("pvalue", "min"),
            min_padj_panel=("padj_panel", "min"),
            consistent=("log2FC", lambda s: bool((s > 0).all() or (s < 0).all())),
        ).reset_index()
        agg["contrast"] = contrast
        agg["contrast_label"] = KEY_CONTRASTS[contrast]
        rows.append(agg[agg["n_regions"] >= 2])
    if not rows:
        return pd.DataFrame()
    return pd.concat(rows, ignore_index=True).sort_values(
        ["contrast", "n_regions", "mean_log2FC"], ascending=[True, False, False]
    )


# --------------------------------------------------------------------------- #
# Figures
# --------------------------------------------------------------------------- #
def figures(mapping_full: pd.DataFrame, tstats: pd.DataFrame, enr: pd.DataFrame,
            cats: pd.DataFrame, att: pd.DataFrame, samples: pd.DataFrame) -> None:
    P.use_style()
    import matplotlib.pyplot as plt

    # --- panel coverage ------------------------------------------------------- #
    cov = (
        mapping_full.groupby("category")
        .agg(n=("analyte", "size"), on_panel=("on_geomx_panel", "sum"))
        .reindex(CATEGORY_ORDER)
        .fillna(0)
    )
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.6),
                             gridspec_kw={"width_ratios": [1.25, 1]})
    ax = axes[0]
    y = np.arange(len(cov))
    ax.barh(y, cov["n"], color="#DDDDDD", height=0.68, label="in target list")
    ax.barh(y, cov["on_panel"], color=[CATEGORY_COLORS[c] for c in cov.index],
            height=0.68, label="measured by GeoMx")
    ax.set_yticks(y)
    ax.set_yticklabels(cov.index, fontsize=6.4)
    ax.invert_yaxis()
    ax.set_xlabel("Analytes")
    ax.legend(fontsize=6.4, loc="lower right")
    ax.grid(axis="x", ls=":", lw=0.5, alpha=0.4)
    for i, (a, b) in enumerate(zip(cov["n"], cov["on_panel"])):
        ax.text(a + 0.4, i, f"{int(b)}/{int(a)}", va="center", fontsize=6.0)
    ax.set_title("CNS / EV target panel coverage", pad=3)
    P.panel_label(ax, "a", dx=-0.42)

    ax = axes[1]
    routes = mapping_full["mapping_route"].value_counts()
    nice = {
        "ensembl_ortholog": "Ensembl ortholog,\non panel",
        "not_on_geomx_panel": "ortholog found,\nnot on panel",
        "no_mouse_ortholog": "no mouse\northolog",
        "symbol_titlecase_fallback": "symbol fallback,\non panel",
    }
    labels = [nice.get(k, k) for k in routes.index]
    ax.bar(range(len(routes)), routes.to_numpy(),
           color=["#55A868", "#BDC3C7", "#C44E52", "#DD8452"][: len(routes)], lw=0)
    ax.set_xticks(range(len(routes)))
    ax.set_xticklabels(labels, fontsize=5.8)
    ax.set_ylabel("Analytes")
    for i, v in enumerate(routes.to_numpy()):
        ax.text(i, v + 1, str(v), ha="center", fontsize=6.4)
    ax.grid(axis="y", ls=":", lw=0.5, alpha=0.4)
    ax.set_title("Human $\\rightarrow$ mouse mapping route", pad=3)
    P.panel_label(ax, "b")
    fig.tight_layout()
    P.save(fig, "Targets_panel_coverage")

    # --- panel-restricted significance --------------------------------------- #
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 3.7),
                             gridspec_kw={"width_ratios": [1, 1.15]})
    ax = axes[0]
    piv = tstats.pivot_table(index="contrast", columns="region_short",
                             values="padj_panel",
                             aggfunc=lambda s: int((s < C.FDR_CUTOFF).sum()))
    piv = piv.reindex(index=list(KEY_CONTRASTS), columns=C.REGION_ORDER).fillna(0)
    im = ax.imshow(piv.to_numpy(), cmap="YlOrRd", vmin=0,
                   vmax=max(1, np.nanmax(piv.to_numpy())))
    ax.set_xticks(range(len(C.REGION_ORDER)))
    ax.set_xticklabels(C.REGION_ORDER)
    ax.set_yticks(range(len(piv)))
    ax.set_yticklabels([KEY_CONTRASTS[c] for c in piv.index], fontsize=6.2)
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            ax.text(j, i, f"{int(piv.to_numpy()[i, j])}", ha="center", va="center",
                    fontsize=7.5)
    cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
    cb.ax.tick_params(labelsize=6.5)
    ax.set_title("Panel targets at BH<0.05\nwithin the panel", pad=3)
    P.panel_label(ax, "a", dx=-0.75)

    ax = axes[1]
    key = "FLT_vs_GC_in_SAL"
    g = tstats[tstats["contrast"] == key].copy()
    for region in C.REGION_ORDER:
        s = g[g["region_short"] == region]
        ax.scatter(s["log2FC"], -np.log10(s["pvalue"]), s=16,
                   color=P.REGION_COLORS[region], lw=0.3, edgecolor="white",
                   label=region, alpha=0.9)
    sig = g[g["padj_panel"] < C.FDR_CUTOFF]
    if len(sig):
        ax.scatter(sig["log2FC"], -np.log10(sig["pvalue"]), s=52, facecolors="none",
                   edgecolors="black", linewidths=0.8, zorder=5)
        texts = []
        for r in sig.itertuples():
            texts.append(ax.text(r.log2FC, -np.log10(r.pvalue), r.analyte,
                                 fontsize=5.8, style="italic"))
        try:
            from adjustText import adjust_text
            adjust_text(texts, ax=ax,
                        arrowprops=dict(arrowstyle="-", color="#666", lw=0.4))
        except Exception:
            pass
    ax.axvline(0, color="#999", lw=0.5, ls=":")
    ax.set_xlabel("log$_2$ fold change")
    ax.set_ylabel("$-$log$_{10}$ nominal $P$")
    ax.legend(fontsize=6.2, title="Region", title_fontsize=6.2)
    ax.set_title("Flight vs Ground (saline), panel only\n"
                 "(ringed = BH<0.05 within panel)", pad=3)
    P.panel_label(ax, "b")
    fig.tight_layout()
    P.save(fig, "Targets_panel_significance")

    # --- category x region heat maps ----------------------------------------- #
    show = ["FLT_vs_GC_in_SAL", "FLT_vs_GC_in_BuOE", "interaction"]
    fig, axes = plt.subplots(1, len(show), figsize=(11.6, 3.9), sharey=True)
    vmax = float(np.nanpercentile(np.abs(cats["mean_t"]), 98)) or 1.0
    im = None
    for ax, contrast in zip(axes, show):
        sub = cats[cats["contrast"] == contrast]
        piv = sub.pivot_table(index="category", columns="region_short", values="mean_t")
        piv = piv.reindex(index=CATEGORY_ORDER, columns=C.REGION_ORDER)
        pq = sub.pivot_table(index="category", columns="region_short",
                             values="p_mean_t_vs_0")
        pq = pq.reindex(index=CATEGORY_ORDER, columns=C.REGION_ORDER)
        im = ax.imshow(piv.to_numpy(), cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                       aspect="auto", interpolation="nearest")
        for i in range(piv.shape[0]):
            for j in range(piv.shape[1]):
                p = pq.to_numpy()[i, j]
                if np.isfinite(p) and p < 0.05:
                    ax.text(j, i, "*", ha="center", va="center", fontsize=8)
        ax.set_xticks(range(len(C.REGION_ORDER)))
        ax.set_xticklabels(C.REGION_ORDER)
        ax.set_title(KEY_CONTRASTS[contrast], fontsize=8, pad=3)
    axes[0].set_yticks(range(len(CATEGORY_ORDER)))
    axes[0].set_yticklabels(CATEGORY_ORDER, fontsize=6.4)
    cb = fig.colorbar(im, ax=axes, fraction=0.016, pad=0.012)
    cb.set_label("Mean $t$ across category members", fontsize=7)
    cb.ax.tick_params(labelsize=6.5)
    fig.suptitle("* marks a category whose mean $t$ differs from zero at $P$<0.05",
                 fontsize=7.6, y=1.01)
    P.save(fig, "Targets_category_heatmap")

    # --- forest plot: the key contrast, grouped by category ------------------ #
    key = "FLT_vs_GC_in_SAL"
    g = tstats[tstats["contrast"] == key]
    piv = g.pivot_table(index=["category", "analyte"], columns="region_short",
                        values="log2FC")
    piv = piv.reindex(columns=C.REGION_ORDER)
    ordered = [i for c in CATEGORY_ORDER for i in piv.index if i[0] == c]
    piv = piv.loc[ordered]
    fig, ax = plt.subplots(figsize=(6.4, 0.145 * len(piv) + 1.4))
    y = np.arange(len(piv))
    for j, region in enumerate(C.REGION_ORDER):
        ax.scatter(piv[region], y + (j - 1.5) * 0.17, s=7,
                   color=P.REGION_COLORS[region], lw=0, label=region)
    ax.axvline(0, color="black", lw=0.7)
    bounds, last, start = [], None, 0
    for i, (cat, _) in enumerate(piv.index):
        if cat != last and last is not None:
            bounds.append((start, i - 1, last))
            start = i
        last = cat
    bounds.append((start, len(piv) - 1, last))
    for a, b, cat in bounds:
        ax.axhspan(a - 0.5, b + 0.5, color=CATEGORY_COLORS[cat], alpha=0.07, lw=0)
        ax.text(1.02, (a + b) / 2, cat, transform=ax.get_yaxis_transform(),
                fontsize=5.4, va="center", ha="left", color=CATEGORY_COLORS[cat])
    ax.set_yticks(y)
    ax.set_yticklabels([a for _, a in piv.index], fontsize=4.6)
    ax.invert_yaxis()
    ax.set_xlabel("log$_2$ fold change, Flight vs Ground (saline)")
    ax.legend(fontsize=6, ncols=4, loc="lower center", bbox_to_anchor=(0.5, 1.0))
    ax.grid(axis="x", ls=":", lw=0.5, alpha=0.4)
    fig.tight_layout()
    P.save(fig, "Targets_forest_flight_saline")

    # --- set enrichment ------------------------------------------------------ #
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 3.6))
    ax = axes[0]
    x = np.arange(len(KEY_CONTRASTS))
    for ri, region in enumerate(C.REGION_ORDER):
        vals = [
            float(enr[(enr["region_short"] == region) & (enr["contrast"] == c)]
                  ["delta_mean_abs_t"].iloc[0])
            for c in KEY_CONTRASTS
        ]
        ax.bar(x + (ri - 1.5) * 0.2, vals, width=0.2,
               color=P.REGION_COLORS[region], label=region, lw=0)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([KEY_CONTRASTS[c] for c in KEY_CONTRASTS], fontsize=5.8,
                       rotation=28, ha="right")
    ax.set_ylabel("Panel mean |$t$| $-$ background mean |$t$|")
    ax.legend(fontsize=6.2, title="Region", title_fontsize=6.2, ncols=2)
    ax.grid(axis="y", ls=":", lw=0.5, alpha=0.4)
    ax.set_title("Is the panel more responsive than the transcriptome?", pad=3)
    P.panel_label(ax, "a")

    ax = axes[1]
    piv = enr.pivot_table(index="contrast", columns="region_short",
                          values="permutation_p_delta_abs_t")
    piv = piv.reindex(index=list(KEY_CONTRASTS), columns=C.REGION_ORDER)
    im = ax.imshow(piv.to_numpy(), cmap="viridis_r", vmin=0, vmax=1)
    ax.set_xticks(range(len(C.REGION_ORDER)))
    ax.set_xticklabels(C.REGION_ORDER)
    ax.set_yticks(range(len(piv)))
    ax.set_yticklabels([KEY_CONTRASTS[c] for c in piv.index], fontsize=6.2)
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            v = piv.to_numpy()[i, j]
            ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=6.5,
                    color="white" if v < 0.5 else "black")
    cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
    cb.set_label("permutation $P$", fontsize=7)
    cb.ax.tick_params(labelsize=6.5)
    ax.set_title(f"Label-permutation test ({N_PERM} permutations)", pad=3)
    P.panel_label(ax, "b", dx=-0.75)
    fig.tight_layout()
    P.save(fig, "Targets_set_enrichment")

    # --- attenuation restricted to the panel --------------------------------- #
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.4))
    r = att.set_index("region_short").reindex(C.REGION_ORDER)
    x = np.arange(len(C.REGION_ORDER))
    ax = axes[0]
    ax.bar(x - 0.2, r["rms_true_saline"], width=0.4, color="#C44E52",
           label="saline", lw=0)
    ax.bar(x + 0.2, r["rms_true_BuOE"], width=0.4, color="#4C72B0",
           label="BuOE", lw=0)
    ax.set_xticks(x)
    ax.set_xticklabels(C.REGION_ORDER)
    ax.set_ylabel("Noise-corrected RMS spaceflight\neffect, panel targets (log$_2$)")
    ax.legend(fontsize=6.5)
    ax.grid(axis="y", ls=":", lw=0.5, alpha=0.4)
    ax.set_title("Panel response magnitude", pad=3)
    P.panel_label(ax, "a")

    ax = axes[1]
    est = r["msq_diff_saline_minus_BuOE"].to_numpy()
    lo = est - r["msq_diff_ci_lo"].to_numpy()
    hi = r["msq_diff_ci_hi"].to_numpy() - est
    ax.errorbar(est, x, xerr=[lo, hi], fmt="o", color="#2C3E50", ms=5, capsize=2.5,
                lw=1.0)
    ax.axvline(0, color="black", ls="--", lw=0.9)
    ax.set_yticks(x)
    ax.set_yticklabels(C.REGION_ORDER)
    ax.set_xlabel("saline $-$ BuOE mean squared effect\n(>0 = attenuation)")
    span = float(np.nanmax(r["msq_diff_ci_hi"]) - np.nanmin(r["msq_diff_ci_lo"]))
    ax.set_xlim(float(np.nanmin(r["msq_diff_ci_lo"])) - 0.1 * span,
                float(np.nanmax(r["msq_diff_ci_hi"])) + 0.4 * span)
    ax.set_ylim(-0.6, len(C.REGION_ORDER) - 0.4)
    for i, p in enumerate(r["bootstrap_p"]):
        txt = "$p$<0.001" if p < 0.001 else f"$p$={p:.3f}"
        ax.text(0.995, i + 0.3, txt, transform=ax.get_yaxis_transform(),
                ha="right", va="center", fontsize=6.0)
    ax.grid(axis="x", ls=":", lw=0.5, alpha=0.4)
    ax.set_title("BuOE attenuation, panel only", pad=3)
    P.panel_label(ax, "b")
    fig.tight_layout()
    P.save(fig, "Targets_attenuation")

    # --- expression of the strongest panel targets --------------------------- #
    key = "FLT_vs_GC_in_SAL"
    g = tstats[tstats["contrast"] == key]
    best = (
        g.groupby(["analyte", "gene"])["pvalue"].min().sort_values().head(8).index.tolist()
    )
    frames = []
    for study in C.STUDIES:
        log2, sub = C.load_expression(study.code, samples)
        keep = [gene for _, gene in best if gene in log2.index]
        sm = log2.loc[keep].T.reset_index(names="sample")
        sm = sub[["sample", "region_short", "group"]].merge(sm, on="sample")
        frames.append(sm.melt(id_vars=["sample", "region_short", "group"],
                              var_name="gene", value_name="log2_q3"))
    expr = pd.concat(frames, ignore_index=True)
    name_of = {gene: analyte for analyte, gene in best}

    fig, axes = plt.subplots(2, 4, figsize=(11.0, 4.6), squeeze=False)
    for k, (analyte, gene) in enumerate(best):
        ax = axes[k // 4][k % 4]
        gg = expr[expr["gene"] == gene]
        for ri, region in enumerate(C.REGION_ORDER):
            for gi, grp in enumerate(C.GROUPS):
                vals = gg[(gg["region_short"] == region) & (gg["group"] == grp)]["log2_q3"]
                xs = np.full(len(vals), ri + (gi - 1.5) * 0.19)
                ax.scatter(xs, vals, s=9, color=C.GROUP_COLORS[grp], lw=0.2,
                           edgecolor="white")
                if len(vals):
                    ax.plot([ri + (gi - 1.5) * 0.19 - 0.07, ri + (gi - 1.5) * 0.19 + 0.07],
                            [vals.mean()] * 2, color=C.GROUP_COLORS[grp], lw=1.1)
        ax.set_xticks(range(len(C.REGION_ORDER)))
        ax.set_xticklabels(C.REGION_ORDER, fontsize=5.8)
        title = analyte if analyte == gene else f"{analyte} ({gene})"
        ax.set_title(title, fontsize=7.2, style="italic", pad=3)
        ax.tick_params(labelsize=6)
        if k % 4 == 0:
            ax.set_ylabel("log$_2$ Q3", fontsize=7)
    handles = [
        plt.Line2D([], [], marker="o", ls="", color=C.GROUP_COLORS[g],
                   label=C.GROUP_LABELS[g], ms=4)
        for g in C.GROUPS
    ]
    fig.legend(handles=handles, loc="lower center", ncols=4, fontsize=6.5,
               bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("Panel targets with the strongest Flight-vs-Ground signal (saline)",
                 fontsize=8.5, y=1.0)
    fig.tight_layout()
    P.save(fig, "Targets_top_expression")


def main() -> int:
    C.ensure_dirs()
    mapping_full = pd.read_csv(C.CONFIG / "cns_ev_targets_mapped.csv")
    mapping = load_mapping()
    print(f"Target panel: {len(mapping)} measurable targets of {len(mapping_full)} analytes")

    de_path = C.CSV_DIR / "DE_all_contrasts_long.csv"
    if not de_path.exists():
        print(f"ERROR: {de_path} missing -- run scripts/analyze_de.py first",
              file=sys.stderr)
        return 1
    long = pd.read_csv(de_path)
    samples = C.load_sample_table()

    tstats = per_target(long, mapping)
    C.write_csv(tstats, "target_stats_all.csv")

    sig = tstats[tstats["padj_panel"] < C.FDR_CUTOFF].sort_values("padj_panel")
    C.write_csv(sig, "target_significant_panelFDR.csv")

    print("\nRunning label-permutation set enrichment ...")
    enr = set_enrichment(samples, mapping)
    C.write_csv(enr, "target_setenrichment.csv")

    cats = category_summary(tstats)
    C.write_csv(cats, "target_category_summary.csv")

    att = panel_attenuation(samples, mapping)
    C.write_csv(att, "target_attenuation.csv")

    rec = recurrence(tstats)
    if len(rec):
        C.write_csv(rec, "target_recurrence.csv")

    figures(mapping_full, tstats, enr, cats, att, samples)

    pd.set_option("display.width", 220)
    pd.set_option("display.max_columns", 60)
    print(f"\n=== Panel-restricted FDR: {len(sig)} significant target x region x contrast "
          f"results ===")
    if len(sig):
        print(
            sig[["analyte", "gene", "category", "region_short", "contrast_label",
                 "log2FC", "pvalue", "padj_panel", "padj_genomewide"]]
            .head(30).to_string(index=False)
        )

    print("\n=== Set enrichment (panel vs transcriptome) ===")
    print(
        enr[["region_short", "contrast", "delta_mean_abs_t", "mannwhitney_p_abs_t",
             "permutation_p_delta_abs_t", "n_panel_sig_panelFDR"]].to_string(index=False)
    )

    print("\n=== Categories with mean t different from zero (P<0.05) ===")
    hot = cats[cats["p_mean_t_vs_0"] < 0.05].sort_values("p_mean_t_vs_0")
    print(
        hot[["region_short", "contrast_label", "category", "n_targets", "mean_t",
             "mean_log2FC", "p_mean_t_vs_0", "q_mean_t"]].head(25).to_string(index=False)
    )

    print("\n=== BuOE attenuation, panel targets only ===")
    print(att.to_string(index=False))

    if len(rec):
        print("\n=== Panel targets nominal in >=2 regions, Flight vs Ground (saline) ===")
        r = rec[rec["contrast"] == "FLT_vs_GC_in_SAL"]
        print(
            r[["analyte", "gene", "category", "n_regions", "regions", "mean_log2FC",
               "min_pvalue", "consistent"]].to_string(index=False)
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
