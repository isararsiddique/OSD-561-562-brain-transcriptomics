#!/usr/bin/env python3
"""Differential expression across the four regions, and the BuOE attenuation test.

Three things happen here.

1. The statistics deposited with GSE239336 are audited.  The exported
   "Adjusted pvalue" column is rounded to two decimals and is smaller than the
   corresponding raw p-value for roughly three quarters of targets, which no
   multiple-testing procedure can produce.  It is therefore reported for
   provenance but never used to call DEGs.

2. Differential expression is re-derived from the Q3-normalised matrix with a
   per-target 2x2 fixed-effects model (Spaceflight x Treatment, n = 3 per cell),
   giving the two main effects, the interaction, and the four simple contrasts,
   with Benjamini-Hochberg FDR at full precision.  This also supplies the
   Flight-vs-Ground comparison within BuOE, which was never deposited.

3. The paper's central claim -- that BuOE attenuates the spaceflight response --
   is tested.  The formal per-target test is the interaction term.  Because
   n = 3 leaves that test underpowered, a noise-corrected global effect-size
   comparison is also computed: the sampling variance of each fold change is
   already known from the model, so the mean squared *true* effect can be
   estimated as mean(log2FC^2 - SE^2) and compared between strata.  A
   selection-based statistic mirroring the conventional approach is reported
   alongside it and explicitly flagged, because conditioning on the saline
   response guarantees apparent attenuation through regression to the mean.

Outputs
  results/tables/csv/deposited_stats_audit.csv
  results/tables/csv/DE_all_contrasts_long.csv     gene-level, every contrast
  results/tables/csv/DEG_summary_all.csv           counts under each definition
  results/tables/csv/attenuation_summary.csv
  results/tables/csv/interaction_top_targets.csv
  results/tables/csv/deposited_vs_derived.csv
  results/figures/Volcano_*.{pdf,png}, DEG_landscape.{pdf,png},
  results/figures/Attenuation_*.{pdf,png}

Usage:  python scripts/analyze_de.py
"""
from __future__ import annotations

import sys
import zlib

import numpy as np
import pandas as pd
from scipy import stats

import common as C
import plotting as P

NEG_PROBE = "NegProbe-WTX"
BOOT_SEED = 20240808
RNG = np.random.default_rng(BOOT_SEED)
N_BOOT = 10000

# Contrasts reported in the landscape, in display order.
LANDSCAPE = [
    ("FLT_vs_GC_in_SAL", "Flight vs Ground (saline)", "Spaceflight"),
    ("FLT_vs_GC_in_BuOE", "Flight vs Ground (BuOE)", "Spaceflight"),
    ("BuOE_vs_SAL_in_GC", "BuOE vs saline (ground)", "Treatment"),
    ("BuOE_vs_SAL_in_FLT", "BuOE vs saline (flight)", "Treatment"),
    ("spaceflight", "Spaceflight main effect", "Spaceflight"),
    ("treatment", "Treatment main effect", "Treatment"),
    ("interaction", "Flight $\\times$ BuOE interaction", "Interaction"),
]


# --------------------------------------------------------------------------- #
# 1. Audit of the deposited statistics
# --------------------------------------------------------------------------- #
def audit_deposited() -> pd.DataFrame:
    rows = []
    for study in C.STUDIES:
        for contrast, spec in C.PROVIDED_CONTRASTS.items():
            de = C.read_provided_de(study.code, contrast)
            p = de["pvalue"].to_numpy()
            adj = de["padj_deposited"].to_numpy()
            ok = np.isfinite(p) & np.isfinite(adj)
            impossible = ok & (adj < p - 1e-12)
            bh = C.bh_fdr(p)
            rows.append(
                {
                    "region_short": study.short,
                    "osd": study.osd,
                    "contrast": contrast,
                    "contrast_label": spec["label"],
                    "n_targets": int(len(de)),
                    "n_adj_below_raw_p": int(impossible.sum()),
                    "frac_adj_below_raw_p": float(impossible.mean()),
                    "n_unique_adj_values": int(pd.Series(adj).nunique()),
                    "n_adj_exactly_zero": int((adj == 0).sum()),
                    "min_raw_p": float(np.nanmin(p)),
                    "n_sig_deposited_adj_lt_0.05": int((adj < 0.05).sum()),
                    "n_sig_recomputed_BH_lt_0.05": int((bh < 0.05).sum()),
                    "min_recomputed_BH": float(np.nanmin(bh)),
                }
            )
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# 2. Re-derived differential expression
# --------------------------------------------------------------------------- #
def derive_de(samples: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    long_rows = []
    per_region: dict[str, pd.DataFrame] = {}

    for study in C.STUDIES:
        log2, sub = C.load_expression(study.code, samples)
        log2 = log2.drop(index=[NEG_PROBE], errors="ignore")
        meta = sub.set_index("sample").loc[log2.columns]
        res = C.two_way_anova(
            log2.to_numpy(), meta["spaceflight"].to_numpy(), meta["treatment"].to_numpy()
        )

        tab = pd.DataFrame({"gene": log2.index, "mean_log2_q3": res["mean_log2"]})
        names = list(C.DERIVED_CONTRASTS) + ["spaceflight", "treatment", "interaction"]
        for name in names:
            tab[f"log2FC_{name}"] = res[f"log2FC_{name}"]
            tab[f"se_{name}"] = res[f"se_{name}"]
            tab[f"p_{name}"] = res[f"p_{name}"]
            tab[f"padj_{name}"] = C.bh_fdr(res[f"p_{name}"])
        tab["df_resid"] = res["df_resid"]
        tab.insert(0, "region_short", study.short)
        tab.insert(1, "osd", study.osd)
        per_region[study.short] = tab

        # Long form, one row per target x contrast.
        for name in names:
            spec = C.DERIVED_CONTRASTS.get(name)
            long_rows.append(
                pd.DataFrame(
                    {
                        "region_short": study.short,
                        "osd": study.osd,
                        "region": study.region,
                        "contrast": name,
                        "factor": spec[2] if spec else name.capitalize(),
                        "stratum": spec[3] if spec else "all",
                        "numerator": spec[0] if spec else "",
                        "denominator": spec[1] if spec else "",
                        "gene": tab["gene"],
                        "mean_log2_q3": tab["mean_log2_q3"],
                        "log2FC": tab[f"log2FC_{name}"],
                        "se": tab[f"se_{name}"],
                        "pvalue": tab[f"p_{name}"],
                        "padj_BH": tab[f"padj_{name}"],
                    }
                )
            )

    return pd.concat(long_rows, ignore_index=True), per_region


def compare_deposited(per_region: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Confirm the re-derived fold changes reproduce the deposited ones."""
    rows = []
    derived_for = {
        "GCvsFLT-SAL": "FLT_vs_GC_in_SAL",
        "SALvsBuOE-GC": "BuOE_vs_SAL_in_GC",
        "SALvsBuOE-FLT": "BuOE_vs_SAL_in_FLT",
    }
    for study in C.STUDIES:
        tab = per_region[study.short].set_index("gene")
        for contrast, dname in derived_for.items():
            de = C.read_provided_de(study.code, contrast).drop_duplicates("gene")
            de = de.set_index("gene")
            shared = tab.index.intersection(de.index)
            a = tab.loc[shared, f"log2FC_{dname}"].to_numpy()
            b = de.loc[shared, "log2FC_deposited"].to_numpy() * C.sign_for(contrast)
            pa = tab.loc[shared, f"p_{dname}"].to_numpy()
            pb = de.loc[shared, "pvalue"].to_numpy()
            rows.append(
                {
                    "region_short": study.short,
                    "deposited_contrast": contrast,
                    "derived_contrast": dname,
                    "n_shared": len(shared),
                    "log2FC_r": float(np.corrcoef(a, b)[0, 1]),
                    "log2FC_max_abs_diff": float(np.nanmax(np.abs(a - b))),
                    "pvalue_r": float(
                        np.corrcoef(pa[np.isfinite(pa) & np.isfinite(pb)],
                                    pb[np.isfinite(pa) & np.isfinite(pb)])[0, 1]
                    ),
                    "pvalue_max_abs_diff": float(np.nanmax(np.abs(pa - pb))),
                }
            )
    return pd.DataFrame(rows)


def summarise_degs(long: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (region, contrast), g in long.groupby(["region_short", "contrast"], sort=False):
        lfc = g["log2FC"].to_numpy()
        padj = g["padj_BH"].to_numpy()
        p = g["pvalue"].to_numpy()
        rec = {
            "region_short": region,
            "osd": g["osd"].iloc[0],
            "contrast": contrast,
            "factor": g["factor"].iloc[0],
            "stratum": g["stratum"].iloc[0],
            "n_targets": len(g),
            "min_padj_BH": float(np.nanmin(padj)),
            "min_pvalue": float(np.nanmin(p)),
        }
        strict = C.deg_counts(padj, lfc)
        rec.update({f"fdr05_lfc1_{k}": v for k, v in strict.items()})
        relaxed = C.deg_counts(padj, lfc, min_lfc=C.LFC_CUTOFF_RELAXED)
        rec.update({f"fdr05_lfc05_{k}": v for k, v in relaxed.items()})

        # Unfiltered nominal counts are the only ones directly comparable with
        # the number expected under the global null (alpha x n_targets); the
        # fold-change-filtered counts below are not, so they are kept separate.
        for pc in (0.01, 0.05):
            sel_p = np.isfinite(p) & (p < pc)
            rec[f"nominal_p{pc:g}_n"] = int(sel_p.sum())
            rec[f"expected_null_p{pc:g}"] = round(pc * len(g), 1)
            rec[f"obs_over_expected_p{pc:g}"] = round(
                float(sel_p.sum()) / (pc * len(g)), 3
            )
            sel = sel_p & (np.abs(lfc) >= C.LFC_CUTOFF)
            rec[f"nominal_p{pc:g}_lfc1_n"] = int(sel.sum())
            rec[f"nominal_p{pc:g}_lfc1_up"] = int((sel & (lfc > 0)).sum())
            rec[f"nominal_p{pc:g}_lfc1_down"] = int((sel & (lfc < 0)).sum())
        sel = np.isfinite(p) & (p < 0.01) & (np.abs(lfc) >= C.LFC_CUTOFF_RELAXED)
        rec["nominal_p0.01_lfc05_n"] = int(sel.sum())
        rows.append(rec)
    return pd.DataFrame(rows)


def sensitivity(per_region: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Smallest fold change this design could have detected.

    With n = 3 per cell the residual variance is known, so the minimum
    detectable effect at 80% power follows directly from the standard error:
    MDE = (t_{1-alpha/2, df} + t_{0.8, df}) * SE.  Reported at the nominal
    threshold and at a genome-wide threshold (alpha / n_targets), which bounds
    what would have been needed to survive multiple-testing control.
    """
    rows = []
    for study in C.STUDIES:
        tab = per_region[study.short]
        df = float(tab["df_resid"].iloc[0])
        n_targets = len(tab)
        t_power = stats.t.ppf(0.8, df)
        for contrast in ["FLT_vs_GC_in_SAL", "FLT_vs_GC_in_BuOE", "spaceflight", "interaction"]:
            se = tab[f"se_{contrast}"].to_numpy()
            rec = {
                "region_short": study.short,
                "contrast": contrast,
                "df_resid": df,
                "median_se_log2": float(np.nanmedian(se)),
            }
            for label, alpha in [
                ("nominal_a0.05", 0.05),
                ("genomewide_a0.05_over_n", 0.05 / n_targets),
            ]:
                mde = (stats.t.ppf(1 - alpha / 2, df) + t_power) * se
                rec[f"mde_log2FC_{label}"] = float(np.nanmedian(mde))
                rec[f"mde_fold_change_{label}"] = float(2 ** np.nanmedian(mde))
            rows.append(rec)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# 3. Attenuation
# --------------------------------------------------------------------------- #
def _msq_true(lfc: np.ndarray, se: np.ndarray) -> float:
    """Unbiased estimate of the mean squared *true* effect.

    E[lfc^2] = true^2 + SE^2 and SE^2 is known from the model, so the difference
    estimates the mean squared true effect.  Deliberately not clipped at zero:
    a negative value is informative, meaning the observed spread of fold changes
    does not exceed sampling noise.
    """
    return float(np.nanmean(lfc ** 2 - se ** 2))


def _rms_true(lfc: np.ndarray, se: np.ndarray) -> float:
    v = _msq_true(lfc, se)
    return float(np.sqrt(v)) if v > 0 else 0.0


def stratum_noise(samples: pd.DataFrame, code: str) -> pd.DataFrame:
    """Per-target within-cell variance estimated separately in each drug arm.

    The pooled 2x2 model assumes one residual variance for all four cells.  When
    the two arms are compared against each other that assumption is load-bearing:
    if the BuOE cells are noisier, a pooled standard error understates their
    noise and overstates the saline noise, which alone would make the BuOE
    response look larger.  Estimating each arm's variance from its own two cells
    (df = 4) removes that dependence.
    """
    log2, sub = C.load_expression(code, samples)
    log2 = log2.drop(index=[NEG_PROBE], errors="ignore")
    meta = sub.set_index("sample").loc[log2.columns]
    out = {"gene": log2.index.to_numpy()}
    for arm in C.TRT_LEVELS:
        ss = np.zeros(log2.shape[0])
        df = 0
        for sf in C.SF_LEVELS:
            cols = meta.index[(meta["treatment"] == arm) & (meta["spaceflight"] == sf)]
            block = log2[list(cols)].to_numpy()
            ss += ((block - block.mean(axis=1, keepdims=True)) ** 2).sum(axis=1)
            df += block.shape[1] - 1
        var = ss / df
        out[f"var_{arm}"] = var
        out[f"df_{arm}"] = np.full(log2.shape[0], df, dtype=float)
        # SE of the Flight-minus-Ground difference inside this arm (n = 3 each).
        out[f"se_within_{arm}"] = np.sqrt(var * (1 / 3 + 1 / 3))
    return pd.DataFrame(out)


def attenuation(
    per_region: dict[str, pd.DataFrame], samples: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    top_rows = []
    for study in C.STUDIES:
        tab = per_region[study.short]
        noise = stratum_noise(samples, study.code)
        if not np.array_equal(noise["gene"].to_numpy(), tab["gene"].to_numpy()):
            raise ValueError(f"{study.short}: target order mismatch in stratum_noise")

        lfc_s = tab["log2FC_FLT_vs_GC_in_SAL"].to_numpy()
        lfc_b = tab["log2FC_FLT_vs_GC_in_BuOE"].to_numpy()
        p_s = tab["p_FLT_vs_GC_in_SAL"].to_numpy()
        # Primary: each arm's own noise. Pooled values retained as a check.
        se_s = noise["se_within_SAL"].to_numpy()
        se_b = noise["se_within_BuOE"].to_numpy()
        se_s_pooled = tab["se_FLT_vs_GC_in_SAL"].to_numpy()
        se_b_pooled = tab["se_FLT_vs_GC_in_BuOE"].to_numpy()
        var_ratio = float(
            np.nanmedian(noise["var_BuOE"].to_numpy() / noise["var_SAL"].to_numpy())
        )

        msq_s = _msq_true(lfc_s, se_s)
        msq_b = _msq_true(lfc_b, se_b)
        rms_s = _rms_true(lfc_s, se_s)
        rms_b = _rms_true(lfc_b, se_b)

        # Bootstrap over targets. The primary quantity is the *difference* in
        # mean squared true effect (saline - BuOE): attenuation means it is
        # positive. The ratio is unstable whenever the saline term sits at the
        # noise floor, so it is recorded but only trusted when the saline
        # magnitude is itself distinguishable from zero.
        # Region-seeded generator so each region's interval is reproducible
        # independently of how many draws other regions consume. zlib.crc32 is
        # used rather than hash(), which is salted per process.
        rng = np.random.default_rng(BOOT_SEED + zlib.crc32(study.code.encode()))
        n = len(lfc_s)
        boot_s = np.empty(N_BOOT)
        boot_b = np.empty(N_BOOT)
        boot_diff = np.empty(N_BOOT)
        for i in range(N_BOOT):
            idx = rng.integers(0, n, n)
            a = _msq_true(lfc_s[idx], se_s[idx])
            b = _msq_true(lfc_b[idx], se_b[idx])
            boot_s[i] = a
            boot_b[i] = b
            boot_diff[i] = a - b
        s_lo, s_hi = np.nanpercentile(boot_s, [2.5, 97.5])
        b_lo, b_hi = np.nanpercentile(boot_b, [2.5, 97.5])
        d_lo, d_hi = np.nanpercentile(boot_diff, [2.5, 97.5])
        # Two-sided bootstrap tail probability: how often the resampled
        # difference lands on the other side of zero. Reported so borderline
        # cases are visible as a number rather than hidden behind a binary.
        frac_le0 = float(np.mean(boot_diff <= 0))
        boot_p = float(min(1.0, 2 * min(frac_le0, 1 - frac_le0)))
        saline_detectable = bool(s_lo > 0)
        ratio = (rms_b / rms_s) if rms_s > 0 else np.nan
        if saline_detectable:
            verdict = (
                "attenuation supported" if d_lo > 0
                else "amplification supported" if d_hi < 0
                else "inconclusive"
            )
            if verdict != "inconclusive" and boot_p > 0.01:
                verdict += " (borderline)"
        else:
            verdict = "untestable: saline spaceflight effect not distinguishable from noise"

        # Same test under the pooled-variance assumption, to expose how much the
        # conclusion depends on it.
        msq_s_pool = _msq_true(lfc_s, se_s_pooled)
        msq_b_pool = _msq_true(lfc_b, se_b_pooled)
        boot_diff_pool = np.empty(N_BOOT)
        for i in range(N_BOOT):
            idx = RNG.integers(0, n, n)
            boot_diff_pool[i] = _msq_true(lfc_s[idx], se_s_pooled[idx]) - _msq_true(
                lfc_b[idx], se_b_pooled[idx]
            )
        dp_lo, dp_hi = np.nanpercentile(boot_diff_pool, [2.5, 97.5])
        verdict_pool = (
            "attenuation supported" if dp_lo > 0
            else "amplification supported" if dp_hi < 0
            else "inconclusive"
        )

        # Reference slopes. OLS is biased low by noise in the predictor; Deming
        # regression with the known error ratio is not.
        ok = np.isfinite(lfc_s) & np.isfinite(lfc_b)
        ols = float(np.polyfit(lfc_s[ok], lfc_b[ok], 1)[0])
        lam = float(np.nanmean(se_b[ok] ** 2) / np.nanmean(se_s[ok] ** 2))
        vx, vy = np.var(lfc_s[ok], ddof=1), np.var(lfc_b[ok], ddof=1)
        cxy = np.cov(lfc_s[ok], lfc_b[ok], ddof=1)[0, 1]
        deming = float(
            (vy - lam * vx + np.sqrt((vy - lam * vx) ** 2 + 4 * lam * cxy ** 2))
            / (2 * cxy)
        )

        inter_padj = tab["padj_interaction"].to_numpy()
        inter_p = tab["p_interaction"].to_numpy()

        # Selection-based statistic: mirrors the conventional approach and is
        # reported only to show how much of the apparent attenuation it manufactures.
        sel = np.argsort(np.where(np.isfinite(p_s), p_s, np.inf))[:500]
        med_sel_s = float(np.nanmedian(np.abs(lfc_s[sel])))
        med_sel_b = float(np.nanmedian(np.abs(lfc_b[sel])))
        try:
            w = stats.wilcoxon(np.abs(lfc_s[sel]), np.abs(lfc_b[sel]))
            w_p = float(w.pvalue)
        except ValueError:
            w_p = np.nan
        # The same statistic computed on a random target set: any gap here is
        # pure selection artefact.
        rnd = RNG.choice(n, size=500, replace=False)
        med_rnd_s = float(np.nanmedian(np.abs(lfc_s[rnd])))
        med_rnd_b = float(np.nanmedian(np.abs(lfc_b[rnd])))

        rows.append(
            {
                "region_short": study.short,
                "osd": study.osd,
                "verdict": verdict,
                "msq_true_saline": msq_s,
                "msq_true_saline_ci_lo": float(s_lo),
                "msq_true_saline_ci_hi": float(s_hi),
                "saline_effect_above_noise": saline_detectable,
                "msq_true_BuOE": msq_b,
                "msq_true_BuOE_ci_lo": float(b_lo),
                "msq_true_BuOE_ci_hi": float(b_hi),
                "msq_diff_saline_minus_BuOE": msq_s - msq_b,
                "msq_diff_ci_lo": float(d_lo),
                "msq_diff_ci_hi": float(d_hi),
                "msq_diff_bootstrap_p": boot_p,
                "rms_true_effect_saline": rms_s,
                "rms_true_effect_BuOE": rms_b,
                "rms_ratio_BuOE_over_saline": ratio,
                "median_var_ratio_BuOE_over_saline": var_ratio,
                "msq_true_saline_pooled_variance": msq_s_pool,
                "msq_true_BuOE_pooled_variance": msq_b_pool,
                "msq_diff_pooled_variance": msq_s_pool - msq_b_pool,
                "verdict_under_pooled_variance": verdict_pool,
                "verdict_depends_on_variance_assumption": bool(
                    verdict_pool.split(":")[0] != verdict.split(":")[0]
                ),
                "slope_ols_biased": ols,
                "slope_deming": deming,
                "n_interaction_fdr05": int(np.nansum(inter_padj < C.FDR_CUTOFF)),
                "min_interaction_padj": float(np.nanmin(inter_padj)),
                "n_interaction_nominal_p01_lfc1": int(
                    np.nansum(
                        (inter_p < 0.01)
                        & (np.abs(tab["log2FC_interaction"].to_numpy()) >= C.LFC_CUTOFF)
                    )
                ),
                "n_nominal_SF_saline_p01_lfc05": int(
                    np.nansum((p_s < 0.01) & (np.abs(lfc_s) >= C.LFC_CUTOFF_RELAXED))
                ),
                "n_nominal_SF_BuOE_p01_lfc05": int(
                    np.nansum(
                        (tab["p_FLT_vs_GC_in_BuOE"].to_numpy() < 0.01)
                        & (np.abs(lfc_b) >= C.LFC_CUTOFF_RELAXED)
                    )
                ),
                "median_abs_lfc_top500_saline_SELECTION_BIASED": med_sel_s,
                "median_abs_lfc_top500_BuOE_SELECTION_BIASED": med_sel_b,
                "wilcoxon_p_selected_SELECTION_BIASED": w_p,
                "median_abs_lfc_random500_saline": med_rnd_s,
                "median_abs_lfc_random500_BuOE": med_rnd_b,
            }
        )

        order = np.argsort(np.where(np.isfinite(inter_p), inter_p, np.inf))[:25]
        top_rows.append(
            pd.DataFrame(
                {
                    "region_short": study.short,
                    "gene": tab["gene"].to_numpy()[order],
                    "log2FC_interaction": tab["log2FC_interaction"].to_numpy()[order],
                    "p_interaction": inter_p[order],
                    "padj_interaction": inter_padj[order],
                    "log2FC_SF_saline": lfc_s[order],
                    "log2FC_SF_BuOE": lfc_b[order],
                    "mean_log2_q3": tab["mean_log2_q3"].to_numpy()[order],
                }
            )
        )
    return pd.DataFrame(rows), pd.concat(top_rows, ignore_index=True)


# --------------------------------------------------------------------------- #
# Figures
# --------------------------------------------------------------------------- #
def figures(long: pd.DataFrame, summary: pd.DataFrame, per_region: dict[str, pd.DataFrame],
            att: pd.DataFrame, audit: pd.DataFrame) -> None:
    P.use_style()
    import matplotlib.pyplot as plt

    # --- DEG landscape -------------------------------------------------------- #
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 3.9),
                             gridspec_kw={"width_ratios": [1.35, 1]})
    ax = axes[0]
    labels = [lab for _, lab, _ in LANDSCAPE]
    width = 0.2
    x = np.arange(len(LANDSCAPE))
    for ri, region in enumerate(C.REGION_ORDER):
        vals = []
        for name, _, _ in LANDSCAPE:
            row = summary[(summary["region_short"] == region) & (summary["contrast"] == name)]
            vals.append(int(row["nominal_p0.01_n"].iloc[0]) if len(row) else 0)
        ax.bar(x + (ri - 1.5) * width, vals, width=width,
               color=P.REGION_COLORS[region], label=region, lw=0)
    exp_null = float(summary["expected_null_p0.01"].iloc[0])
    ax.axhline(exp_null, color="black", ls="--", lw=0.9)
    ax.text(len(LANDSCAPE) - 0.45, exp_null * 1.04,
            f"expected under the global null ({exp_null:.0f})",
            fontsize=6.3, ha="right", va="bottom")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=6.2, rotation=30, ha="right")
    ax.set_ylabel("Targets at nominal $P$ < 0.01")
    ax.legend(title="Region", fontsize=6.5, title_fontsize=6.5, ncols=2)
    ax.grid(axis="y", ls=":", lw=0.5, alpha=0.4)
    ax.set_title("Nominal signal against the null expectation", pad=4)

    ax = axes[1]
    piv = summary.pivot_table(index="contrast", columns="region_short",
                              values="fdr05_lfc1_n_deg", aggfunc="sum")
    piv = piv.reindex(index=[n for n, _, _ in LANDSCAPE], columns=C.REGION_ORDER)
    im = ax.imshow(piv.to_numpy(), cmap="Reds", vmin=0, vmax=max(1, np.nanmax(piv.to_numpy())))
    ax.set_xticks(range(len(C.REGION_ORDER)))
    ax.set_xticklabels(C.REGION_ORDER)
    ax.set_yticks(range(len(piv)))
    ax.set_yticklabels([lab for _, lab, _ in LANDSCAPE], fontsize=6.4)
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            v = piv.to_numpy()[i, j]
            ax.text(j, i, "0" if v == 0 else f"{int(v)}", ha="center", va="center",
                    fontsize=7.5, color="black")
    ax.set_title("FDR-significant targets\n(BH<0.05, |log$_2$FC|$\\geq$1)", pad=4)
    fig.tight_layout()
    P.save(fig, "DEG_landscape")

    # --- deposited vs correct FDR -------------------------------------------- #
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.5))
    ax = axes[0]
    idx = np.arange(len(audit))
    ax.bar(idx - 0.2, audit["n_sig_deposited_adj_lt_0.05"], width=0.4,
           color="#C44E52", label="deposited 'adjusted P' < 0.05", lw=0)
    ax.bar(idx + 0.2, audit["n_sig_recomputed_BH_lt_0.05"], width=0.4,
           color="#4C72B0", label="Benjamini-Hochberg < 0.05", lw=0)
    ax.set_xticks(idx)
    ax.set_xticklabels(
        [f"{r}\n{c.replace('SALvsBuOE','BuOE').replace('GCvsFLT','Flight')}"
         for r, c in zip(audit["region_short"], audit["contrast"])],
        fontsize=5.6,
    )
    ax.set_ylabel("Significant targets")
    ax.legend(fontsize=6.5)
    ax.grid(axis="y", ls=":", lw=0.5, alpha=0.4)
    ax.set_title("Deposited significance does not survive correct FDR control", pad=4)

    ax = axes[1]
    de = C.read_provided_de("CA", "GCvsFLT-SAL")
    ax.scatter(de["pvalue"], de["padj_deposited"], s=3, c="#7F8C8D", lw=0,
               alpha=0.4, rasterized=True)
    lim = [0, 1]
    ax.plot(lim, lim, color="black", lw=0.9, ls="--")
    ax.set_xlabel("Deposited raw $P$")
    ax.set_ylabel("Deposited 'adjusted $P$'")
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.02, 1)
    frac = audit.loc[audit["region_short"] == "CA1"].iloc[0]["frac_adj_below_raw_p"]
    ax.text(0.03, 0.97,
            f"CA1, Flight vs Ground (saline)\n{frac*100:.0f}% of targets fall below the\n"
            "identity line, which no correction\ncan produce",
            transform=ax.transAxes, fontsize=6.3, va="top")
    ax.set_title("Adjusted $P$ vs raw $P$ as deposited", pad=4)
    fig.tight_layout()
    P.save(fig, "Deposited_statistics_audit")

    # --- volcanoes: spaceflight in each stratum, per region ------------------- #
    for contrast, nice in [
        ("FLT_vs_GC_in_SAL", "Flight vs Ground, saline"),
        ("FLT_vs_GC_in_BuOE", "Flight vs Ground, BuOE"),
        ("interaction", "Flight x BuOE interaction"),
    ]:
        fig, axes = plt.subplots(1, 4, figsize=(12.4, 3.3))
        for ax, region in zip(axes, C.REGION_ORDER):
            g = long[(long["region_short"] == region) & (long["contrast"] == contrast)]
            P.volcano(
                ax,
                g["log2FC"].to_numpy(),
                g["pvalue"].to_numpy(),
                padj=g["padj_BH"].to_numpy(),
                genes=g["gene"].to_numpy(),
                lfc_cut=C.LFC_CUTOFF_RELAXED,
                title=region,
                label_n=6,
            )
        fig.suptitle(nice, fontsize=9.5, y=1.02)
        fig.tight_layout()
        P.save(fig, f"Volcano_{contrast}")

    # --- attenuation ---------------------------------------------------------- #
    fig, axes = plt.subplots(1, 3, figsize=(11.4, 3.5))
    ax = axes[0]
    region = "CA1"
    tab = per_region[region]
    ax.scatter(tab["log2FC_FLT_vs_GC_in_SAL"], tab["log2FC_FLT_vs_GC_in_BuOE"],
               s=3, c="#7F8C8D", lw=0, alpha=0.35, rasterized=True)
    lim = np.nanpercentile(
        np.abs(np.r_[tab["log2FC_FLT_vs_GC_in_SAL"], tab["log2FC_FLT_vs_GC_in_BuOE"]]), 99.5
    )
    ax.plot([-lim, lim], [-lim, lim], color="black", lw=0.9, ls="--", label="no attenuation")
    slope = float(
        att.loc[att["region_short"] == region, "rms_ratio_BuOE_over_saline"].iloc[0]
    )
    ax.plot([-lim, lim], [-lim * slope, lim * slope], color="#C44E52", lw=1.1,
            label=f"noise-corrected ratio {slope:.2f}")
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("log$_2$FC Flight vs Ground (saline)")
    ax.set_ylabel("log$_2$FC Flight vs Ground (BuOE)")
    ax.set_title(f"{region}: spaceflight response, saline vs BuOE", pad=4)
    ax.legend(fontsize=6.3, loc="upper left")

    ax = axes[1]
    x = np.arange(len(C.REGION_ORDER))
    r = att.set_index("region_short").reindex(C.REGION_ORDER)
    for off, key, color, lab in [
        (-0.16, "saline", "#C44E52", "saline"),
        (0.16, "BuOE", "#4C72B0", "BuOE"),
    ]:
        est = r[f"msq_true_{key}"].to_numpy()
        lo = est - r[f"msq_true_{key}_ci_lo"].to_numpy()
        hi = r[f"msq_true_{key}_ci_hi"].to_numpy() - est
        ax.errorbar(x + off, est, yerr=[lo, hi], fmt="o", color=color, ms=4.5,
                    capsize=2.5, lw=1.0, label=lab)
    ax.axhline(0, color="black", ls="--", lw=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels(C.REGION_ORDER)
    ax.set_ylabel("Noise-corrected mean squared\nspaceflight effect (log$_2$ units$^2$)")
    ax.legend(fontsize=6.8)
    ax.grid(axis="y", ls=":", lw=0.5, alpha=0.4)
    ax.set_title("Zero means indistinguishable from noise", pad=4)

    ax = axes[2]
    est = r["msq_diff_saline_minus_BuOE"].to_numpy()
    lo = est - r["msq_diff_ci_lo"].to_numpy()
    hi = r["msq_diff_ci_hi"].to_numpy() - est
    ax.errorbar(est, x, xerr=[lo, hi], fmt="o", color="#2C3E50", ms=5,
                capsize=2.5, lw=1.0)
    ax.axvline(0.0, color="black", ls="--", lw=0.9)
    ax.set_yticks(x)
    ax.set_yticklabels(C.REGION_ORDER)
    ax.set_xlabel("Mean squared effect, saline $-$ BuOE\n(>0 would indicate attenuation)")
    ax.set_title("Attenuation test, bootstrap 95% CI", pad=4)
    ax.grid(axis="x", ls=":", lw=0.5, alpha=0.4)
    fig.tight_layout()
    P.save(fig, "Attenuation_summary")

    # --- how much the verdict depends on the variance assumption -------------- #
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.3))
    ax = axes[0]
    ax.bar(x, r["median_var_ratio_BuOE_over_saline"], width=0.55, color="#8172B3", lw=0)
    ax.axhline(1.0, color="black", ls="--", lw=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels(C.REGION_ORDER)
    ax.set_ylabel("Median within-arm variance ratio,\nBuOE / saline")
    ax.set_ylim(0, max(1.6, float(r["median_var_ratio_BuOE_over_saline"].max()) * 1.15))
    ax.grid(axis="y", ls=":", lw=0.5, alpha=0.4)
    ax.set_title("The BuOE arm is the noisier arm", pad=4)

    ax = axes[1]
    ax.bar(x - 0.19, r["msq_diff_pooled_variance"], width=0.36,
           color="#BDC3C7", label="one pooled variance", lw=0)
    ax.bar(x + 0.19, r["msq_diff_saline_minus_BuOE"], width=0.36,
           color="#2C3E50", label="each arm's own variance", lw=0)
    ax.axhline(0, color="black", ls="--", lw=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels(C.REGION_ORDER)
    ax.set_ylabel("Mean squared effect,\nsaline $-$ BuOE")
    ax.legend(fontsize=6.5)
    ax.grid(axis="y", ls=":", lw=0.5, alpha=0.4)
    ax.set_title("Above zero = attenuation:\nthe sign depends on the assumption", pad=4)
    fig.tight_layout()
    P.save(fig, "Attenuation_variance_assumption")

    # --- selection-bias demonstration ---------------------------------------- #
    fig, ax = plt.subplots(figsize=(4.8, 3.3))
    x = np.arange(len(C.REGION_ORDER))
    ax.bar(x - 0.31, r["median_abs_lfc_top500_saline_SELECTION_BIASED"], width=0.19,
           color="#C44E52", label="saline, top 500 by saline $P$", lw=0)
    ax.bar(x - 0.11, r["median_abs_lfc_top500_BuOE_SELECTION_BIASED"], width=0.19,
           color="#E8A0A0", label="BuOE, same 500 targets", lw=0)
    ax.bar(x + 0.11, r["median_abs_lfc_random500_saline"], width=0.19,
           color="#4C72B0", label="saline, 500 random targets", lw=0)
    ax.bar(x + 0.31, r["median_abs_lfc_random500_BuOE"], width=0.19,
           color="#9DB8D8", label="BuOE, same random targets", lw=0)
    ax.set_xticks(x)
    ax.set_xticklabels(C.REGION_ORDER)
    ax.set_ylabel("Median |log$_2$FC| for Flight vs Ground")
    ax.legend(fontsize=5.9, loc="upper right")
    ax.grid(axis="y", ls=":", lw=0.5, alpha=0.4)
    ax.set_title("Selecting targets on the saline response\nmanufactures attenuation", pad=4)
    fig.tight_layout()
    P.save(fig, "Attenuation_selection_bias")


def main() -> int:
    C.ensure_dirs()
    samples = C.load_sample_table()

    print("Auditing the deposited statistics ...")
    audit = audit_deposited()
    C.write_csv(audit, "deposited_stats_audit.csv")
    pd.set_option("display.width", 220)
    pd.set_option("display.max_columns", 60)
    print(
        audit[
            ["region_short", "contrast", "frac_adj_below_raw_p",
             "n_sig_deposited_adj_lt_0.05", "n_sig_recomputed_BH_lt_0.05",
             "min_recomputed_BH"]
        ].to_string(index=False)
    )

    print("\nRe-deriving differential expression from the Q3 matrix ...")
    long, per_region = derive_de(samples)
    C.write_csv(long, "DE_all_contrasts_long.csv")

    cmp = compare_deposited(per_region)
    C.write_csv(cmp, "deposited_vs_derived.csv")
    print("\nReproduction of the deposited fold changes:")
    print(cmp.to_string(index=False))

    summary = summarise_degs(long)
    C.write_csv(summary, "DEG_summary_all.csv")
    print("\nDEG summary (FDR-controlled vs nominal, against the null expectation):")
    print(
        summary[
            ["region_short", "contrast", "fdr05_lfc1_n_deg", "min_padj_BH",
             "nominal_p0.01_n", "expected_null_p0.01", "obs_over_expected_p0.01"]
        ].to_string(index=False)
    )

    sig = long[long["padj_BH"] < C.FDR_CUTOFF]
    C.write_csv(sig.sort_values("padj_BH"), "significant_targets_FDR05.csv")

    sens = sensitivity(per_region)
    C.write_csv(sens, "sensitivity_minimum_detectable_effect.csv")
    print("\nSmallest effect this design could detect at 80% power:")
    print(
        sens[sens["contrast"] == "FLT_vs_GC_in_SAL"][
            ["region_short", "df_resid", "median_se_log2",
             "mde_log2FC_nominal_a0.05", "mde_log2FC_genomewide_a0.05_over_n",
             "mde_fold_change_genomewide_a0.05_over_n"]
        ].to_string(index=False)
    )

    print("\nTesting BuOE attenuation of the spaceflight response ...")
    att, top = attenuation(per_region, samples)
    C.write_csv(att, "attenuation_summary.csv")
    C.write_csv(top, "interaction_top_targets.csv")
    print(
        att[
            ["region_short", "msq_true_saline", "msq_true_saline_ci_lo",
             "msq_true_saline_ci_hi", "saline_effect_above_noise",
             "msq_diff_saline_minus_BuOE", "msq_diff_ci_lo", "msq_diff_ci_hi",
             "msq_diff_bootstrap_p", "n_interaction_fdr05", "verdict"]
        ].to_string(index=False)
    )
    print("\nWithin-arm variance check (pooled model assumes this ratio is 1):")
    print(
        att[
            ["region_short", "median_var_ratio_BuOE_over_saline",
             "msq_true_saline", "msq_true_saline_pooled_variance",
             "msq_true_BuOE", "msq_true_BuOE_pooled_variance"]
        ].to_string(index=False)
    )

    figures(long, summary, per_region, att, audit)
    return 0


if __name__ == "__main__":
    sys.exit(main())
