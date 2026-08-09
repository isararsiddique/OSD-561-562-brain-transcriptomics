#!/usr/bin/env python3
"""Cohort description and quality control for the four GeoMx DSP regions.

Outputs
  results/tables/csv/sample_metadata.csv     one row per ROI (48)
  results/tables/csv/cohort_summary.csv      design balance per region
  results/tables/csv/qc_metrics.csv          sequencing + ROI + background QC
  results/tables/csv/design_confounds.csv    explicit confound audit
  results/figures/QC_*.{pdf,png}

Usage:  python scripts/analyze_qc.py
"""
from __future__ import annotations

import sys

import numpy as np
import pandas as pd

import common as C
import plotting as P

NEG_PROBE = "NegProbe-WTX"


def build_qc(samples: pd.DataFrame) -> pd.DataFrame:
    """Per-ROI QC, combining sequencing metrics with expression-based background."""
    rows = []
    for study in C.STUDIES:
        log2, sub = C.load_expression(study.code, samples)
        genes = log2.drop(index=[NEG_PROBE], errors="ignore")
        neg = log2.loc[NEG_PROBE] if NEG_PROBE in log2.index else None
        for s in log2.columns:
            col = genes[s]
            rec = {
                "sample": s,
                "region_short": study.short,
                "osd": study.osd,
                "median_log2_q3": float(col.median()),
                "iqr_log2_q3": float(col.quantile(0.75) - col.quantile(0.25)),
            }
            if neg is not None:
                nb = float(neg[s])
                rec["negprobe_log2_q3"] = nb
                rec["snr_median_over_negprobe_log2"] = float(col.median()) - nb
                rec["n_targets_above_negprobe"] = int((col > nb).sum())
                rec["frac_targets_above_negprobe"] = float((col > nb).mean())
                rec["n_targets_above_2x_negprobe"] = int((col > nb + 1.0).sum())
            rows.append(rec)
    qc = pd.DataFrame(rows)
    keep = [
        "sample", "title", "region_short", "region", "osd", "glds", "slide", "roi",
        "group", "spaceflight", "treatment", "rep",
        "raw_reads", "trimmed_reads", "stitched_reads", "aligned_reads",
        "dedup_reads", "align_rate", "dedup_rate", "seq_saturation",
        "area_um2", "nuclei", "roi_x", "roi_y", "isa_rrna_pct",
    ]
    keep = [k for k in keep if k in samples.columns]
    return samples[keep].merge(qc.drop(columns=["region_short", "osd"]), on="sample")


def cohort_summary(samples: pd.DataFrame, n_targets: int) -> pd.DataFrame:
    rows = []
    for study in C.STUDIES:
        sub = samples.loc[samples["region_code"] == study.code]
        rec = {
            "osd": study.osd,
            "glds": study.glds,
            "region": study.region,
            "region_short": study.short,
            "n_roi": len(sub),
            "n_slides": sub["slide"].nunique(),
            "n_targets": n_targets,
            "sex": "/".join(sorted(sub["sex"].dropna().unique())),
            "strain": "/".join(sorted(sub["strain"].dropna().unique())),
        }
        for g in C.GROUPS:
            rec[f"n_{g}"] = int((sub["group"] == g).sum())
        if "age_at_launch_wk" in sub:
            rec["age_at_launch_wk"] = ", ".join(
                str(int(v)) for v in sorted(sub["age_at_launch_wk"].dropna().unique())
            )
        if "mission_duration_d" in sub:
            rec["mission_duration_d"] = ", ".join(
                str(int(v)) for v in sorted(sub["mission_duration_d"].dropna().unique())
            )
        rows.append(rec)
    total = {
        "osd": "all four",
        "glds": "-",
        "region": "CA1 + DG + FCtx + Ctx",
        "region_short": "combined",
        "n_roi": len(samples),
        "n_slides": samples["slide"].nunique(),
        "n_targets": n_targets,
        "sex": "/".join(sorted(samples["sex"].dropna().unique())),
        "strain": "/".join(sorted(samples["strain"].dropna().unique())),
    }
    for g in C.GROUPS:
        total[f"n_{g}"] = int((samples["group"] == g).sum())
    rows.append(total)
    return pd.DataFrame(rows)


def confound_audit(samples: pd.DataFrame) -> pd.DataFrame:
    """Record the structural confounds that limit what this design can support."""
    rows = []
    xt = samples.groupby(["group", "slide"], observed=True).size().reset_index(name="n")
    per_group_slides = xt.groupby("group")["slide"].nunique()
    rows.append(
        {
            "issue": "group vs slide",
            "detail": "; ".join(
                f"{g}={s}" for g, s in zip(xt["group"], xt["slide"])
            ),
            "n_levels_per_group": int(per_group_slides.max()),
            "confounded": bool((per_group_slides == 1).all()),
            "consequence": (
                "Each experimental group sits on its own GeoMx slide, so slide / "
                "batch effects cannot be separated from treatment effects."
            ),
        }
    )
    per_region = samples.groupby("region_code", observed=True).apply(
        lambda d: d.groupby("group", observed=True)["rep"].nunique().min(),
        include_groups=False,
    )
    rows.append(
        {
            "issue": "biological replication",
            "detail": "; ".join(f"{k}={v}" for k, v in per_region.items()),
            "n_levels_per_group": int(per_region.min()),
            "confounded": False,
            "consequence": (
                "Three tissue sections (rep1-3) per group per region. Replication "
                "is at the section/ROI level; with one slide per group the "
                "effective independent unit is ambiguous."
            ),
        }
    )
    shared = (
        samples.groupby(["slide", "rep"], observed=True)["region_code"]
        .nunique()
        .value_counts()
        .to_dict()
    )
    rows.append(
        {
            "issue": "region is a within-section factor",
            "detail": f"regions per (slide, rep): {shared}",
            "n_levels_per_group": len(C.REGION_CODES),
            "confounded": False,
            "consequence": (
                "All four regions are sampled from the same sections, so "
                "cross-region comparisons are paired rather than independent."
            ),
        }
    )
    return pd.DataFrame(rows)


def figures(qc: pd.DataFrame) -> None:
    P.use_style()
    import matplotlib.pyplot as plt

    order = C.REGION_ORDER
    gcolors = [C.GROUP_COLORS[g] for g in C.GROUPS]

    # --- sequencing depth + saturation + background -------------------------- #
    fig, axes = plt.subplots(2, 3, figsize=(10.2, 5.6))
    metrics = [
        ("raw_reads", "Raw reads (millions)", 1e-6),
        ("align_rate", "Aligned / raw reads", 1.0),
        ("seq_saturation", "Sequencing saturation (%)", 1.0),
        ("dedup_reads", "Deduplicated reads (thousands)", 1e-3),
        ("nuclei", "Nuclei per ROI", 1.0),
        ("area_um2", "ROI area (10$^3$ $\\mu$m$^2$)", 1e-3),
    ]
    for ax, (col, label, scale) in zip(axes.ravel(), metrics):
        if col not in qc:
            ax.set_visible(False)
            continue
        for gi, g in enumerate(C.GROUPS):
            xs, ys = [], []
            for ri, region in enumerate(order):
                sel = qc[(qc["region_short"] == region) & (qc["group"] == g)]
                jitter = (gi - 1.5) * 0.17
                xs.extend(np.full(len(sel), ri + jitter))
                ys.extend(sel[col].to_numpy() * scale)
            ax.scatter(xs, ys, s=13, color=C.GROUP_COLORS[g], lw=0.3,
                       edgecolor="white", label=C.GROUP_LABELS[g] if col == "raw_reads" else None)
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels(order)
        ax.set_ylabel(label)
        ax.grid(axis="y", ls=":", lw=0.5, alpha=0.4)
    axes[0, 0].legend(loc="upper left", fontsize=6.2, ncols=2, columnspacing=0.8,
                      handletextpad=0.3)
    fig.tight_layout()
    P.save(fig, "QC_sequencing_and_roi_metrics")

    # --- expression background ---------------------------------------------- #
    if "negprobe_log2_q3" in qc:
        fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.1))
        for ax, (col, label) in zip(
            axes,
            [
                ("negprobe_log2_q3", "NegProbe-WTX (log$_2$ Q3)"),
                ("snr_median_over_negprobe_log2", "Median target $-$ NegProbe (log$_2$)"),
                ("frac_targets_above_negprobe", "Fraction of targets > NegProbe"),
            ],
        ):
            for gi, g in enumerate(C.GROUPS):
                xs, ys = [], []
                for ri, region in enumerate(order):
                    sel = qc[(qc["region_short"] == region) & (qc["group"] == g)]
                    xs.extend(np.full(len(sel), ri + (gi - 1.5) * 0.17))
                    ys.extend(sel[col].to_numpy())
                ax.scatter(xs, ys, s=14, color=C.GROUP_COLORS[g], lw=0.3, edgecolor="white",
                           label=C.GROUP_LABELS[g] if col == "negprobe_log2_q3" else None)
            ax.set_xticks(range(len(order)))
            ax.set_xticklabels(order)
            ax.set_ylabel(label)
            ax.grid(axis="y", ls=":", lw=0.5, alpha=0.4)
        axes[0].legend(loc="best", fontsize=6.2, ncols=2, columnspacing=0.8,
                       handletextpad=0.3)
        fig.tight_layout()
        P.save(fig, "QC_background_and_detection")

    # --- ROI map -------------------------------------------------------------- #
    fig, axes = plt.subplots(1, 4, figsize=(11.0, 3.2), sharey=True)
    for ax, slide in zip(axes, ["GC_SAL_2", "GC_TX_2", "FLT_SAL_2", "FLT_TX_2"]):
        sel = qc[qc["slide"] == slide]
        for region in order:
            s2 = sel[sel["region_short"] == region]
            ax.scatter(s2["roi_x"] / 1e3, s2["roi_y"] / 1e3, s=34,
                       color=P.REGION_COLORS[region], lw=0.4, edgecolor="white",
                       label=region if slide == "GC_SAL_2" else None)
        ax.set_title(slide.replace("_2", "").replace("TX", "BuOE"), pad=3)
        ax.set_xlabel("ROI x (10$^3$ px)")
        ax.invert_yaxis()
    axes[0].set_ylabel("ROI y (10$^3$ px)")
    axes[0].legend(loc="best", fontsize=6.5, title="Region", title_fontsize=6.5)
    fig.tight_layout()
    P.save(fig, "QC_roi_layout")


def main() -> int:
    C.ensure_dirs()
    samples = C.load_sample_table()
    print(f"Loaded {len(samples)} ROIs across {samples['region_code'].nunique()} regions")

    n_targets = C.read_q3(C.STUDIES[0].code).shape[0]
    qc = build_qc(samples)

    C.write_csv(samples, "sample_metadata.csv")
    C.write_csv(cohort_summary(samples, n_targets - 1), "cohort_summary.csv")
    C.write_csv(qc, "qc_metrics.csv")
    C.write_csv(confound_audit(samples), "design_confounds.csv")

    figures(qc)

    print("\nDesign balance (region x group):")
    print(
        pd.crosstab(samples["region_short"], samples["group"])
        .reindex(index=C.REGION_ORDER, columns=C.GROUPS)
        .to_string()
    )
    print("\nSlide x group (all groups on a single slide => confounded):")
    print(pd.crosstab(samples["slide"], samples["group"]).to_string())
    print(
        f"\nSequencing depth: {qc['raw_reads'].min()/1e6:.1f}-"
        f"{qc['raw_reads'].max()/1e6:.1f} M raw reads; "
        f"saturation {qc['seq_saturation'].min():.1f}-{qc['seq_saturation'].max():.1f}%"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
