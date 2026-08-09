#!/usr/bin/env python3
"""Which targets recur across brain regions, and is that more than chance?

No target survives FDR control in any single region, so the only way to get a
usable candidate list from this experiment is to ask which targets move in the
same direction in several regions at once.

An important caveat is built into the output: the four regions are sampled from
the *same* tissue sections, so regions are repeated measures on the same three
animals per group, not independent replicates.  Recurrence across regions is
therefore weaker evidence than it would be across independent cohorts.  To keep
that honest, the observed recurrence is compared against a permutation null in
which the flight/ground labels are shuffled within each arm, preserving the
region correlation structure.

Outputs
  results/tables/csv/recurrent_targets.csv
  results/tables/csv/recurrence_null_test.csv
  results/figures/Recurrence.{pdf,png}

Usage:  python scripts/analyze_recurrence.py
"""
from __future__ import annotations

import itertools
import sys

import numpy as np
import pandas as pd

import common as C
import plotting as P

NEG_PROBE = "NegProbe-WTX"
P_CUT = 0.01
LFC_CUT = C.LFC_CUTOFF_RELAXED
N_PERM = 200

CONTRASTS = {
    "FLT_vs_GC_in_SAL": "Flight vs Ground (saline)",
    "FLT_vs_GC_in_BuOE": "Flight vs Ground (BuOE)",
    "interaction": "Flight x BuOE interaction",
}


def observed(long: pd.DataFrame, contrast: str) -> pd.DataFrame:
    sel = long[
        (long["contrast"] == contrast)
        & (long["pvalue"] < P_CUT)
        & (long["log2FC"].abs() >= LFC_CUT)
    ]
    if sel.empty:
        return pd.DataFrame(columns=["gene", "n_regions", "regions", "mean_log2FC",
                                     "consistent_direction", "contrast"])
    g = sel.groupby("gene").agg(
        n_regions=("region_short", "nunique"),
        regions=("region_short", lambda s: ",".join(sorted(set(s), key=C.REGION_ORDER.index))),
        mean_log2FC=("log2FC", "mean"),
        min_pvalue=("pvalue", "min"),
        consistent_direction=("log2FC", lambda s: bool((s > 0).all() or (s < 0).all())),
    ).reset_index()
    g["contrast"] = contrast
    return g.sort_values(["n_regions", "mean_log2FC"], ascending=[False, False])


def permutation_null(samples: pd.DataFrame, contrast: str, rng) -> np.ndarray:
    """Distribution of the number of targets recurring in >= 3 regions.

    Flight/Ground labels are permuted within each drug arm, independently per
    region but using the same relabelling across regions for a given section, so
    the between-region correlation induced by shared animals is preserved.
    """
    mats, metas = {}, {}
    for study in C.STUDIES:
        log2, sub = C.load_expression(study.code, samples)
        mats[study.code] = log2.drop(index=[NEG_PROBE], errors="ignore")
        metas[study.code] = sub.set_index("sample").loc[mats[study.code].columns]

    genes = mats[C.STUDIES[0].code].index
    counts = np.empty(N_PERM, dtype=int)

    for k in range(N_PERM):
        hits = pd.Series(0, index=genes, dtype=int)
        # One relabelling of the three sections per arm, reused across regions.
        perm = {arm: rng.permutation(3) for arm in C.TRT_LEVELS}
        for study in C.STUDIES:
            meta = metas[study.code].copy()
            new_sf = meta["spaceflight"].to_numpy().copy()
            for arm in C.TRT_LEVELS:
                idx = np.where(meta["treatment"].to_numpy() == arm)[0]
                labels = meta["spaceflight"].to_numpy()[idx]
                # Shuffle the six labels within the arm using the section permutation.
                order = np.concatenate([perm[arm], perm[arm] + 3])
                new_sf[idx] = labels[order % len(labels)]
            res = C.two_way_anova(
                mats[study.code].to_numpy(), new_sf, meta["treatment"].to_numpy()
            )
            key = contrast if contrast in ("interaction",) else contrast
            p = res[f"p_{key}"]
            lfc = res[f"log2FC_{key}"]
            hits += ((p < P_CUT) & (np.abs(lfc) >= LFC_CUT)).astype(int)
        counts[k] = int((hits >= 3).sum())
    return counts


def main() -> int:
    C.ensure_dirs()
    path = C.CSV_DIR / "DE_all_contrasts_long.csv"
    if not path.exists():
        print(f"ERROR: {path} missing -- run scripts/analyze_de.py first", file=sys.stderr)
        return 1
    long = pd.read_csv(path)
    samples = C.load_sample_table()

    tables, null_rows = [], []
    rng = np.random.default_rng(11)
    for contrast, label in CONTRASTS.items():
        obs = observed(long, contrast)
        obs.insert(1, "contrast_label", label)
        tables.append(obs)
        n_obs3 = int((obs["n_regions"] >= 3).sum()) if len(obs) else 0
        n_obs2 = int((obs["n_regions"] >= 2).sum()) if len(obs) else 0
        null = permutation_null(samples, contrast, rng)
        null_rows.append(
            {
                "contrast": contrast,
                "contrast_label": label,
                "criterion": f"nominal P<{P_CUT}, |log2FC|>={LFC_CUT}",
                "n_recurring_2plus_regions": n_obs2,
                "n_recurring_3plus_regions": n_obs3,
                "null_mean_3plus": float(null.mean()),
                "null_p95_3plus": float(np.percentile(null, 95)),
                "permutation_p": float((null >= n_obs3).mean()),
                "n_permutations": N_PERM,
            }
        )
        print(
            f"  {label:<28} recurring in >=3 regions: {n_obs3:>4}  "
            f"(null mean {null.mean():.1f}, perm p = {(null >= n_obs3).mean():.3f})"
        )

    rec = pd.concat(tables, ignore_index=True)
    nulls = pd.DataFrame(null_rows)
    C.write_csv(rec, "recurrent_targets.csv")
    C.write_csv(nulls, "recurrence_null_test.csv")

    P.use_style()
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.4))
    ax = axes[0]
    x = np.arange(len(nulls))
    ax.bar(x - 0.2, nulls["n_recurring_3plus_regions"], width=0.4, color="#2C3E50",
           label="observed", lw=0)
    ax.bar(x + 0.2, nulls["null_mean_3plus"], width=0.4, color="#BDC3C7",
           label="permutation null (mean)", lw=0)
    ax.set_xticks(x)
    ax.set_xticklabels(nulls["contrast_label"], fontsize=6.2, rotation=20, ha="right")
    ax.set_ylabel(f"Targets nominal in $\\geq$3 regions")
    ax.legend(fontsize=6.5)
    ax.grid(axis="y", ls=":", lw=0.5, alpha=0.4)
    for i, p in enumerate(nulls["permutation_p"]):
        ax.text(i, max(nulls["n_recurring_3plus_regions"][i], nulls["null_mean_3plus"][i]),
                f" perm $p$={p:.3f}", fontsize=6.0, ha="center", va="bottom")
    ax.set_title("Cross-region recurrence vs a label-permutation null", pad=4)

    ax = axes[1]
    top = rec[(rec["contrast"] == "FLT_vs_GC_in_SAL") & (rec["n_regions"] >= 3)]
    top = top.reindex(top["mean_log2FC"].abs().sort_values(ascending=False).index).head(18)
    if len(top):
        colors = ["#C0392B" if v > 0 else "#2471A3" for v in top["mean_log2FC"]]
        ax.barh(range(len(top)), top["mean_log2FC"], color=colors, height=0.7)
        ax.set_yticks(range(len(top)))
        ax.set_yticklabels(top["gene"], fontsize=6.2, style="italic")
        ax.invert_yaxis()
        ax.axvline(0, color="black", lw=0.7)
        ax.set_xlabel("Mean log$_2$FC across regions")
        ax.grid(axis="x", ls=":", lw=0.5, alpha=0.4)
    else:
        ax.text(0.5, 0.5, "no targets recur in $\\geq$3 regions", ha="center",
                va="center", transform=ax.transAxes, fontsize=7)
        ax.set_axis_off()
    ax.set_title("Flight vs Ground (saline):\nmost consistent candidates", pad=4)
    fig.tight_layout()
    P.save(fig, "Recurrence")

    print("\nTop recurring candidates, Flight vs Ground (saline):")
    show = rec[(rec["contrast"] == "FLT_vs_GC_in_SAL") & (rec["n_regions"] >= 3)]
    if len(show):
        print(
            show[["gene", "n_regions", "regions", "mean_log2FC",
                  "consistent_direction", "min_pvalue"]]
            .head(20).to_string(index=False)
        )
    else:
        print("  none")
    return 0


if __name__ == "__main__":
    sys.exit(main())
