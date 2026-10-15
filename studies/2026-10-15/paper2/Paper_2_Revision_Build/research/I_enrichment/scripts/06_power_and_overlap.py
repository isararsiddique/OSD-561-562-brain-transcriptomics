#!/usr/bin/env python3
"""How large an effect could these tests have detected?

A null result is only informative if the test had power to see something. For
every cell class x gene set, this finds the smallest overlap k that would have
reached an uncorrected one-sided Fisher p <= 0.05 and the odds ratio that
implies -- the minimum detectable effect. The same is done for the functional
terms, as a minimum detectable fold-enrichment by term size.

Also reports the descriptive overlap between the flight and housing gene sets.
That overlap is NOT a clean test: FLT-GC and GC-VIV share the ground-control
arm, so genes with a noisy GC estimate enter both lists under the null. It is
given as a description with that dependence stated, not as an inference.
"""
from __future__ import annotations
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, hypergeom

HERE = Path(__file__).resolve().parents[1]
TAB = HERE / "tables"


def read_list(name: str) -> list[str]:
    return [g for g in (TAB / name).read_text().split("\n") if g.strip()]


def min_detectable(n_q: int, n_m: int, N: int, alpha: float = 0.05):
    """Smallest overlap reaching one-sided p<=alpha, and its Haldane OR."""
    for k in range(0, min(n_q, n_m) + 1):
        b, c, d = n_q - k, n_m - k, N - n_q - (n_m - k)
        _, p = fisher_exact([[k, b], [c, d]], alternative="greater")
        if p <= alpha:
            orr = ((k + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
            return k, p, orr
    return None, None, None


def main() -> None:
    background = read_list("background_measured_genes.txt")
    N = len(background)
    bg = set(background)
    queries = {
        "FLT_GC_BY": set(read_list("set_FLT_GC_BY.txt")) & bg,
        "FLT_GC_top100": set(read_list("set_FLT_GC_top100.txt")) & bg,
        "GC_VIV_top100": set(read_list("set_GC_VIV_top100.txt")) & bg,
    }

    ct = pd.read_csv(TAB / "celltype_enrichment_results.csv")
    rows = []
    for r in ct.itertuples():
        k, p, orr = min_detectable(r.n_query_in_background, r.n_markers_in_background, N)
        rows.append({
            "marker_definition": r.marker_definition,
            "gene_set": r.gene_set,
            "cell_class": r.cell_class,
            "n_query": r.n_query_in_background,
            "n_markers": r.n_markers_in_background,
            "observed_overlap": r.overlap,
            "expected_overlap": r.expected_overlap,
            "min_detectable_overlap_p05": k,
            "min_detectable_odds_ratio_p05": orr,
            "detectable_at_all": k is not None,
        })
    pw = pd.DataFrame(rows)
    pw.to_csv(TAB / "celltype_minimum_detectable_effect.csv", index=False)

    # Minimum detectable fold-enrichment for functional terms, by term size.
    gp = pd.read_csv(TAB / "gprofiler_recomputed_with_BH.csv.gz")
    term_rows = []
    for gs, sub in gp.groupby("gene_set"):
        nq = int(sub.query_size.iloc[0])
        Ne = int(sub.effective_domain_size.iloc[0])
        for ts in (5, 10, 25, 50, 100, 250, 500):
            k = int(hypergeom.isf(0.05, Ne, ts, nq)) + 1
            exp = ts * nq / Ne
            term_rows.append({
                "gene_set": gs, "query_size": nq, "effective_domain_size": Ne,
                "term_size": ts, "expected_overlap": exp,
                "min_overlap_for_uncorrected_p05": k,
                "min_detectable_fold_enrichment": k / exp if exp else np.nan,
            })
    pd.DataFrame(term_rows).to_csv(TAB / "term_minimum_detectable_effect.csv", index=False)

    # Descriptive flight/housing set overlap (shared GC arm - see docstring).
    f, g = queries["FLT_GC_top100"], queries["GC_VIV_top100"]
    by = queries["FLT_GC_BY"]
    ov = {}
    for lab, A, Bs in (("FLT_top_vs_GCVIV_top", f, g), ("FLT_BY_vs_GCVIV_top", by, g)):
        k = len(A & Bs)
        b, c, d = len(A) - k, len(Bs) - k, N - len(A) - (len(Bs) - k)
        orr, p = fisher_exact([[k, b], [c, d]], alternative="greater")
        ov[lab] = {
            "n_a": len(A), "n_b": len(Bs), "overlap": k,
            "expected_overlap_if_independent": len(A) * len(Bs) / N,
            "odds_ratio": orr, "p_one_sided_NOT_A_VALID_TEST": p,
            "caveat": "FLT-GC and GC-VIV share the ground-control arm, so the two "
                      "rankings are positively dependent under the null; this p "
                      "value is reported for transparency and must not be read as "
                      "evidence of a shared biological programme.",
            "shared_genes": sorted(A & Bs),
        }

    summary = {
        "background_n": N,
        "celltype_power": {
            "n_tests": int(len(pw)),
            "n_classes_where_no_overlap_could_ever_reach_p05": int((~pw.detectable_at_all).sum()),
            "median_min_detectable_odds_ratio": float(
                pw.min_detectable_odds_ratio_p05.median(skipna=True)
            ),
            "by_gene_set": {
                gs: {
                    "min_detectable_OR_range": [
                        float(d.min_detectable_odds_ratio_p05.min()),
                        float(d.min_detectable_odds_ratio_p05.max()),
                    ],
                    "median": float(d.min_detectable_odds_ratio_p05.median()),
                }
                for gs, d in pw[pw.marker_definition == "all_markers"].groupby("gene_set")
            },
        },
        "set_overlap_descriptive": ov,
    }
    (TAB / "power_and_overlap_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps({k: v for k, v in summary.items() if k != "set_overlap_descriptive"}, indent=2))
    for lab, v in ov.items():
        print(f"\n{lab}: overlap {v['overlap']} (expected if independent "
              f"{v['expected_overlap_if_independent']:.2f}), OR={v['odds_ratio']:.2f}")


if __name__ == "__main__":
    main()
