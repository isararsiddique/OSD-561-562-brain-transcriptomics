#!/usr/bin/env python3
"""Sensitivity analysis on the g:Profiler term tables.

g:Profiler's g:SCS correction returned nothing significant for any set, and
with all_results=True every reported p is the *corrected* value capped at 1.0,
so the table alone cannot say how close anything came. Each returned term
carries the full 2x2 contingency (term_size, query_size, intersection_size,
effective_domain_size), so the uncorrected one-sided hypergeometric p and the
odds ratio are recomputed here exactly, and a Benjamini-Hochberg correction --
deliberately less conservative than g:SCS -- is applied within each gene set.

This can only make significance easier to reach, so a null result here is a
stronger null than g:SCS alone.
"""
from __future__ import annotations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

HERE = Path(__file__).resolve().parents[1]
TAB = HERE / "tables"

MIN_TERM, MAX_TERM = 5, 500   # standard term-size window; huge terms are uninformative


def bh(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p values (step-up, monotone)."""
    p = np.asarray(p, dtype=float)
    m = p.size
    order = np.argsort(p, kind="mergesort")
    ranked = p[order] * m / np.arange(1, m + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty(m, dtype=float)
    q[order] = np.minimum(ranked, 1.0)
    return q


def main() -> None:
    a = pd.read_csv(TAB / "gprofiler_all_terms.csv.gz")
    N = a.effective_domain_size
    K = a.term_size
    n = a.query_size
    k = a.intersection_size
    # Uncorrected one-sided (over-representation) hypergeometric p.
    a["p_hypergeom_uncorrected"] = hypergeom.sf(k - 1, N, K, n)
    # Haldane-Anscombe corrected odds ratio from the 2x2 table.
    a_ = k + 0.5
    b_ = n - k + 0.5
    c_ = K - k + 0.5
    d_ = N - K - n + k + 0.5
    a["odds_ratio"] = (a_ * d_) / (b_ * c_)

    out = []
    for gs, sub in a.groupby("gene_set"):
        w = sub[(sub.term_size >= MIN_TERM) & (sub.term_size <= MAX_TERM)].copy()
        rej, q = None, bh(w.p_hypergeom_uncorrected.to_numpy())
        rej = q <= 0.05
        w["q_BH_within_set"] = q
        w["significant_BH_05"] = rej
        w["n_terms_in_BH_family"] = len(w)
        out.append(w)
        best = w.nsmallest(1, "p_hypergeom_uncorrected")
        print(
            f"{gs}: terms tested={len(w)}  min uncorrected p={w.p_hypergeom_uncorrected.min():.3g}"
            f"  min BH q={q.min():.3g}  n BH-significant={int(rej.sum())}"
        )
        if len(best):
            r = best.iloc[0]
            print(f"   strongest term: {r['source']} {r['native']} {r['name']}"
                  f" | {int(r.intersection_size)}/{int(r.term_size)} OR={r.odds_ratio:.2f}")

    res = pd.concat(out, ignore_index=True).sort_values(
        ["gene_set", "p_hypergeom_uncorrected"]
    )
    res.to_csv(TAB / "gprofiler_recomputed_with_BH.csv.gz", index=False)

    # Top 20 per set, for the supplementary table.
    top = pd.concat(
        [d.nsmallest(20, "p_hypergeom_uncorrected") for _, d in res.groupby("gene_set")],
        ignore_index=True,
    )
    cols = ["gene_set", "source", "native", "name", "term_size", "query_size",
            "intersection_size", "expected_overlap", "fold_enrichment",
            "odds_ratio", "p_hypergeom_uncorrected", "q_BH_within_set",
            "p_value", "significant"]
    top[[c for c in cols if c in top.columns]].rename(
        columns={"p_value": "p_gSCS_corrected", "significant": "significant_gSCS_05"}
    ).to_csv(TAB / "gprofiler_top20_per_set.csv", index=False)
    print("\nwrote gprofiler_recomputed_with_BH.csv.gz and gprofiler_top20_per_set.csv")


if __name__ == "__main__":
    main()
