#!/usr/bin/env python3
"""Resolve the orientation of the deposited log2 fold-change columns.

The NanoString DSP exports in GSE239336 are named ``GCvsFLT``/``SALvsBuOE`` but
the files themselves never state which level is the numerator, so the sign of
the ``Log2`` column is ambiguous.  Guessing it would silently invert every
biological conclusion, so it is determined here from the data: for each region
and contrast the deposited column is correlated against the group-mean
difference computed from the Q3-normalised matrix in both orientations.

Writes ``config/sign_convention.json`` (consumed by common.sign_for) and prints
a reconciliation report including the resolution of the deposited adjusted
p-values, which are rounded to two decimals in the export.

Usage:  python scripts/check_conventions.py
"""
from __future__ import annotations

import json
import sys

import numpy as np
import pandas as pd

import common as C


def main() -> int:
    samples = C.load_sample_table()
    report_rows = []
    votes: dict[str, list[int]] = {k: [] for k in C.PROVIDED_CONTRASTS}
    coverage_rows = []

    for study in C.STUDIES:
        log2, sub = C.load_expression(study.code, samples)
        group_of = dict(zip(sub["sample"], sub["group"]))
        by_group = {
            g: log2.loc[:, [s for s in log2.columns if group_of[s] == g]]
            for g in C.GROUPS
        }

        for contrast, spec in C.PROVIDED_CONTRASTS.items():
            de = C.read_provided_de(study.code, contrast)
            de = de.drop_duplicates("gene").set_index("gene")
            shared = [g for g in log2.index if g in de.index]

            num, den = spec["numerator"], spec["denominator"]
            forward = (
                by_group[num].loc[shared].mean(axis=1)
                - by_group[den].loc[shared].mean(axis=1)
            ).to_numpy()
            deposited = de.loc[shared, "log2FC_deposited"].to_numpy()

            ok = np.isfinite(forward) & np.isfinite(deposited)
            r_fwd = float(np.corrcoef(forward[ok], deposited[ok])[0, 1])
            r_rev = float(np.corrcoef(-forward[ok], deposited[ok])[0, 1])
            sign = 1 if r_fwd >= r_rev else -1
            votes[contrast].append(sign)

            slope = float(
                np.polyfit(deposited[ok] * sign, forward[ok], 1)[0]
            )
            report_rows.append(
                {
                    "region": study.short,
                    "osd": study.osd,
                    "contrast": contrast,
                    "dataset_name": de.attrs.get("dataset_name", ""),
                    "stated_orientation": f"{num} - {den}",
                    "r_as_deposited": round(r_fwd, 4),
                    "r_negated": round(r_rev, 4),
                    "resolved_sign": sign,
                    "slope_q3_vs_deposited": round(slope, 4),
                    "n_shared_genes": len(shared),
                }
            )

            if contrast == "GCvsFLT-SAL":
                padj = de.loc[shared, "padj_deposited"]
                pval = de.loc[shared, "pvalue"]
                coverage_rows.append(
                    {
                        "region": study.short,
                        "n_q3_targets": int(log2.shape[0]),
                        "n_de_targets": int(de.shape[0]),
                        "in_q3_not_in_de": ", ".join(
                            sorted(set(log2.index) - set(de.index))
                        ),
                        "n_unique_padj_deposited": int(padj.nunique()),
                        "n_unique_pvalue": int(pval.nunique()),
                        "min_padj_deposited": float(padj.min()),
                        "n_padj_exactly_zero": int((padj == 0).sum()),
                    }
                )

    report = pd.DataFrame(report_rows)
    coverage = pd.DataFrame(coverage_rows)

    pd.set_option("display.width", 200)
    pd.set_option("display.max_columns", 50)
    print("=" * 100)
    print("LOG2 ORIENTATION CHECK  (deposited column vs Q3-normalised group means)")
    print("=" * 100)
    print(report.to_string(index=False))

    print()
    print("=" * 100)
    print("TARGET COVERAGE AND P-VALUE RESOLUTION")
    print("=" * 100)
    print(coverage.to_string(index=False))

    signs: dict[str, int] = {}
    for contrast, vs in votes.items():
        if len(set(vs)) != 1:
            print(
                f"\nERROR: inconsistent sign across regions for {contrast}: {vs}",
                file=sys.stderr,
            )
            return 1
        signs[contrast] = vs[0]

    weak = report.loc[report[["r_as_deposited", "r_negated"]].abs().max(axis=1) < 0.9]
    if len(weak):
        print("\nWARNING: weak agreement for:")
        print(weak.to_string(index=False))

    C.CONFIG.mkdir(parents=True, exist_ok=True)
    out = {
        "description": (
            "Multiplier applied to the deposited NanoString 'Log2' column so that "
            "positive values mean higher in the numerator recorded in "
            "common.PROVIDED_CONTRASTS. Resolved empirically against the "
            "Q3-normalised matrix by scripts/check_conventions.py."
        ),
        "signs": signs,
        "evidence": report.to_dict(orient="records"),
        "coverage": coverage.to_dict(orient="records"),
    }
    path = C.CONFIG / "sign_convention.json"
    path.write_text(json.dumps(out, indent=2) + "\n")
    print(f"\nResolved signs: {signs}")
    print(f"Wrote {path.relative_to(C.ROOT)}")

    C.ensure_dirs()
    C.write_csv(report, "log2_orientation_check.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
