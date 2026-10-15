"""Independent calculations without importing the analysis implementation."""
from pathlib import Path
import argparse
from collections import defaultdict
import csv
import gzip
import hashlib
import itertools
import json
import math
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

B = Path(__file__).resolve().parents[1]
T = B / "tables"


def read_raw_counts(path):
    pairs, alternatives = {}, {"start": defaultdict(list), "end": defaultdict(list)}
    with gzip.open(path, "rt") as f:
        reader = csv.reader(f)
        samples = next(reader)[1:]
        for row in reader:
            chrom, strand, start, end = row[0].split(":")
            vec = np.array([int(x) for x in row[1:]], dtype=np.int64)
            assert row[0] not in pairs
            pairs[row[0]] = vec
            alternatives["start"][(chrom, strand, start)].append(vec)
            alternatives["end"][(chrom, strand, end)].append(vec)
    sums = {bound: {key: np.sum(vecs, axis=0) for key, vecs in alt.items()} for bound, alt in alternatives.items()}
    return samples, pairs, sums


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, required=True, help="Beta analysis/package root containing inputs/ and results/")
    args = parser.parse_args()
    RR = args.base.resolve()
    checks = []
    summary = json.loads((T / "summary.json").read_text())
    manifest = pd.read_csv(T / "input_manifest_sha256.csv")
    for row in manifest.itertuples():
        p = (RR if row.scope == "beta" else B) / row.path
        current = hashlib.sha256(p.read_bytes()).hexdigest()
        assert current == row.sha256_before == row.sha256_after
    checks.append(dict(check="All source SHA256 unchanged", passed=True, detail=f"{len(manifest)} inputs including frozen plan"))
    definitions = pd.read_csv(RR / "inputs/primary_all26_event_definitions.csv.gz", dtype={"chr": str}).set_index("event_id")
    for tag, path, stored_frac, stored_pair, stored_dom in [
        ("RRRM2", B / "inputs/OSD-561_animal_junction_counts.csv.gz",
         T / "RRRM2_all26_fractions.csv.gz", T / "RRRM2_all26_pair_counts.csv.gz", T / "RRRM2_all26_pair_dominance.csv.gz"),
        ("RR23", RR / "inputs/all_animal_junction_counts.csv.gz",
         RR / "results/all26_animal_fractions.csv.gz", None, None)]:
        samples, counts, totals = read_raw_counts(path)
        frac = pd.read_csv(stored_frac, index_col=0).loc[definitions.index, samples].to_numpy()
        pair = pd.read_csv(stored_pair, index_col=0).loc[definitions.index, samples].to_numpy() if stored_pair else None
        dom = pd.read_csv(stored_dom, index_col=0).loc[definitions.index, samples].to_numpy() if stored_dom else None
        zero = np.zeros(len(samples), dtype=int)
        for i, row in enumerate(definitions.itertuples()):
            a, z = counts.get(row.junction_a, zero), counts.get(row.junction_b, zero)
            p = a + z
            f = np.divide(a, p, out=np.full(len(samples), np.nan), where=p > 0)
            np.testing.assert_allclose(f, frac[i], atol=1e-14, rtol=0, equal_nan=True)
            if pair is not None:
                np.testing.assert_array_equal(p, pair[i])
                key = (row.chr, str(row.strand), str(row.boundary))
                tot = totals[row.boundary_type].get(key, zero)
                d = np.divide(p, tot, out=np.zeros(len(samples), dtype=float), where=tot > 0)
                np.testing.assert_allclose(d, dom[i], atol=1e-14, rtol=0)
        checks.append(dict(check=f"{tag} raw CSV reconstruction", passed=True,
                           detail=f"All {len(definitions)} rows × {len(samples)} animals checked with independent dictionary boundary sums"))
    m = pd.read_csv(T / "RRRM2_heldout_animal_scores.csv")
    s = m.learned_direction_score.to_numpy()
    event = pd.read_csv(T / "complete_common_panel_effects.csv.gz").set_index("event_id")
    x = pd.read_csv(T / "RRRM2_all26_fractions.csv.gz", index_col=0).loc[event.index, m["sample"]].to_numpy()
    projected = ((x - event.RR23_GC_center.to_numpy()[:, None]) * event.RR23_learned_direction.to_numpy()[:, None]).mean(0)
    np.testing.assert_allclose(projected, s, atol=1e-14, rtol=0)
    assert len(m) == 19 and m["Source Name"].is_unique
    r = spearmanr(event.RR23_delta_FLT_minus_GC, event.RRRM2_delta_stratified_FLT_minus_GC).statistic
    assert np.isclose(r, summary["spearman_rho"], atol=1e-14)
    checks.append(dict(check="Fixed full-panel projection and profile correlation", passed=True,
                       detail=f"All {len(event)} eligible events; every held-out animal score and rho reproduced"))
    # Independent direct combination sums; no permutation-weight matrices.
    parts, observed, event_observed = [], 0.0, np.zeros(len(event))
    for st in sorted(m.stratum.unique()):
        ids = np.flatnonzero(m.stratum.eq(st))
        nf = int(m.iloc[ids].exposure.eq("FLT").sum())
        ng = len(ids) - nf
        flt_ids = ids[m.iloc[ids].exposure.eq("FLT").to_numpy()]
        ground_ids = ids[m.iloc[ids].exposure.eq("GC").to_numpy()]
        w = len(ids) / len(m)
        observed += w * (s[flt_ids].mean() - s[ground_ids].mean())
        event_observed += w * (x[:, flt_ids].mean(1) - x[:, ground_ids].mean(1))
        values = []
        for selected in itertools.combinations(ids, nf):
            flt_sum = s[list(selected)].sum()
            ground_sum = s[ids].sum() - flt_sum
            values.append(w * (flt_sum / nf - ground_sum / ng))
        assert len(values) == math.comb(len(ids), nf)
        parts.append(values)
    null = np.array([sum(parts_here) for parts_here in itertools.product(*parts)])
    assert len(null) == 11200 == summary["RRRM2_exact_label_assignments"]
    np.testing.assert_allclose(event_observed, event.RRRM2_delta_stratified_FLT_minus_GC, atol=1e-14, rtol=0)
    np.testing.assert_allclose(np.sort(null), np.sort(pd.read_csv(T / "RRRM2_exact_score_null.csv.gz").score_delta), atol=1e-14, rtol=0)
    p = float(np.mean(null >= observed - 1e-12))
    p2 = float(np.mean(np.abs(null) >= abs(observed) - 1e-12))
    assert np.isclose(p, summary["learned_direction_score_exact_one_sided_p"])
    assert np.isclose(p2, summary["learned_direction_score_exact_two_sided_p"])
    assert np.isclose(observed, summary["learned_direction_score_delta"])
    checks.append(dict(check="Independent exact stratified score enumeration", passed=True,
                       detail=f"All 11,200 assignments reproduced directly; one-sided p={p}, two-sided p={p2}"))
    # Deterministic event checks spread across row order, independent of significance.
    check_rows = np.unique(np.linspace(0, len(event) - 1, 13, dtype=int))
    event_parts = []
    for st in sorted(m.stratum.unique()):
        ids = np.flatnonzero(m.stratum.eq(st))
        nf = int(m.iloc[ids].exposure.eq("FLT").sum())
        ng = len(ids) - nf
        w = len(ids) / len(m)
        vals = []
        for selected in itertools.combinations(ids, nf):
            a = x[check_rows][:, list(selected)].sum(1)
            b = x[check_rows][:, ids].sum(1) - a
            vals.append(w * (a / nf - b / ng))
        event_parts.append(vals)
    event_null = np.array([np.sum(v, axis=0) for v in itertools.product(*event_parts)])
    sign = event.RR23_learned_direction.to_numpy()[check_rows]
    obs = event_observed[check_rows]
    ep = ((event_null * sign) >= (obs * sign - 1e-12)).mean(0)
    np.testing.assert_allclose(ep, event.RRRM2_exact_directional_p.iloc[check_rows], atol=1e-14, rtol=0)
    # Independently derive Holm as max of ordered scaled p values.
    p_event = event.RRRM2_exact_directional_p.to_numpy()
    ordered = np.argsort(p_event)
    previous, h = 0.0, np.zeros(len(event))
    for position, j in enumerate(ordered):
        previous = max(previous, p_event[j] * (len(event) - position))
        h[j] = min(1.0, previous)
    np.testing.assert_allclose(h, event.RRRM2_directional_Holm_full_common_panel, atol=1e-14, rtol=0)
    checks.append(dict(check="Secondary event p values and Holm family", passed=True,
                       detail="13 deterministic event exact nulls and all 1,512 Holm adjusted values independently reproduced"))
    boots = pd.read_csv(T / "animal_bootstrap_draws.csv.gz")
    assert len(boots) == 2000
    for col, key in [("spearman_rho", "spearman_animal_bootstrap_CI95"),
                     ("same_direction_fraction", "same_direction_animal_bootstrap_CI95"),
                     ("fixed_score_delta", "learned_direction_score_animal_bootstrap_CI95")]:
        np.testing.assert_allclose(np.quantile(boots[col], [.025, .975]), summary[key], atol=1e-14, rtol=0)
    checks.append(dict(check="Saved bootstrap interval consistency", passed=True, detail="All three percentile intervals reproduced from 2,000 saved animal draws"))
    out = dict(all_checks_passed=True, checks=checks, independent_checker_imports_analysis_code=False,
               exact_null_one_sided_extreme_count=int((null >= observed - 1e-12).sum()),
               exact_null_two_sided_extreme_count=int((np.abs(null) >= abs(observed) - 1e-12).sum()))
    pd.DataFrame(checks).to_csv(T / "independent_checks.csv", index=False)
    (T / "independent_checks.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
