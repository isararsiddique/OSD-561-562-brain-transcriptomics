#!/usr/bin/env python3
"""Design calculator for rodent spaceflight junction studies: how many animals
per group are needed before a LOCAL (per-cluster) discovery is attainable at
all, and how many for 80% power.

Written for Paper 2.  Self-contained and seeded.  Three parts.

PART 1 -- ATTAINABILITY (exact arithmetic, no simulation)
  A balanced n vs n two-sided exact allocation test cannot return a p value
  below p_min = 2 / C(2n, n).  Equivalently, with mirror-pair reduction the
  test enumerates R = C(2n-1, n-1) = C(2n, n)/2 distinct allocation
  statistics and p = count / R, so p_min = 1 / R.

  A block of k tests tied at that floor attains Benjamini-Yekutieli adjusted
  value q(k) = p_min * m * H_m / k, which DECREASES in k.  So no BY rejection
  is possible at level alpha unless the family holds at least

        k*(m, n) = ceil( m * H_m / (R * alpha) )

  tests at the floor.  k* depends only on the design and the family size --
  never on the effect size, the sequencing depth or the biology.  Two derived
  design thresholds are reported:

    n_BY(m)   smallest n with k* = 1, i.e. the smallest design at which ONE
              isolated cluster can ever be reported under BY;
    n_Holm(m) smallest n with p_min * m <= alpha, the same question under Holm.

  The same arithmetic constrains the IMPLEMENTATION, not only the design: a
  random-permutation test with R drawn allocations has floor 1/R too, so the
  permutation budget must satisfy

        R >= R_min(m) = ceil( m * H_m / alpha )

  before a single test can clear BY, however many animals were flown.

PART 2 -- POWER (Monte Carlo)
  Dirichlet-multinomial compositional counts for one contrast of 923
  three-alternative clusters, S = 80 of them carrying a designated group
  separation, corrected inside the manuscript's real family m = 2,769
  (923 clusters x 3 contrasts).  Allocation p values are EXACT by complete
  mirror-reduced enumeration where C(2n, n) <= 200,000 (n <= 10), and
  estimated from R_MC = 480,000 random balanced allocations INCLUDING THE
  OBSERVED ONE for n = 12 and n = 15, where enumeration is infeasible.
  R_MC = 480,000 is chosen to exceed R_min(2769) = 470,932, so the resampling
  floor is NOT the binding constraint at those designs and the power curve
  stays design-limited throughout.

  The 80 designated clusters are permuted from counts.  The remaining
  m - 80 = 2,689 tests of the family are drawn from their exact discrete
  permutation null, which is uniform on {1/R, 2/R, ..., 1} under label
  exchangeability.

PART 3 -- FIGURE AND RECOMMENDATION TABLE
  Figure_11_design.pdf/.png, 183 mm wide, and a table of the n per group
  needed for attainability and for 80% per-effect power at the calibrated
  values.

Usage
  python3 design_calculator.py --reps 60                  # full run
  python3 design_calculator.py --reps 10 --tag check      # reproducibility
  python3 design_calculator.py --bench                    # timing probe
  python3 design_calculator.py --figure --tag main        # figure only
"""

from __future__ import annotations

# Thread pinning must happen before numpy imports its BLAS.
import os

for _v in (
    "OPENBLAS_NUM_THREADS",
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ[_v] = "1"

import argparse
import csv
import hashlib
import json
import math
import platform
import resource
import sys
import time
from itertools import combinations
from pathlib import Path

import numpy as np

# --------------------------------------------------------------------------
# Fixed constants.  Calibration values are quoted from the cohort, not tuned.
# --------------------------------------------------------------------------
ALPHA = 0.05
M_FAMILY = 2769           # 923 clusters x 3 contrasts -- the real family
N_CLUSTERS = 923
N_ALT = 3
S_DESIGNATED = 80         # designated true effects among the 923 clusters

# Real calibration (see D_summary.md and the cohort tables):
TV_CAL = 17.4             # candidate effect total variation, median, pp
THETA_CAL = 0.017         # DM overdispersion, median (IQR 0.01024-0.02753)
DEPTH_CAL = 140           # median selected-alternative reads per animal-cluster
DEPTH_LOGSD = 0.70
DEPTH_MIN = 20            # manuscript coverage rule

# Grids
N_GRID = (6, 8, 9, 10, 12, 15)
M_GRID = (500, 923, M_FAMILY, 5000)
N_ATTAIN = tuple(range(4, 16))
THETA_GRID = (0.010, 0.017, 0.028)
DEPTH_GRID = (70, 140, 280)
TV_GRID = (10.0, 17.4, 25.0)

# Permutation machinery
EXACT_MAX = 200_000       # enumerate exactly while C(2n, n) <= this
R_MC = 480_000            # random mirror-pair allocations otherwise
CHUNK = 8192
TIE_TOL = 1e-12

SEED = 20261015
SEED_SIGNS = 777_20261015

OKABE = {
    "orange": "#E69F00",
    "skyblue": "#56B4E9",
    "green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "purple": "#CC79A7",
    "grey": "#767676",
    "black": "#000000",
}


# --------------------------------------------------------------------------
# PART 1 -- exact attainability arithmetic
# --------------------------------------------------------------------------
def harmonic(m: int) -> float:
    """H_m = sum_{i=1..m} 1/i, summed ascending in float64."""
    return float(np.sum(1.0 / np.arange(1, m + 1, dtype=np.float64)))


def n_rows_exact(n: int) -> int:
    """Distinct mirror-reduced allocation statistics for a balanced n vs n."""
    return math.comb(2 * n - 1, n - 1)


def p_floor(n: int) -> float:
    """Smallest attainable two-sided exact allocation p value, 2/C(2n, n)."""
    return 2.0 / math.comb(2 * n, n)


def k_star_from_rows(m: int, rows: int, hm: float | None = None,
                     alpha: float = ALPHA) -> dict:
    """Floor ties needed for any BY rejection, given `rows` attainable levels."""
    hm = harmonic(m) if hm is None else hm
    pmin = 1.0 / rows
    raw = pmin * m * hm / alpha
    k = max(1, math.ceil(raw))
    return {
        "m": m,
        "n_levels": rows,
        "p_min": pmin,
        "harmonic_H_m": hm,
        "k_star_real": raw,
        "k_star": k,
        "attainable_within_family": bool(k <= m),
        "q_at_k_star": pmin * m * hm / k,
        "q_at_k_star_minus_1": (pmin * m * hm / (k - 1)) if k > 1 else float("nan"),
    }


def attainability_table() -> list:
    """Exact grid over n = 4..15 and the four family sizes."""
    hms = {m: harmonic(m) for m in M_GRID}
    rows = []
    for m in M_GRID:
        for n in N_ATTAIN:
            r = n_rows_exact(n)
            rec = k_star_from_rows(m, r, hms[m])
            holm_single = rec["p_min"] * m <= ALPHA
            rows.append({
                "m": m,
                "n_per_group": n,
                "n_total_animals": 2 * n,
                "n_allocations_C2nn": math.comb(2 * n, n),
                "n_mirror_levels": r,
                "p_min_two_sided": rec["p_min"],
                "harmonic_H_m": rec["harmonic_H_m"],
                "k_star_real": rec["k_star_real"],
                "k_star": rec["k_star"],
                "k_star_exceeds_family": not rec["attainable_within_family"],
                "single_test_passes_BY": rec["k_star"] == 1,
                "bonferroni_holm_single_p": rec["p_min"] * m,
                "single_test_passes_Holm": bool(holm_single),
            })
    return rows


def design_thresholds(table: list) -> dict:
    """Smallest n with k* = 1 (BY) and with a single Holm rejection, per m."""
    out = {}
    for m in M_GRID:
        sub = sorted([r for r in table if r["m"] == m], key=lambda r: r["n_per_group"])
        by = next((r["n_per_group"] for r in sub if r["single_test_passes_BY"]), None)
        holm = next((r["n_per_group"] for r in sub if r["single_test_passes_Holm"]), None)
        hm = harmonic(m)
        out[str(m)] = {
            "m": m,
            "harmonic_H_m": hm,
            "n_min_single_BY": by,
            "n_min_single_Holm": holm,
            "animals_total_single_BY": None if by is None else 2 * by,
            "animals_total_single_Holm": None if holm is None else 2 * holm,
            "k_star_at_9v9": next(r["k_star"] for r in sub if r["n_per_group"] == 9),
            "permutation_budget_R_min": math.ceil(m * hm / ALPHA),
        }
    return out


def by_discoveries(p_sorted: np.ndarray, m: int, hm: float,
                   alpha: float = ALPHA) -> int:
    """Benjamini-Yekutieli step-up count on an ascending p vector of length m."""
    thresh = np.arange(1, m + 1, dtype=np.float64) * alpha / (m * hm)
    ok = np.nonzero(p_sorted <= thresh)[0]
    return int(ok[-1] + 1) if ok.size else 0


# --------------------------------------------------------------------------
# PART 2 -- allocation sets
# --------------------------------------------------------------------------
_SIGN_CACHE: dict[int, tuple[np.ndarray, str]] = {}


def _exact_signs(n: int) -> np.ndarray:
    """All mirror-reduced balanced +-1 vectors; row 0 is the observed labelling."""
    total = 2 * n
    combos = list(combinations(range(1, total), n - 1))
    s = np.full((len(combos), total), -1, dtype=np.int8)
    s[:, 0] = 1
    for row, cols in enumerate(combos):
        s[row, list(cols)] = 1
    assert s.shape[0] == n_rows_exact(n)
    assert np.all(s.sum(axis=1) == 0)
    # Put the observed labelling (+1 for animals 0..n-1) in row 0.
    obs = np.r_[np.ones(n, dtype=np.int8), -np.ones(n, dtype=np.int8)]
    hit = np.nonzero((s == obs).all(axis=1))[0]
    assert hit.size == 1
    if hit[0] != 0:
        s[[0, hit[0]]] = s[[hit[0], 0]]
    assert np.array_equal(s[0], obs)
    return s


def _mc_signs(n: int, rows: int) -> np.ndarray:
    """`rows` random mirror-pair representatives, row 0 = the observed labelling.

    Every row fixes animal 0 to +1 (one representative per mirror pair) and
    places the remaining n-1 plus signs uniformly among the other 2n-1 animals.
    Row 0 is the observed labelling, so the returned p = count / rows is the
    standard Monte Carlo permutation p value with the observed allocation
    included; its floor is 1 / rows.
    """
    total = 2 * n
    rng = np.random.default_rng(np.random.SeedSequence([SEED_SIGNS, n, rows]))
    s = np.full((rows, total), -1, dtype=np.int8)
    s[:, 0] = 1
    s[0, 1:n] = 1                                  # observed labelling
    step = 20000
    for a in range(1, rows, step):
        b = min(a + step, rows)
        u = rng.random((b - a, total - 1))
        idx = np.argpartition(u, n - 2, axis=1)[:, : n - 1] + 1
        np.put_along_axis(s[a:b], idx, np.int8(1), axis=1)
    assert np.all(s.sum(axis=1) == 0)
    assert np.array_equal(s[0], np.r_[np.ones(n, dtype=np.int8),
                                      -np.ones(n, dtype=np.int8)])
    return s


def signs_for(n: int) -> tuple[np.ndarray, str]:
    """Cached allocation set for design size n, plus the mode used."""
    if n not in _SIGN_CACHE:
        if math.comb(2 * n, n) <= EXACT_MAX:
            _SIGN_CACHE[n] = (_exact_signs(n), "exact")
        else:
            _SIGN_CACHE[n] = (_mc_signs(n, R_MC), "monte_carlo")
    return _SIGN_CACHE[n]


def allocation_p(signs: np.ndarray, z: np.ndarray, n_feat: int,
                 chunk: int = CHUNK) -> np.ndarray:
    """Two-sided allocation p values, p = count / rows, one per feature.

    z has shape (2n, 3 * n_feat): the within-animal square-root composition
    vectors used by the manuscript's equal-cluster kernel.  The per-cluster
    statistic is T_f(s) = ||sum_i s_i z_fi||^2, whose ranking over allocations
    is identical to the squared distance between transformed group means, and
    which is invariant to a global sign flip (hence the mirror reduction).
    """
    rows = signs.shape[0]
    counts = np.zeros(n_feat, dtype=np.int64)
    t_obs = None
    for a in range(0, rows, chunk):
        blk = signs[a : a + chunk].astype(np.float64)
        v = blk @ z
        t = np.einsum("rfj,rfj->rf", v.reshape(blk.shape[0], n_feat, N_ALT),
                      v.reshape(blk.shape[0], n_feat, N_ALT))
        if t_obs is None:                     # row 0 is the observed labelling
            t_obs = t[0].copy()
        counts += np.count_nonzero(t >= t_obs * (1.0 - TIE_TOL), axis=0)
        del blk, v, t
    return counts / float(rows)


# --------------------------------------------------------------------------
# Data generation
# --------------------------------------------------------------------------
def base_compositions(rng: np.random.Generator, n_feat: int) -> np.ndarray:
    pi = rng.dirichlet(np.array([5.0, 3.0, 2.0]), size=n_feat)
    pi = np.clip(pi, 0.06, 0.88)
    return pi / pi.sum(axis=1, keepdims=True)


def split_group_means(rng: np.random.Generator, pi: np.ndarray, tv_pp: float):
    """Sum-zero separation of the two group means with total variation tv_pp."""
    d = rng.normal(size=pi.shape)
    d -= d.mean(axis=1, keepdims=True)
    l1 = np.abs(d).sum(axis=1, keepdims=True)
    d = d / l1 * (2.0 * tv_pp / 100.0)        # TV = 0.5 * sum|delta|
    a = pi + d / 2.0
    b = pi - d / 2.0
    bad = (a.min(axis=1) < 0.02) | (b.min(axis=1) < 0.02)
    a[bad] = pi[bad]
    b[bad] = pi[bad]
    return a, b, bad


def draw_counts(rng: np.random.Generator, mean_a: np.ndarray, mean_b: np.ndarray,
                n: int, theta: float, depth_median: float) -> np.ndarray:
    """Dirichlet-multinomial counts, shape (2n, n_feat, 3)."""
    n_feat = mean_a.shape[0]
    conc = (1.0 - theta) / theta
    out = np.empty((2 * n, n_feat, N_ALT), dtype=np.int64)
    depth = np.maximum(
        DEPTH_MIN,
        np.round(depth_median * np.exp(rng.normal(0.0, DEPTH_LOGSD,
                                                  size=(2 * n, n_feat)))),
    ).astype(np.int64)
    for i in range(2 * n):
        means = mean_a if i < n else mean_b
        p = rng.standard_gamma(means * conc)
        p /= p.sum(axis=1, keepdims=True)
        out[i] = rng.multinomial(depth[i], p)
    return out


def fractions_to_z(counts: np.ndarray) -> np.ndarray:
    tot = counts.sum(axis=2, keepdims=True).astype(np.float64)
    tot[tot == 0] = 1.0
    return np.sqrt(counts / tot)


# --------------------------------------------------------------------------
# One replicate of one design cell
# --------------------------------------------------------------------------
def run_replicate(seed_seq, n: int, tv_pp: float, theta: float,
                  depth_median: float, hm: float) -> dict:
    rng = np.random.default_rng(seed_seq)
    signs, mode = signs_for(n)
    rows = signs.shape[0]

    pi = base_compositions(rng, S_DESIGNATED)
    mean_a, mean_b, dropped = split_group_means(rng, pi, tv_pp)
    is_true = ~dropped
    s_real = int(is_true.sum())

    counts = draw_counts(rng, mean_a, mean_b, n, theta, depth_median)
    z = fractions_to_z(counts).reshape(2 * n, S_DESIGNATED * N_ALT)
    p_loc = allocation_p(signs, z, S_DESIGNATED)

    # Remaining family members: exact discrete permutation null.
    n_extra = M_FAMILY - S_DESIGNATED
    p_extra = rng.integers(1, rows + 1, size=n_extra) / float(rows)

    p_family = np.sort(np.concatenate([p_loc, p_extra]))
    k = by_discoveries(p_family, M_FAMILY, hm)
    pmin = 1.0 / rows
    floor_tol = pmin * (1.0 + TIE_TOL)

    if k > 0:
        cut = p_family[k - 1]
        n_true_rej = int(np.count_nonzero(p_loc[is_true] <= cut))
        n_false_rej = k - n_true_rej
    else:
        n_true_rej = 0
        n_false_rej = 0

    return {
        "n_per_group": n,
        "tv_pp": tv_pp,
        "theta": theta,
        "depth_median": depth_median,
        "perm_mode": mode,
        "n_levels": rows,
        "s_designated": S_DESIGNATED,
        "s_realized": s_real,
        "discoveries": k,
        "discoveries_true": n_true_rej,
        "discoveries_false": n_false_rej,
        "fdp": (n_false_rej / k) if k > 0 else 0.0,
        "per_effect_rejected": n_true_rej,
        "floor_ties_total": int(np.count_nonzero(p_loc <= floor_tol)
                                + np.count_nonzero(p_extra <= floor_tol)),
        "floor_ties_true": int(np.count_nonzero(p_loc[is_true] <= floor_tol)),
        "min_p_family": float(p_family[0]),
        "median_depth": float(np.median(counts.sum(axis=2))),
        "mean_tv_realized_pp": float(
            np.mean(0.5 * np.abs(mean_a[is_true] - mean_b[is_true]).sum(1)) * 100.0
        ) if s_real else float("nan"),
    }


# --------------------------------------------------------------------------
# Cell list -- fixed and deterministic, so reduced re-runs reproduce row for row
# --------------------------------------------------------------------------
SCENARIOS = (
    ("favourable", 0.010, 280),
    ("calibrated", THETA_CAL, DEPTH_CAL),
    ("adverse", 0.028, 70),
)


def build_cells() -> list:
    """Reduced grid.  The reduction is stated in H_summary.md.

    sweep_band  n x {favourable, calibrated, adverse} at TV = 17.4 pp.  The
                three scenarios are the extremes and the centre of the
                theta x depth grid, so they bracket all nine combinations.
    sweep_tv    n x TV in {10, 25} pp at the calibrated theta and depth.
    sweep_grid  the full theta x depth 3 x 3 at n = 9 and n = 10 (cheap, exact
                enumeration), which verifies that the nine combinations do fall
                inside the band drawn from the three scenarios.
    """
    cells = []
    for name, theta, depth in SCENARIOS:
        for n in N_GRID:
            cells.append({"sweep": "band", "scenario": name, "n_per_group": n,
                          "tv_pp": TV_CAL, "theta": theta, "depth_median": depth})
    for tv in (10.0, 25.0):
        for n in N_GRID:
            cells.append({"sweep": "tv", "scenario": f"tv{tv:g}", "n_per_group": n,
                          "tv_pp": tv, "theta": THETA_CAL,
                          "depth_median": DEPTH_CAL})
    for theta in THETA_GRID:
        for depth in DEPTH_GRID:
            for n in (9, 10):
                cells.append({"sweep": "grid",
                              "scenario": f"th{theta:g}_d{depth:g}",
                              "n_per_group": n, "tv_pp": TV_CAL,
                              "theta": theta, "depth_median": depth})
    return cells


def aggregate(cell: dict, recs: list) -> dict:
    d = np.array([r["discoveries"] for r in recs])
    dt = np.array([r["discoveries_true"] for r in recs])
    df = np.array([r["discoveries_false"] for r in recs])
    sr = np.array([r["s_realized"] for r in recs])
    fdp = np.array([r["fdp"] for r in recs])
    ties = np.array([r["floor_ties_total"] for r in recs])
    rows = recs[0]["n_levels"]
    ks = k_star_from_rows(M_FAMILY, rows)
    return {
        **{k: cell[k] for k in ("sweep", "scenario", "n_per_group", "tv_pp",
                                "theta", "depth_median")},
        "replicates": len(recs),
        "perm_mode": recs[0]["perm_mode"],
        "n_levels": rows,
        "k_star": ks["k_star"],
        "p_min": ks["p_min"],
        "s_realized_mean": float(sr.mean()),
        "expected_discoveries": float(d.mean()),
        "discoveries_sd": float(d.std(ddof=1)) if len(d) > 1 else 0.0,
        "discoveries_median": float(np.median(d)),
        "discoveries_max": int(d.max()),
        "prob_at_least_one": float(np.mean(d > 0)),
        "prob_at_least_one_true": float(np.mean(dt > 0)),
        "per_effect_power": float(dt.sum() / sr.sum()),
        "fdr_mean_fdp": float(fdp.mean()),
        "false_fraction_of_calls": float(df.sum() / d.sum()) if d.sum() else 0.0,
        "floor_ties_mean": float(ties.mean()),
        "median_depth_mean": float(np.mean([r["median_depth"] for r in recs])),
        "tv_realized_pp_mean": float(np.nanmean([r["mean_tv_realized_pp"]
                                                 for r in recs])),
    }


# --------------------------------------------------------------------------
# Recommendation table
# --------------------------------------------------------------------------
def recommend(agg: list, thresholds: dict) -> list:
    """n per group needed for attainability and for 80% per-effect power."""
    out = []
    for name, theta, depth in SCENARIOS:
        sub = sorted([a for a in agg if a["sweep"] == "band"
                      and a["scenario"] == name], key=lambda a: a["n_per_group"])
        out.append(_rec_row(f"TV 17.4 pp, {name} (theta {theta:g}, depth {depth:g})",
                            sub, thresholds))
    for tv in (10.0, 25.0):
        sub = sorted([a for a in agg if a["sweep"] == "tv"
                      and a["scenario"] == f"tv{tv:g}"],
                     key=lambda a: a["n_per_group"])
        out.append(_rec_row(f"TV {tv:g} pp, calibrated theta/depth", sub, thresholds))
    return out


def _interp_n(sub: list, key: str, target: float):
    """Smallest n on the grid reaching `target`, plus a guarded interpolation.

    Returns (grid_n, interp_n, interp_is_meaningful).  The interpolation is
    reported as meaningful ONLY when the bracketing designs share the same k*.
    Across a k* change the power curve has a genuine discontinuity -- at n = 8
    nothing is attainable at m = 2,769 whatever the effect size -- so a linear
    crossing between n = 8 and n = 9 is an artefact of the straight line, not
    an estimate of the design that reaches the target.
    """
    ns = [a["n_per_group"] for a in sub]
    ys = [a[key] for a in sub]
    ks = [a["k_star"] for a in sub]
    hit = next((n for n, y in zip(ns, ys) if y >= target), None)
    interp, ok = None, False
    for i in range(1, len(ns)):
        if ys[i - 1] < target <= ys[i]:
            f = (target - ys[i - 1]) / (ys[i] - ys[i - 1])
            interp = round(ns[i - 1] + f * (ns[i] - ns[i - 1]), 1)
            ok = ks[i - 1] == ks[i]
            break
    return hit, interp, bool(ok)


def _rec_row(label: str, sub: list, thresholds: dict) -> dict:
    t = thresholds[str(M_FAMILY)]
    n80, n80i, ok80 = _interp_n(sub, "per_effect_power", 0.80)
    n50, n50i, ok50 = _interp_n(sub, "per_effect_power", 0.50)
    n95, n95i, ok95 = _interp_n(sub, "prob_at_least_one", 0.95)
    return {
        "scenario": label,
        "n_attainable_single_BY": t["n_min_single_BY"],
        "n_attainable_single_Holm": t["n_min_single_Holm"],
        "k_star_at_9v9": t["k_star_at_9v9"],
        "n_grid_prob_any_ge_95": n95,
        "n_interp_prob_any_ge_95": n95i,
        "interp_any95_crosses_kstar_step": (not ok95) if n95i is not None else None,
        "n_grid_per_effect_power_ge_50": n50,
        "n_interp_per_effect_power_ge_50": n50i,
        "interp_power50_crosses_kstar_step": (not ok50) if n50i is not None else None,
        "n_grid_per_effect_power_ge_80": n80,
        "n_interp_per_effect_power_ge_80": n80i,
        "interp_power80_crosses_kstar_step": (not ok80) if n80i is not None else None,
        "per_effect_power_at_9v9": next(
            (a["per_effect_power"] for a in sub if a["n_per_group"] == 9), None),
        "expected_discoveries_at_9v9": next(
            (a["expected_discoveries"] for a in sub if a["n_per_group"] == 9), None),
        "permutation_budget_R_min": t["permutation_budget_R_min"],
    }


# --------------------------------------------------------------------------
# PART 3 -- figure
# --------------------------------------------------------------------------
def make_figure(out: Path, tag: str = "main") -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7,
        "axes.labelsize": 7.5,
        "axes.titlesize": 8,
        "legend.fontsize": 7,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "mathtext.fontset": "custom",
        "mathtext.rm": "Arial",
        "mathtext.it": "Arial:italic",
        "mathtext.bf": "Arial:bold",
        "mathtext.default": "it",
        "axes.linewidth": 0.6,
        "lines.linewidth": 1.0,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })

    res = json.loads((out / f"design_calculator_results_{tag}.json").read_text())
    att = res["attainability"]["grid"]
    thr = res["attainability"]["thresholds"][str(M_FAMILY)]
    agg = res["power"]["cells"]

    def band(key):
        ns = list(N_GRID)
        g = {s: {a["n_per_group"]: a[key] for a in agg
                 if a["sweep"] == "band" and a["scenario"] == s}
             for s, _, _ in SCENARIOS}
        return (ns,
                [g["calibrated"][n] for n in ns],
                [g["adverse"][n] for n in ns],
                [g["favourable"][n] for n in ns])

    fig, axs = plt.subplots(2, 2, figsize=(183 / 25.4, 132 / 25.4))
    fig.subplots_adjust(left=0.085, right=0.985, top=0.925, bottom=0.075,
                        hspace=0.42, wspace=0.255)

    def panel(ax, letter, title):
        ax.set_title(f"{letter})  {title}", loc="left", fontweight="bold",
                     fontsize=8, pad=5)

    # ---- a: k* vs n -------------------------------------------------------
    ax = axs[0, 0]
    panel(ax, "a", "Local discoveries are unattainable below a design size")
    mcols = {500: OKABE["skyblue"], 923: OKABE["green"],
             M_FAMILY: OKABE["blue"], 5000: OKABE["purple"]}
    for m in M_GRID:
        g = sorted([r for r in att if r["m"] == m], key=lambda r: r["n_per_group"])
        xs = [r["n_per_group"] for r in g]
        ys = [r["k_star"] for r in g]
        bad = [r["k_star_exceeds_family"] for r in g]
        ax.plot(xs, ys, "-", color=mcols[m], lw=1.0, zorder=2,
                label=f"m = {m:,}".replace(",", "\u2009"))
        ax.plot([x for x, b in zip(xs, bad) if not b],
                [y for y, b in zip(ys, bad) if not b],
                "o", color=mcols[m], ms=3.0, mec="none", zorder=3)
        ax.plot([x for x, b in zip(xs, bad) if b],
                [y for y, b in zip(ys, bad) if b],
                "x", color=mcols[m], ms=4.0, mew=0.9, zorder=3)
    ax.axhline(1.0, color=OKABE["black"], ls="--", lw=0.7, zorder=4)
    k9 = thr["k_star_at_9v9"]
    ax.plot([9], [k9], "*", color=OKABE["black"], ms=9, mec="none", zorder=6)
    ax.annotate(f"this cohort\n9v9, m = 2\u2009769\nk* = {k9}",
                xy=(9.15, k9), xytext=(10.4, 55), fontsize=7, ha="left",
                va="center",
                arrowprops=dict(arrowstyle="-", lw=0.5, color=OKABE["black"]))
    ax.text(15.3, 1.12, "k* = 1",
            fontsize=7, ha="right", va="bottom", color=OKABE["grey"])
    ax.set_yscale("log")
    ax.set_xlabel("Animals per group (balanced, exact permutation)")
    ax.set_ylabel("k*  (minimum size of any BY discovery list)")
    ax.set_xlim(3.6, 15.5)
    ax.set_xticks([4, 6, 8, 10, 12, 15])
    h, l = ax.get_legend_handles_labels()
    h.append(Line2D([], [], ls="none", marker="x", color=OKABE["grey"],
                    ms=4.0, mew=0.9))
    l.append("k* > m: no rejection possible")
    ax.legend(h, l, frameon=False, loc="upper right", handlelength=1.3,
              borderpad=0.1, labelspacing=0.22, fontsize=7)

    # ---- b: P(at least one discovery) ------------------------------------
    ax = axs[0, 1]
    panel(ax, "b", "Chance of any local discovery (S = 80)")
    ns, cal, adv, fav = band("prob_at_least_one")
    ax.fill_between(ns, adv, fav, color=OKABE["blue"], alpha=0.16, lw=0,
                    zorder=1, label="θ 0.010-0.028, depth 70-280")
    ax.plot(ns, cal, "-o", color=OKABE["blue"], ms=3.2, mec="none", zorder=3,
            label="calibrated: θ 0.017, depth 140")
    gx = [a["n_per_group"] for a in agg if a["sweep"] == "grid"]
    gy = [a["prob_at_least_one"] for a in agg if a["sweep"] == "grid"]
    ax.plot(gx, gy, ".", color=OKABE["grey"], ms=3.0, zorder=2,
            label="full θ × depth grid (n = 9, 10)")
    ax.axhline(0.95, color=OKABE["black"], ls="--", lw=0.7, zorder=4)
    ax.text(5.75, 0.958, "0.95", fontsize=7, ha="left", va="bottom")
    ax.set_ylim(-0.03, 1.06)
    ax.set_xlim(5.6, 15.4)
    ax.set_xticks(list(N_GRID))
    ax.set_xlabel("Animals per group")
    ax.set_ylabel("P(≥ 1 BY discovery)")
    ax.legend(frameon=False, loc="lower right", handlelength=1.3, borderpad=0.1,
              labelspacing=0.22, fontsize=7)

    # ---- c: per-effect power ---------------------------------------------
    ax = axs[1, 0]
    panel(ax, "c", "Per-effect power at the calibrated effect size")
    ns, cal, adv, fav = band("per_effect_power")
    ax.fill_between(ns, adv, fav, color=OKABE["vermillion"], alpha=0.16, lw=0,
                    zorder=1, label="θ, depth range")
    ax.plot(ns, cal, "-o", color=OKABE["vermillion"], ms=3.2, mec="none",
            zorder=4, label="TV 17.4 pp (median)")
    for tv, c, ls in ((10.0, OKABE["skyblue"], "--"), (25.0, OKABE["purple"], "-.")):
        sub = sorted([a for a in agg if a["sweep"] == "tv"
                      and a["scenario"] == f"tv{tv:g}"],
                     key=lambda a: a["n_per_group"])
        ax.plot([a["n_per_group"] for a in sub],
                [a["per_effect_power"] for a in sub], ls, color=c, lw=1.0,
                zorder=3, label=f"TV {tv:g} pp")
    ax.axhline(0.80, color=OKABE["black"], ls="--", lw=0.7, zorder=5)
    ax.text(15.3, 0.815, "80% power", fontsize=7, ha="right", va="bottom")
    ax.set_ylim(-0.03, 1.06)
    ax.set_xlim(5.6, 15.4)
    ax.set_xticks(list(N_GRID))
    ax.set_xlabel("Animals per group")
    ax.set_ylabel("Per-effect power (BY, m = 2 769)")
    leg = ax.legend(frameon=True, loc="upper left", handlelength=1.5, borderpad=0.25,
                    labelspacing=0.22, bbox_to_anchor=(0.0, 0.70), fancybox=False)
    leg.get_frame().set_facecolor("white"); leg.get_frame().set_edgecolor("none")
    leg.get_frame().set_alpha(0.9); leg.set_zorder(10)

    # ---- d: expected discoveries ------------------------------------------
    ax = axs[1, 1]
    panel(ax, "d", "Expected number of BY discoveries")
    ns, cal, adv, fav = band("expected_discoveries")
    ax.fill_between(ns, np.maximum(adv, 1e-2), np.maximum(fav, 1e-2),
                    color=OKABE["green"], alpha=0.16, lw=0, zorder=1,
                    label="θ 0.010-0.028, depth 70-280")
    ax.plot(ns, np.maximum(cal, 1e-2), "-o", color=OKABE["green"], ms=3.2,
            mec="none", zorder=3, label="calibrated: θ 0.017, depth 140")
    sr = np.mean([a["s_realized_mean"] for a in agg if a["sweep"] == "band"])
    ax.axhline(sr, color=OKABE["black"], ls=":", lw=0.7, zorder=4)
    ax.text(5.55, sr * 1.10, f"realized true effects, mean {sr:.1f}",
            fontsize=7, ha="left", va="bottom")
    ks = {n: k_star_from_rows(M_FAMILY, signs_rows(n))["k_star"] for n in N_GRID}
    for i, n in enumerate(N_GRID):
        ax.annotate(f"k* = {ks[n]}", xy=(n, 1500 if i % 2 == 0 else 450), ha="center",
                    va="center", color=OKABE["grey"])
    ax.set_yscale("log")
    ax.set_ylim(1.5e-2, 4e3)
    ax.set_xlim(4.9, 15.6)
    ax.set_xticks(list(N_GRID))
    ax.set_xlabel("Animals per group")
    ax.set_ylabel("Expected BY discoveries")
    ax.legend(frameon=False, loc="lower right", handlelength=1.3, borderpad=0.1,
              labelspacing=0.22, fontsize=7)

    for ext in ("pdf", "png"):
        fig.savefig(out / f"Figure_11_design.{ext}",
                    dpi=600 if ext == "png" else None, facecolor="white")
    plt.close(fig)
    print(f"[figure] Figure_11_design.pdf/.png written to {out}")


def signs_rows(n: int) -> int:
    """Allocation-level count actually used for design size n."""
    return n_rows_exact(n) if math.comb(2 * n, n) <= EXACT_MAX else R_MC


# --------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=60)
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--tag", default="main")
    ap.add_argument("--bench", action="store_true")
    ap.add_argument("--figure", action="store_true")
    ap.add_argument("--no-figure", action="store_true")
    ap.add_argument("--outdir", default=str(Path(__file__).resolve().parent))
    args = ap.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    hm = harmonic(M_FAMILY)

    if args.figure:
        make_figure(out, args.tag)
        return 0

    if args.bench:
        ss = np.random.SeedSequence(1)
        for n in N_GRID:
            t0 = time.time()
            signs_for(n)
            t_sign = time.time() - t0
            t0 = time.time()
            run_replicate(ss.spawn(1)[0], n, TV_CAL, THETA_CAL, DEPTH_CAL, hm)
            print(f"n={n:3d} rows={signs_rows(n):>7,} "
                  f"signs={t_sign:6.2f}s rep={time.time() - t0:6.3f}s "
                  f"rss={resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6:7.1f} MB")
        return 0

    t0 = time.time()
    att = attainability_table()
    thr = design_thresholds(att)
    with (out / "attainability_grid.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(att[0].keys()))
        w.writeheader()
        w.writerows(att)
    ref = next(r for r in att if r["m"] == M_FAMILY and r["n_per_group"] == 9)
    assert ref["k_star"] == 20, ref            # agrees with Figure 9 / D_simulation
    assert thr[str(M_FAMILY)]["n_min_single_BY"] == 12, thr
    assert thr[str(M_FAMILY)]["n_min_single_Holm"] == 10, thr
    assert R_MC >= thr[str(M_FAMILY)]["permutation_budget_R_min"], R_MC
    print(f"[attainability] {len(att)} cells; k*(2769, 9v9) = 20 verified; "
          f"n_BY = 12, n_Holm = 10, R_min = "
          f"{thr[str(M_FAMILY)]['permutation_budget_R_min']:,}")

    cells = build_cells()
    root = np.random.SeedSequence(args.seed)
    cell_seqs = root.spawn(len(cells))
    records, agg = [], []
    for ci, (cell, seq) in enumerate(zip(cells, cell_seqs)):
        kids = seq.spawn(args.reps)
        recs = []
        for ri, kid in enumerate(kids):
            r = run_replicate(kid, cell["n_per_group"], cell["tv_pp"],
                              cell["theta"], cell["depth_median"], hm)
            r["sweep"] = cell["sweep"]
            r["scenario"] = cell["scenario"]
            r["replicate"] = ri
            recs.append(r)
        records.extend(recs)
        a = aggregate(cell, recs)
        agg.append(a)
        print(f"[{ci + 1:3d}/{len(cells)}] {cell['sweep']:5s} {cell['scenario']:14s} "
              f"n={cell['n_per_group']:2d} k*={a['k_star']:4d} "
              f"E[D]={a['expected_discoveries']:7.2f} "
              f"P(>=1)={a['prob_at_least_one']:.3f} "
              f"power={a['per_effect_power']:.3f} "
              f"FDR={a['fdr_mean_fdp']:.4f} [{time.time() - t0:.0f}s]",
              flush=True)

    rec_fields = list(records[0].keys())
    with (out / f"power_replicates_{args.tag}.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=rec_fields)
        w.writeheader()
        w.writerows(records)
    with (out / f"power_cells_{args.tag}.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(agg[0].keys()))
        w.writeheader()
        w.writerows(agg)

    rec_table = recommend(agg, thr)
    with (out / f"design_recommendation_{args.tag}.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rec_table[0].keys()))
        w.writeheader()
        w.writerows(rec_table)

    summary = {
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "seed": args.seed,
        "replicates_per_cell": args.reps,
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "platform": platform.platform(),
        "peak_rss_mb": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6, 1),
        "design": {
            "m_family": M_FAMILY,
            "n_clusters": N_CLUSTERS,
            "n_alternatives": N_ALT,
            "alpha": ALPHA,
            "s_designated": S_DESIGNATED,
            "tv_calibrated_pp": TV_CAL,
            "theta_calibrated": THETA_CAL,
            "depth_calibrated": DEPTH_CAL,
            "depth_logsd": DEPTH_LOGSD,
            "depth_min": DEPTH_MIN,
            "exact_enumeration_max_allocations": EXACT_MAX,
            "monte_carlo_allocations": R_MC,
            "n_grid": list(N_GRID),
            "theta_grid": list(THETA_GRID),
            "depth_grid": list(DEPTH_GRID),
            "tv_grid": list(TV_GRID),
            "scenarios": [{"name": s, "theta": t, "depth": d}
                          for s, t, d in SCENARIOS],
            "allocation_levels_by_n": {str(n): signs_rows(n) for n in N_GRID},
            "perm_mode_by_n": {str(n): ("exact" if math.comb(2 * n, n) <= EXACT_MAX
                                        else "monte_carlo") for n in N_GRID},
        },
        "attainability": {"grid": att, "thresholds": thr},
        "power": {"cells": agg},
        "recommendation": rec_table,
    }
    summary["runtime_seconds"] = round(time.time() - t0, 1)
    txt = json.dumps(summary, indent=2, sort_keys=True)
    (out / f"design_calculator_results_{args.tag}.json").write_text(txt)
    print(f"[done] {summary['runtime_seconds']}s  peak RSS "
          f"{summary['peak_rss_mb']} MB  sha256(json)="
          f"{hashlib.sha256(txt.encode()).hexdigest()[:16]}")
    if not args.no_figure:
        make_figure(out, args.tag)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
