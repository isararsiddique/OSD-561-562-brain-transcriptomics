#!/usr/bin/env python3
"""Discrete p-value floor and Benjamini-Yekutieli discovery counts in small
exact-permutation compositional designs.

Self-contained, seeded simulation written for Paper 2. It has two parts.

PART A (analytic, exact arithmetic)
  For family size m and balanced design size n vs n, the two-sided exact
  allocation p value cannot fall below p_min = 2 / C(2n, n).  The BY adjusted
  value of a block of k tests tied at that floor is

        q(k) = p_min * m * H_m / k ,        H_m = sum_{i=1..m} 1/i

  which is decreasing in k, so BY can reject nothing at level alpha unless the
  family contains at least

        k* = ceil(p_min * m * H_m / alpha)

  tests at the floor.  k* is a property of the design and the family size only:
  it does not depend on effect size.  The script tabulates k* over a grid of
  m and n and checks the m = 2769, n = 9 cell against the manuscript.

PART B (Monte Carlo)
  Dirichlet-multinomial compositional counts for 18 animals (9 vs 9), three
  alternatives per cluster, a fixed set of S true effects held IDENTICAL across
  read-depth conditions, and exact enumeration of all C(18,9) = 48,620 label
  allocations.  Each replicate is evaluated at full depth and on the two
  complementary halves obtained by binomial thinning (A ~ Binom(N, 1/2),
  B = N - A), which is how the manuscript's read splits were built.  Reported
  per condition: the BY discovery count, the number of floor-level ties, the
  probability of zero discoveries, the split-half cosine of the native effect
  vectors, and the aggregate (global) permutation test p value.

  Only S true clusters plus N_PERM_NULL null clusters are permuted.  The
  remaining tests of the m = 2769 family are accounted for analytically: under
  label exchangeability the exact two-sided permutation p value is uniform on
  the discrete support {2/48620, 4/48620, ...}, so their p values are drawn
  from that exact null distribution rather than simulated from counts.

Usage
  python3 sim_by_floor.py --reps 200                 # full run
  python3 sim_by_floor.py --reps 20 --tag check      # reproducibility re-run
  python3 sim_by_floor.py --calibrate                # effect/depth grid probe
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
import sys
import time
from itertools import combinations
from pathlib import Path

import numpy as np

# --------------------------------------------------------------------------
# Fixed design constants (taken from the manuscript, not tuned)
# --------------------------------------------------------------------------
M_FAMILY = 2769          # 923 clusters x 3 contrasts
N_CLUSTERS = 923         # clusters in the simulated contrast
N_PER_GROUP = 9          # FLT n = 9 vs GC n = 9
N_ALT = 3                # alternatives per cluster
ALPHA = 0.05
THETA = 0.017            # Dirichlet-multinomial overdispersion
N_PERM_NULL = 500        # null clusters permuted explicitly
S_GRID = (20, 40, 80)    # number of true effects
SEED = 20261015

# Calibrated against the real cohort (see --calibrate and D_summary.md):
# the real FLT-GC BY candidates have median total variation 17.4 pp (q75 19.6),
# and real per-animal selected-alternative depths reconstruct to a median near
# 80-160 reads.  These values are chosen so that the half-depth floor-tie count
# lands on the k* = 20 boundary; they are NOT tuned to any outcome below.
TV_PP = 18.0             # total-variation separation of the two group means, pp
DEPTH_MEDIAN = 140       # median selected-alternative reads per animal-cluster
DEPTH_LOGSD = 0.70
DEPTH_MIN = 20           # manuscript coverage rule

CHUNK = 1500             # allocation rows per matmul block (bounds RAM < ~150 MB)
TIE_TOL = 1e-12          # relative tie tolerance, as in the manuscript

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
# Exact BY arithmetic
# --------------------------------------------------------------------------
def harmonic(m: int) -> float:
    """H_m = sum_{i=1..m} 1/i, summed ascending in float64."""
    return float(np.sum(1.0 / np.arange(1, m + 1, dtype=np.float64)))


def p_floor(n: int) -> float:
    """Smallest attainable two-sided exact allocation p for n vs n."""
    return 2.0 / math.comb(2 * n, n)


def k_star(m: int, n: int, alpha: float = ALPHA) -> dict:
    """Minimum number of floor-level ties needed for any BY rejection."""
    hm = harmonic(m)
    pmin = p_floor(n)
    raw = pmin * m * hm / alpha
    k = math.ceil(raw)
    attainable = k <= m
    return {
        "m": m,
        "n_per_group": n,
        "n_allocations": math.comb(2 * n, n),
        "harmonic_H_m": hm,
        "p_min": pmin,
        "k_star_real": raw,
        "k_star": k,
        "q_at_k_star": pmin * m * hm / k if k > 0 else float("nan"),
        "q_at_k_star_minus_1": pmin * m * hm / (k - 1) if k > 1 else float("nan"),
        "attainable_within_family": bool(attainable),
    }


def by_discoveries(p_sorted: np.ndarray, m: int, hm: float, alpha: float = ALPHA) -> int:
    """Benjamini-Yekutieli step-up count on an ascending p vector of length m."""
    thresh = np.arange(1, m + 1, dtype=np.float64) * alpha / (m * hm)
    ok = np.nonzero(p_sorted <= thresh)[0]
    return int(ok[-1] + 1) if ok.size else 0


# --------------------------------------------------------------------------
# Allocation enumeration
# --------------------------------------------------------------------------
def sign_matrix(n: int) -> np.ndarray:
    """All balanced +-1 label vectors of length 2n, one per mirror pair.

    Row 0 is the observed labelling (+1 for animals 0..n-1).  The statistic
    ||sum_i s_i z_i||^2 is invariant to a global sign flip, so enumerating
    C(2n-1, n-1) rows covers all C(2n, n) allocations exactly twice.
    """
    total = 2 * n
    combos = list(combinations(range(1, total), n - 1))
    assert len(combos) == math.comb(total - 1, n - 1)
    s = np.full((len(combos), total), -1.0, dtype=np.float64)
    s[:, 0] = 1.0
    for row, cols in enumerate(combos):
        s[row, list(cols)] = 1.0
    assert np.all(s.sum(axis=1) == 0)
    assert np.array_equal(s[0], np.r_[np.ones(n), -np.ones(n)])
    return s


def permute_exact(signs: np.ndarray, z: np.ndarray, n_feat: int, n_groups: int = 1,
                  chunk: int = CHUNK):
    """Exact two-sided allocation p values plus one aggregate test per group.

    z has shape (2n, 3 * n_feat): per feature the within-animal sqrt-composition
    vector used by the manuscript's equal-cluster kernel.  The per-feature
    statistic is T_f(s) = ||sum_i s_i z_fi||^2, whose ranking over allocations
    is identical to the squared distance between transformed group means.  The
    n_feat features are read as n_groups contiguous blocks (here: full depth,
    half A, half B); a block's aggregate statistic sums T_f over it, which is
    monotone equivalent to the manuscript's cluster-averaged kernel statistic.
    """
    assert n_feat % n_groups == 0
    per = n_feat // n_groups
    n_rows = signs.shape[0]
    counts = np.zeros(n_feat, dtype=np.int64)
    g_counts = np.zeros(n_groups, dtype=np.int64)
    t_obs = None
    g_obs = None
    for start in range(0, n_rows, chunk):
        blk = signs[start : start + chunk]
        v = (blk @ z).reshape(blk.shape[0], n_feat, N_ALT)
        t = np.einsum("rfj,rfj->rf", v, v)
        g = t.reshape(blk.shape[0], n_groups, per).sum(axis=2)
        if t_obs is None:                      # row 0 is the observed labelling
            t_obs = t[0].copy()
            g_obs = g[0].copy()
        counts += np.count_nonzero(t >= t_obs * (1.0 - TIE_TOL), axis=0)
        g_counts += np.count_nonzero(g >= g_obs * (1.0 - TIE_TOL), axis=0)
        del v, t, g
    # p = 2 * counts / C(2n, n); the mirror pair makes the floor 2/C(2n, n).
    return counts / float(n_rows), g_counts / float(n_rows)


# --------------------------------------------------------------------------
# Data generation
# --------------------------------------------------------------------------
def base_compositions(rng: np.random.Generator, n_feat: int) -> np.ndarray:
    """Pooled mean compositions obeying the manuscript's selection bounds."""
    pi = rng.dirichlet(np.array([5.0, 3.0, 2.0]), size=n_feat)
    pi = np.clip(pi, 0.06, 0.88)
    return pi / pi.sum(axis=1, keepdims=True)


def split_group_means(rng: np.random.Generator, pi: np.ndarray, tv_pp: float):
    """Sum-zero separation of two group means with total variation tv_pp."""
    d = rng.normal(size=pi.shape)
    d -= d.mean(axis=1, keepdims=True)
    l1 = np.abs(d).sum(axis=1, keepdims=True)
    d = d / l1 * (2.0 * tv_pp / 100.0)         # TV = 0.5 * sum|delta|
    a = pi + d / 2.0
    b = pi - d / 2.0
    bad = (a.min(axis=1) < 0.02) | (b.min(axis=1) < 0.02)
    a[bad] = pi[bad]
    b[bad] = pi[bad]
    return a, b, bad


def draw_counts(rng: np.random.Generator, mean_a: np.ndarray, mean_b: np.ndarray,
                n: int) -> np.ndarray:
    """Dirichlet-multinomial counts, shape (2n, n_feat, 3).

    Animal compositions are Dirichlet with concentration (1-theta)/theta, so the
    marginal variance carries overdispersion theta relative to a multinomial;
    reads are then multinomial at the animal-cluster depth.
    """
    n_feat = mean_a.shape[0]
    conc = (1.0 - THETA) / THETA
    out = np.empty((2 * n, n_feat, N_ALT), dtype=np.int64)
    depth = np.maximum(
        DEPTH_MIN,
        np.round(DEPTH_MEDIAN * np.exp(rng.normal(0.0, DEPTH_LOGSD, size=(2 * n, n_feat)))),
    ).astype(np.int64)
    for i in range(2 * n):
        means = mean_a if i < n else mean_b
        p = rng.standard_gamma(means * conc)       # (n_feat, 3) independent gammas
        p /= p.sum(axis=1, keepdims=True)          # -> Dirichlet(means * conc)
        out[i] = rng.multinomial(depth[i], p)
    return out


def fractions_to_z(counts: np.ndarray) -> np.ndarray:
    """Within-animal sqrt compositions; shape (2n, n_feat, 3) -> same shape."""
    tot = counts.sum(axis=2, keepdims=True).astype(np.float64)
    tot[tot == 0] = 1.0
    return np.sqrt(counts / tot)


def native_effect(z: np.ndarray, n: int) -> np.ndarray:
    """Mean transformed FLT vector minus mean transformed GC vector."""
    return z[:n].mean(axis=0) - z[n:].mean(axis=0)


def cosines(d1: np.ndarray, d2: np.ndarray) -> np.ndarray:
    num = np.einsum("fj,fj->f", d1, d2)
    den = np.linalg.norm(d1, axis=1) * np.linalg.norm(d2, axis=1)
    den[den == 0] = np.nan
    return num / den


# --------------------------------------------------------------------------
# One replicate
# --------------------------------------------------------------------------
def run_replicate(seed_seq: np.random.SeedSequence, s_true: int, signs: np.ndarray,
                  hm: float, n: int = N_PER_GROUP) -> dict:
    rng = np.random.default_rng(seed_seq)
    n_feat = s_true + N_PERM_NULL
    pi = base_compositions(rng, n_feat)

    mean_a = pi.copy()
    mean_b = pi.copy()
    is_true = np.zeros(n_feat, dtype=bool)
    is_true[:s_true] = True
    ta, tb, dropped = split_group_means(rng, pi[:s_true], TV_PP)
    mean_a[:s_true] = ta
    mean_b[:s_true] = tb
    is_true[:s_true] &= ~dropped

    counts_full = draw_counts(rng, mean_a, mean_b, n)
    half_a = rng.binomial(counts_full, 0.5)
    half_b = counts_full - half_a

    rec = {
        "S": s_true,
        "n_true_effects": int(is_true.sum()),
        "rms_fraction_diff_pp": float(
            np.sqrt(np.mean((mean_a[is_true] - mean_b[is_true]) ** 2)) * 100.0
        ),
        "tv_pp": float(np.mean(0.5 * np.abs(mean_a[is_true] - mean_b[is_true]).sum(1)) * 100.0),
        "median_depth_full": float(np.median(counts_full.sum(axis=2))),
    }

    zs, effects = {}, {}
    conditions = ("full", "half_A", "half_B")
    for name, cnt in (("full", counts_full), ("half_A", half_a), ("half_B", half_b)):
        z = fractions_to_z(cnt)
        zs[name] = z
        effects[name] = native_effect(z, n)

    # One sweep over all 24,310 mirror-reduced allocations covers the three
    # depth conditions for both the local family and the aggregate test.
    stacked = np.concatenate(
        [zs[k].reshape(2 * n, n_feat * N_ALT) for k in conditions], axis=1
    )
    p_all, p_glob = permute_exact(signs, stacked, n_feat * 3, n_groups=3)
    n_alloc_half = signs.shape[0]

    pmin = 1.0 / n_alloc_half
    n_analytic = M_FAMILY - n_feat
    for i, name in enumerate(conditions):
        p_loc = p_all[i * n_feat : (i + 1) * n_feat]
        p_extra = rng.integers(1, n_alloc_half + 1, size=n_analytic) / float(n_alloc_half)
        p_family = np.sort(np.concatenate([p_loc, p_extra]))
        k = by_discoveries(p_family, M_FAMILY, hm)
        at_floor = int(np.count_nonzero(p_loc <= pmin * (1.0 + TIE_TOL)))
        at_floor_extra = int(np.count_nonzero(p_extra <= pmin * (1.0 + TIE_TOL)))
        rec[f"{name}_discoveries"] = k
        rec[f"{name}_floor_ties"] = at_floor + at_floor_extra
        rec[f"{name}_floor_ties_true"] = int(
            np.count_nonzero(p_loc[is_true] <= pmin * (1.0 + TIE_TOL))
        )
        rec[f"{name}_p_global"] = float(p_glob[i])
        rec[f"{name}_min_p"] = float(p_family[0])
        # Composition of the rejected set, and a null-calibration check.
        if k > 0:
            cut = p_family[k - 1]
            n_true_rej = int(np.count_nonzero(p_loc[is_true] <= cut))
            rec[f"{name}_discoveries_true"] = n_true_rej
            rec[f"{name}_discoveries_false"] = k - n_true_rej
        else:
            rec[f"{name}_discoveries_true"] = 0
            rec[f"{name}_discoveries_false"] = 0
        rec[f"{name}_null_frac_p_le_05"] = float(np.mean(p_loc[~is_true] <= 0.05))

    cos_ab = cosines(effects["half_A"], effects["half_B"])
    rec["cos_true_median"] = float(np.nanmedian(cos_ab[is_true]))
    rec["cos_true_min"] = float(np.nanmin(cos_ab[is_true]))
    rec["cos_true_pos_frac"] = float(np.mean(cos_ab[is_true] > 0))
    rec["cos_null_median"] = float(np.nanmedian(cos_ab[~is_true]))
    rec["_cos_true"] = cos_ab[is_true]
    rec["_cos_null"] = cos_ab[~is_true]
    return rec


# --------------------------------------------------------------------------
# Figure
# --------------------------------------------------------------------------
def make_figure(out: Path, tag: str = "main") -> None:
    """Figure 9: discrete p-value floor governs BY discovery counts."""
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
        "legend.fontsize": 6.5,
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

    res = json.loads((out / f"sim_by_floor_results_{tag}.json").read_text())
    with (out / f"mc_replicates_{tag}.csv").open() as fh:
        rows = list(csv.DictReader(fh))
    cz = np.load(out / f"cosines_{tag}.npz")
    cos_true, cos_null = cz["true"], cz["null"]

    def col(name, cast=float):
        return np.array([cast(r[name]) for r in rows])

    s_col = col("S", int)
    conds = ("full", "half_A", "half_B")
    cmap = {"full": OKABE["blue"], "half_A": OKABE["vermillion"], "half_B": OKABE["orange"]}
    clab = {"full": "Full depth", "half_A": "Half A", "half_B": "Half B (complement)"}
    disc = {c: col(f"{c}_discoveries", int) for c in conds}
    ties = {c: col(f"{c}_floor_ties", int) for c in conds}

    ref = res["analytic"]["reference_cell"]
    kstar, hm, pmin = ref["k_star"], ref["harmonic_H_m"], ref["p_min"]

    fig, axs = plt.subplots(2, 3, figsize=(183 / 25.4, 126 / 25.4))
    fig.subplots_adjust(left=0.060, right=0.990, top=0.905, bottom=0.085,
                        hspace=0.70, wspace=0.26)

    def panel(ax, letter, title):
        ax.set_title(f"{letter})  {title}", loc="left", fontweight="bold", fontsize=8, pad=5)

    # ---- a: analytic k* grid ------------------------------------------------
    ax = axs[0, 0]
    panel(ax, "a", "Ties needed for any BY rejection")
    grid = res["analytic"]["grid"]
    ms = sorted({g["m"] for g in grid})
    mcols = [OKABE["skyblue"], OKABE["green"], OKABE["blue"], OKABE["purple"]]
    for m, c in zip(ms, mcols):
        g = sorted([x for x in grid if x["m"] == m], key=lambda x: x["n_per_group"])
        xs = [x["n_per_group"] for x in g]
        ys = [x["k_star"] for x in g]
        okmask = [x["attainable_within_family"] for x in g]
        ax.plot(xs, ys, "-", color=c, lw=1.0, zorder=2,
                label=f"m = {m:,}".replace(",", " "))
        ax.plot([x for x, o in zip(xs, okmask) if o], [y for y, o in zip(ys, okmask) if o],
                "o", color=c, ms=3.2, mec="none", zorder=3)
        ax.plot([x for x, o in zip(xs, okmask) if not o],
                [y for y, o in zip(ys, okmask) if not o],
                "x", color=c, ms=4.2, mew=0.9, zorder=3)
    ax.plot([9], [kstar], "*", color=OKABE["black"], ms=9, zorder=5, mec="none")
    ax.annotate("this study\nm = 2 769, 9v9\nk* = 20", xy=(9.05, kstar * 1.25),
                xytext=(9.35, 700), fontsize=6.3, ha="left", va="center",
                arrowprops=dict(arrowstyle="-", lw=0.5, color=OKABE["black"]))
    ax.set_yscale("log")
    ax.set_xticks([5, 8, 9, 10])
    ax.set_xticklabels(["5v5", "8v8", "9v9", "10v10"])
    ax.set_xlabel("Exact design size (animals per group)")
    ax.set_ylabel("k*  (tests at the p floor)")
    ax.set_xlim(4.6, 11.0)
    h, l = ax.get_legend_handles_labels()
    h.append(Line2D([], [], ls="none", marker="x", color=OKABE["grey"], ms=4.2, mew=0.9))
    l.append("k* > m: no rejection possible")
    ax.legend(h, l, frameon=False, loc="lower left", handlelength=1.3,
              borderpad=0.1, labelspacing=0.25, fontsize=6.0)

    # ---- b: full discovery-count distributions ------------------------------
    ax = axs[0, 1]
    panel(ax, "b", "Discovery counts are bimodal")
    ax.axhspan(0.5, kstar - 0.5, color="#EDEDED", zorder=0)
    rng = np.random.default_rng(1)
    pos, xt, xl = [], [], []
    for si, s in enumerate(S_GRID):
        for ci, c in enumerate(conds):
            x0 = si * 3.4 + ci
            d = disc[c][s_col == s]
            ax.scatter(x0 + rng.uniform(-0.3, 0.3, d.size), d, s=1.6,
                       color=cmap[c], alpha=0.30, lw=0, zorder=2, rasterized=True)
            ax.plot([x0 - 0.42, x0 + 0.42], [np.median(d)] * 2, color=OKABE["black"],
                    lw=1.1, zorder=4, solid_capstyle="butt")
            ax.text(x0, 12.5, f"{100 * np.mean(d == 0):.0f}", fontsize=5.8,
                    color=OKABE["grey"], ha="center", va="center", zorder=5)
            pos.append(x0)
        xt.append(si * 3.4 + 1)
        xl.append(f"S = {s}")
    ax.axhline(kstar, color=OKABE["black"], ls="--", lw=0.7, zorder=3)
    ax.text(-0.72, kstar + 1.6, "k* = 20", fontsize=6.3, ha="left", zorder=5)
    ax.text(-0.72, 17.2, "P(zero discoveries), %", fontsize=5.8,
            color=OKABE["grey"], ha="left", va="center", zorder=5)
    ax.text(-0.72, 5.0, "1-19 discoveries\nare unreachable", fontsize=6.0,
            color=OKABE["grey"], ha="left", va="center", zorder=1)
    ax.set_xticks(xt)
    ax.set_xticklabels(xl)
    ax.set_xlim(-0.85, pos[-1] + 0.6)
    ax.set_ylim(-2.5, 56)
    ax.set_ylabel("BY discoveries (m = 2 769)")
    ax.set_xlabel("True effects present (same at every depth)")
    ax.legend(handles=[Line2D([], [], marker="o", ls="none", ms=3, color=cmap[c],
                              label=clab[c]) for c in conds],
              frameon=False, loc="upper left", handlelength=1.0, borderpad=0.1,
              labelspacing=0.22, fontsize=6.0, ncol=1)

    # ---- c: the BY boundary curve -------------------------------------------
    ax = axs[0, 2]
    panel(ax, "c", "The boundary sits at k* = 20")
    ylo, yhi = 0.016, 0.118
    kk = np.linspace(8, 46, 400)
    ax.plot(kk, pmin * M_FAMILY * hm / kk, color=OKABE["grey"], lw=1.2, zorder=4,
            label=r"$q(k)=p_{\min}\,m\,H_m/k$")
    ax.axhline(ALPHA, color=OKABE["black"], ls="--", lw=0.7, zorder=4,
               label="BY level 0.05")
    ax.axvline(kstar, color=OKABE["black"], ls=":", lw=0.7, zorder=4)
    ax.plot([kstar], [pmin * M_FAMILY * hm / kstar], "o", ms=3.5,
            color=OKABE["green"], zorder=6, mec="none")
    ax.annotate(f"k = 20: q = {pmin * M_FAMILY * hm / kstar:.5f}\n"
                f"k = 19: q = {pmin * M_FAMILY * hm / (kstar - 1):.5f}",
                xy=(kstar, ALPHA), xytext=(23.5, 0.086), fontsize=6.3,
                arrowprops=dict(arrowstyle="-", lw=0.5, color=OKABE["black"]))
    # Observed floor-tie densities, drawn in the lower third of the same axes
    # (no second y axis, so nothing is clipped at the figure edge).
    sel = s_col == 80
    band = (ylo, ylo + 0.33 * (yhi - ylo))
    dens_max = 0.0
    hists = {}
    for c in ("full", "half_A"):
        t = ties[c][sel]
        bins = np.arange(t.min() - 0.5, t.max() + 1.5)
        dn, be = np.histogram(t, bins=bins, density=True)
        hists[c] = (dn, be)
        dens_max = max(dens_max, dn.max())
    for c in ("full", "half_A"):
        dn, be = hists[c]
        y = band[0] + (band[1] - band[0]) * dn / dens_max
        ax.stairs(y, be, baseline=band[0], color=cmap[c], lw=0.9, fill=False, zorder=3)
        ax.stairs(y, be, baseline=band[0], color=cmap[c], lw=0, fill=True,
                  alpha=0.18, zorder=2)
    ax.text(45.5, band[1] + 0.004, "observed floor ties, S = 80\n(full depth / half A)",
            fontsize=6.0, ha="right", va="bottom", color=OKABE["grey"])
    ax.set_xlabel("Tests tied at the exact p floor, k")
    ax.set_ylabel("BY adjusted value of the floor block")
    ax.set_ylim(ylo, yhi)
    ax.set_xlim(8, 46)
    ax.legend(frameon=False, loc="upper right", handlelength=1.3, borderpad=0.1,
              labelspacing=0.25, fontsize=6.0)

    # ---- d: floor ties determine the count ----------------------------------
    ax = axs[1, 0]
    panel(ax, "d", "Floor ties set the discovery count")
    for c in conds:
        t, d = ties[c][sel], disc[c][sel]
        ax.scatter(t + rng.uniform(-0.22, 0.22, t.size),
                   d + rng.uniform(-0.22, 0.22, d.size), s=2.2, color=cmap[c],
                   alpha=0.45, lw=0, zorder=2, rasterized=True)
    lim = np.arange(0, 50)
    ax.plot(lim, np.where(lim >= kstar, lim, 0), color=OKABE["black"], lw=0.8, zorder=3)
    ax.axvline(kstar, color=OKABE["black"], ls=":", lw=0.7, zorder=3)
    ax.set_xlabel("Tests tied at the exact p floor")
    ax.set_ylabel("BY discoveries")
    ax.set_xlim(2, 46)
    ax.set_ylim(-3, 56)
    ax.text(0.97, 0.30, "step: reject k if k ≥ k*, else 0",
            transform=ax.transAxes, fontsize=6.2, ha="right", va="center")
    ax.text(0.97, 0.15, "above the step: a second\ndiscrete level also clears BY",
            transform=ax.transAxes, fontsize=6.0, ha="right", va="center",
            color=OKABE["grey"])

    # ---- e: effect alignment stays near 1 -----------------------------------
    ax = axs[1, 1]
    panel(ax, "e", "Effect direction stays stable")
    for v, c, lab in ((cos_true, OKABE["green"], "True effects"),
                      (cos_null, OKABE["grey"], "Null clusters")):
        v = np.sort(v[np.isfinite(v)])
        ax.plot(v, np.arange(1, v.size + 1) / v.size, color=c, lw=1.1, label=lab)
    al = res["split_half_alignment"]
    ax.axvline(0, color=OKABE["black"], ls="--", lw=0.6)
    ax.set_xlim(-1.03, 1.03)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("Cosine between complementary-half effect vectors")
    ax.set_ylabel("Cumulative fraction of instances")
    ax.text(0.04, 0.60, f"True effects, n = {al['n_true_effect_instances']:,}\n"
                        f"median {al['cosine_true_median']:.4f}\n"
                        f"5th percentile {al['cosine_true_q05']:.4f}\n"
                        f"{100 * al['cosine_true_positive_fraction']:.3f}% positive"
            .replace(",", " "),
            transform=ax.transAxes, fontsize=6.2, va="center")
    ax.legend(frameon=False, loc="upper left", handlelength=1.3, borderpad=0.1,
              labelspacing=0.25, fontsize=6.0)

    # ---- f: complementary halves disagree -----------------------------------
    ax = axs[1, 2]
    panel(ax, "f", "Halves disagree categorically")
    a, b = disc["half_A"][sel], disc["half_B"][sel]
    ax.scatter(a + rng.uniform(-0.5, 0.5, a.size), b + rng.uniform(-0.5, 0.5, b.size),
               s=3.0, color=OKABE["purple"], alpha=0.45, lw=0, zorder=2, rasterized=True)
    ax.axvline(kstar / 2, color=OKABE["black"], ls=":", lw=0.7)
    ax.axhline(kstar / 2, color=OKABE["black"], ls=":", lw=0.7)
    cp = res["monte_carlo"]["S=80"]["complementary_pairs"]
    ax.text(0.03, 0.97, f"both halves 0\n{100 * cp['prob_both_zero']:.1f}%",
            transform=ax.transAxes, fontsize=6.3, va="top", ha="left")
    ax.text(0.97, 0.97, f"both halves > 0\n{100 * cp['prob_both_positive']:.1f}%",
            transform=ax.transAxes, fontsize=6.3, va="top", ha="right")
    ax.text(0.97, 0.47, f"exactly one half\npositive: "
                        f"{100 * cp['prob_exactly_one_half_positive']:.1f}% of pairs\n"
                        f"mean gap {cp['mean_within_pair_gap']:.1f}, max "
                        f"{cp['max_within_pair_gap']}",
            transform=ax.transAxes, fontsize=6.3, ha="right", va="center",
            color=OKABE["vermillion"], fontweight="bold")
    gp = res["monte_carlo"]["S=80"]["half_A"]["global_power_holm3_0.05"]
    ax.text(0.5, 0.20, f"the aggregate Holm test is positive\nin {100 * gp:.0f}% of these same halves",
            transform=ax.transAxes, fontsize=6.2, ha="center", va="center",
            color=OKABE["blue"])
    ax.set_xlabel("BY discoveries, half A")
    ax.set_ylabel("BY discoveries, half B")
    ax.set_xlim(-3.5, 34)
    ax.set_ylim(-3.5, 34)

    for ext in ("pdf", "png"):
        fig.savefig(out / f"Figure_9_simulation.{ext}", dpi=600 if ext == "png" else None,
                    facecolor="white")
    plt.close(fig)
    print(f"[figure] Figure_9_simulation.pdf/.png written to {out}")


# --------------------------------------------------------------------------
# Drivers
# --------------------------------------------------------------------------
def analytic_table(out: Path) -> list:
    rows = []
    for m in (500, 1000, M_FAMILY, 5000):
        for n in (5, 8, 9, 10):
            rows.append(k_star(m, n))
    ref = k_star(M_FAMILY, N_PER_GROUP)
    assert ref["k_star"] == 20, ref
    assert abs(ref["harmonic_H_m"] - 8.50363774780702) < 1e-11, ref
    assert abs(ref["p_min"] - 2.0 / 48620) < 1e-15, ref
    assert abs(ref["q_at_k_star"] - 0.04842980856371377) < 1e-12, ref
    assert abs(ref["q_at_k_star_minus_1"] - 0.05097874585654081) < 1e-12, ref
    with (out / "analytic_floor_grid.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    return rows


def calibrate(signs, hm):
    print("tv_pp depth  floor_prob_full  floor_prob_half")
    global TV_PP, DEPTH_MEDIAN
    for tv in (10.0, 12.0, 14.0, 16.0):
        for dep in (70, 110, 160):
            TV_PP, DEPTH_MEDIAN = tv, dep
            ss = np.random.SeedSequence(99)
            full, half = [], []
            for child in ss.spawn(6):
                r = run_replicate(child, 40, signs, hm)
                nt = r["n_true_effects"]
                full.append(r["full_floor_ties_true"] / nt)
                half.append(0.5 * (r["half_A_floor_ties_true"] + r["half_B_floor_ties_true"]) / nt)
            print(f"{tv:5.1f} {dep:5d}  {np.mean(full):15.3f}  {np.mean(half):15.3f}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=200)
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--tag", default="main")
    ap.add_argument("--calibrate", action="store_true")
    ap.add_argument("--figure", action="store_true",
                    help="rebuild the figure from existing outputs and exit")
    ap.add_argument("--no-figure", action="store_true")
    ap.add_argument("--outdir", default=str(Path(__file__).resolve().parent))
    args = ap.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    if args.figure:
        make_figure(out, args.tag)
        return 0

    signs = sign_matrix(N_PER_GROUP)
    hm = harmonic(M_FAMILY)

    if args.calibrate:
        calibrate(signs, hm)
        return 0

    t0 = time.time()
    grid = analytic_table(out)
    print(f"[analytic] {len(grid)} cells; m=2769 n=9 -> k*=20 verified")

    root = np.random.SeedSequence(args.seed)
    per_s = dict(zip(S_GRID, root.spawn(len(S_GRID))))
    records, cos_true, cos_null = [], [], []
    for s_true in S_GRID:
        children = per_s[s_true].spawn(args.reps)
        for r_i, child in enumerate(children):
            rec = run_replicate(child, s_true, signs, hm)
            cos_true.append(rec.pop("_cos_true"))
            cos_null.append(rec.pop("_cos_null"))
            rec["replicate"] = r_i
            records.append(rec)
        done = [r for r in records if r["S"] == s_true]
        print(
            f"[mc] S={s_true:3d} reps={len(done)} "
            f"full={np.mean([r['full_discoveries'] for r in done]):6.2f} "
            f"A={np.mean([r['half_A_discoveries'] for r in done]):6.2f} "
            f"B={np.mean([r['half_B_discoveries'] for r in done]):6.2f} "
            f"[{time.time() - t0:.0f}s]"
        )

    fields = [k for k in records[0] if not k.startswith("_")]
    with (out / f"mc_replicates_{args.tag}.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows({k: r[k] for k in fields} for r in records)

    cos_true = np.concatenate(cos_true)
    cos_null = np.concatenate(cos_null)
    np.savez_compressed(out / f"cosines_{args.tag}.npz", true=cos_true, null=cos_null)

    summary = {
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "seed": args.seed,
        "replicates_per_S": args.reps,
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "platform": platform.platform(),
        "design": {
            "m_family": M_FAMILY,
            "n_clusters_simulated": N_CLUSTERS,
            "n_per_group": N_PER_GROUP,
            "n_allocations": math.comb(2 * N_PER_GROUP, N_PER_GROUP),
            "n_alternatives": N_ALT,
            "alpha": ALPHA,
            "dirichlet_multinomial_theta": THETA,
            "permuted_null_clusters": N_PERM_NULL,
            "analytic_null_tests": M_FAMILY - (min(S_GRID) + N_PERM_NULL),
            "tv_pp_target": TV_PP,
            "depth_median": DEPTH_MEDIAN,
            "depth_logsd": DEPTH_LOGSD,
            "depth_min": DEPTH_MIN,
        },
        "analytic": {
            "reference_cell": k_star(M_FAMILY, N_PER_GROUP),
            "grid": grid,
        },
        "monte_carlo": {},
        "effect_sizes": {
            "rms_fraction_diff_pp_mean": float(np.mean([r["rms_fraction_diff_pp"] for r in records])),
            "tv_pp_mean": float(np.mean([r["tv_pp"] for r in records])),
            "median_depth_full_mean": float(np.mean([r["median_depth_full"] for r in records])),
        },
        "split_half_alignment": {
            "n_true_effect_instances": int(cos_true.size),
            "cosine_true_median": float(np.nanmedian(cos_true)),
            "cosine_true_q05": float(np.nanpercentile(cos_true, 5)),
            "cosine_true_q95": float(np.nanpercentile(cos_true, 95)),
            "cosine_true_positive_fraction": float(np.mean(cos_true > 0)),
            "cosine_true_min": float(np.nanmin(cos_true)),
            "n_null_instances": int(cos_null.size),
            "cosine_null_median": float(np.nanmedian(cos_null)),
            "cosine_null_positive_fraction": float(np.mean(cos_null > 0)),
        },
    }
    for s_true in S_GRID:
        sub = [r for r in records if r["S"] == s_true]
        blk = {}
        for name in ("full", "half_A", "half_B"):
            d = np.array([r[f"{name}_discoveries"] for r in sub])
            ft = np.array([r[f"{name}_floor_ties"] for r in sub])
            pg = np.array([r[f"{name}_p_global"] for r in sub])
            dt = np.array([r[f"{name}_discoveries_true"] for r in sub])
            df_ = np.array([r[f"{name}_discoveries_false"] for r in sub])
            blk[name] = {
                "discoveries_mean": float(d.mean()),
                "discoveries_median": float(np.median(d)),
                "discoveries_min": int(d.min()),
                "discoveries_max": int(d.max()),
                "prob_zero_discoveries": float(np.mean(d == 0)),
                "prob_at_least_20": float(np.mean(d >= 20)),
                "discoveries_nonzero_mean": float(d[d > 0].mean()) if np.any(d > 0) else 0.0,
                "floor_ties_mean": float(ft.mean()),
                "floor_ties_min": int(ft.min()),
                "floor_ties_max": int(ft.max()),
                "discoveries_true_total": int(dt.sum()),
                "discoveries_false_total": int(df_.sum()),
                "false_discovery_fraction": (
                    float(df_.sum() / d.sum()) if d.sum() > 0 else 0.0
                ),
                "null_frac_p_le_05_mean": float(
                    np.mean([r[f"{name}_null_frac_p_le_05"] for r in sub])
                ),
                "p_global_median": float(np.median(pg)),
                "p_global_max": float(pg.max()),
                "global_power_holm3_0.05": float(np.mean(np.minimum(1.0, 3.0 * pg) <= ALPHA)),
                "global_power_raw_0.05": float(np.mean(pg <= ALPHA)),
            }
        d_a = np.array([r["half_A_discoveries"] for r in sub])
        d_b = np.array([r["half_B_discoveries"] for r in sub])
        blk["complementary_pairs"] = {
            "prob_exactly_one_half_positive": float(np.mean((d_a == 0) != (d_b == 0))),
            "prob_both_zero": float(np.mean((d_a == 0) & (d_b == 0))),
            "prob_both_positive": float(np.mean((d_a > 0) & (d_b > 0))),
            "max_within_pair_gap": int(np.max(np.abs(d_a - d_b))),
            "mean_within_pair_gap": float(np.mean(np.abs(d_a - d_b))),
        }
        summary["monte_carlo"][f"S={s_true}"] = blk

    summary["runtime_seconds"] = round(time.time() - t0, 1)
    txt = json.dumps(summary, indent=2, sort_keys=True)
    (out / f"sim_by_floor_results_{args.tag}.json").write_text(txt)
    print(f"[done] {summary['runtime_seconds']}s  sha256(json)="
          f"{hashlib.sha256(txt.encode()).hexdigest()[:16]}")
    if not args.no_figure:
        make_figure(out, args.tag)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
