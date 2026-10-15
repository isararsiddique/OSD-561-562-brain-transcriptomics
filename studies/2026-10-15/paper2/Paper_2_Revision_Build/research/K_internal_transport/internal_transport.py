#!/usr/bin/env python3
"""Internal (held-out RR23) transport benchmark; see PLAN_2026-10-15.md. Read-only on the package."""
import os
os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
import itertools, json, pathlib
import numpy as np, pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
PKG = HERE.parents[2] / "Paper_2_Complete_Submission_Package"
MS = PKG / "extensions/multiway_splicing/tables"; CT = PKG / "extensions/comparator_transport/results"
assert (HERE / "PLAN_2026-10-15.md").exists()

tr = pd.read_csv(CT / "RR23_only_frozen_directions.csv.gz"); tr = tr[tr.direction == "FLT_minus_GC"]
ids = tr.cluster_id.drop_duplicates().to_numpy()
comp = pd.read_csv(MS / "animal_alternative_compositions.csv.gz")
grp = pd.read_csv(MS / "analysis_animals.csv").set_index("sample").group
P = comp.pivot_table(index=["cluster_id", "alternative_index"], columns="sample", values="proportion")
P = P.loc[list(zip(tr.cluster_id, tr.alternative_index))]          # frozen alternative order
X = P.to_numpy(); S = np.array(P.columns)
F = np.flatnonzero([grp[s] == "FLT" for s in S]); G = np.flatnonzero([grp[s] == "GC" for s in S])
assert len(F) == 9 and len(G) == 9
off = np.r_[0, np.flatnonzero(np.diff(pd.factorize(tr.cluster_id)[0])) + 1, len(tr)]
sizes = np.diff(off)

# reproduce the frozen RR23 direction exactly (guards against any alignment error)
d_full = X[:, F].mean(1) - X[:, G].mean(1)
assert np.allclose(d_full, tr.RR23_native_delta.to_numpy(), atol=1e-12)

def unit(d):
    n = np.sqrt(np.add.reduceat(d * d, off[:-1])); n = np.repeat(n, sizes)
    return np.divide(d, n, out=np.zeros_like(d), where=n > 1e-12)

def score(u, dtest, mask):  # equal-cluster mean projection, in percentage points
    return 100 * np.add.reduceat(u * dtest, off[:-1])[mask].mean()

panels = {"all923": np.ones(len(ids), bool)}
for tag in ["RR10", "RRRM2"]:
    z = np.load(CT / f"{tag}_frozen_measurements.npz", allow_pickle=True); assert list(z["cluster_ids"]) == list(ids)
    panels[tag] = z["gate"]

def split_half(Fi, Gi, rng=None, nmax=None):
    out = {k: [] for k in panels}
    combos = [(a, b) for a in itertools.combinations(range(9), 5) for b in itertools.combinations(range(9), 5)]
    if nmax: combos = [combos[i] for i in rng.choice(len(combos), nmax, replace=False)]
    for a, b in combos:
        ftr, gtr = Fi[list(a)], Gi[list(b)]; fte, gte = np.setdiff1d(Fi, ftr), np.setdiff1d(Gi, gtr)
        u = unit(X[:, ftr].mean(1) - X[:, gtr].mean(1)); dt = X[:, fte].mean(1) - X[:, gte].mean(1)
        for k, m in panels.items(): out[k].append(score(u, dt, m))
    return {k: np.array(v) for k, v in out.items()}

sh = split_half(F, G)
lopo = {k: [] for k in panels}
for f in F:
    for g in G:
        u = unit(X[:, np.setdiff1d(F, f)].mean(1) - X[:, np.setdiff1d(G, g)].mean(1)); dt = X[:, f] - X[:, g]
        for k, m in panels.items(): lopo[k].append(score(u, dt, m))
lopo = {k: np.array(v) for k, v in lopo.items()}

# label-permuted reference: same procedure after shuffling FLT/GC labels within RR23
rng = np.random.default_rng(2026101599); perm = {k: [] for k in panels}; FG = np.r_[F, G]
for _ in range(2000):
    p = rng.permutation(FG); r = split_half(p[:9], p[9:], rng, nmax=50)
    for k in panels: perm[k].append(np.median(r[k]))
perm = {k: np.array(v) for k, v in perm.items()}

# ---- Amendment 1: same-date GC-vs-VIV control and matched per-split permutation reference
V = np.flatnonzero([grp[s] == "VIV" for s in S]); assert len(V) == 8
gv = {k: [] for k in panels}
for a in itertools.combinations(range(9), 5):
    for b in itertools.combinations(range(8), 4):
        gtr, vtr = G[list(a)], V[list(b)]; gte, vte = np.setdiff1d(G, gtr), np.setdiff1d(V, vtr)
        u = unit(X[:, gtr].mean(1) - X[:, vtr].mean(1)); dt = X[:, gte].mean(1) - X[:, vte].mean(1)
        for k, m in panels.items(): gv[k].append(score(u, dt, m))
gv = {k: np.array(v) for k, v in gv.items()}
rng2 = np.random.default_rng(2026101598); pooled = {k: [] for k in panels}
for _ in range(2000):
    p = rng2.permutation(FG); r = split_half(p[:9], p[9:], rng2, nmax=50)
    for k in panels: pooled[k].append(r[k])
pooled = {k: np.concatenate(v) for k, v in pooled.items()}

ext = pd.read_csv(CT / "all_comparator_global_transport.csv")
ext = ext[(ext.direction == "FLT_minus_GC") & (ext.panel == "cohort_measured")].set_index("cohort")
q = lambda a: [float(np.quantile(a, .025)), float(np.median(a)), float(np.quantile(a, .975))]
res = {}
for k in panels:
    res[k] = dict(n_clusters=int(panels[k].sum()), split_half_q025_median_q975=q(sh[k]), split_half_frac_positive=float((sh[k] > 0).mean()),
                  lopo_q025_median_q975=q(lopo[k]), lopo_frac_positive=float((lopo[k] > 0).mean()),
                  permuted_reference_median_q025_median_q975=q(perm[k]),
                  internal_vs_permuted_exceedance=float((perm[k] >= np.median(sh[k])).mean()))
    if k in ext.index:
        e = ext.loc[k]; res[k].update(external_effect_pp=float(e.effect_pp), external_CI95=[float(e.conditional_CI95_low_pp), float(e.conditional_CI95_high_pp)],
                                      external_over_internal_split_median=float(e.effect_pp / np.median(sh[k])),
                                      internal_split_frac_at_or_below_external=float((sh[k] <= e.effect_pp).mean()))
for k in panels:
    res[k].update(GC_VIV_same_date_split_q025_median_q975=q(gv[k]), GC_VIV_frac_positive=float((gv[k] > 0).mean()),
                  GC_VIV_over_FLT_GC_internal_median=float(np.median(gv[k]) / np.median(sh[k])),
                  pooled_perm_split_frac_ge_internal_median=float((pooled[k] >= np.median(sh[k])).mean()),
                  pooled_perm_split_q975=float(np.quantile(pooled[k], .975)))
    if "external_effect_pp" in res[k]:
        res[k]["external_point_over_internal_median"] = float(res[k]["external_effect_pp"] / np.median(sh[k]))
res["n_GC_VIV_splits"] = int(len(gv["all923"]))
res["n_split_half"] = int(len(sh["all923"])); res["n_lopo"] = int(len(lopo["all923"])); res["n_permutations"] = 2000
np.savez_compressed(HERE / "internal_transport_draws.npz", **{f"split_{k}": v for k, v in sh.items()}, **{f"lopo_{k}": v for k, v in lopo.items()}, **{f"perm_{k}": v for k, v in perm.items()}, **{f"gcviv_{k}": v for k, v in gv.items()})
(HERE / "summary.json").write_text(json.dumps(res, indent=1)); print(json.dumps(res, indent=1))
