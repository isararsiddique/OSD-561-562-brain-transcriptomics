#!/usr/bin/env python3
"""Evidence-graded transport (see PLAN_2026-10-15.md). Exact, label-allocation null.

Reproduces the published per-cluster external effects first (assert), then tests whether
they rise with RR23 discovery evidence. Read-only on the submission package.
"""
import os
os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
import itertools, json, pathlib
import numpy as np, pandas as pd
from scipy.stats import rankdata, spearmanr

HERE = pathlib.Path(__file__).resolve().parent
PKG = HERE.parents[2] / "Paper_2_Complete_Submission_Package"
CT = PKG / "extensions/comparator_transport/results"
MS = PKG / "extensions/multiway_splicing/tables"
assert (HERE / "PLAN_2026-10-15.md").exists()
TOL = 1e-12

def weights(meta):  # identical to comparator_transport/scripts/analyse.py
    observed = np.zeros(len(meta)); alloc = np.zeros((1, len(meta)))
    for st in sorted(meta.stratum.unique()):
        ix = np.flatnonzero(meta.stratum.to_numpy() == st); nf = int(meta.iloc[ix].exposure.eq("FLT").sum()); ng = len(ix) - nf
        w = len(ix) / len(meta); combos = list(itertools.combinations(ix, nf))
        D = np.zeros((len(combos), len(meta))); D[:, ix] = -w / ng
        for j, co in enumerate(combos): D[j, list(co)] = w / nf
        observed[ix] = np.where(meta.iloc[ix].exposure.eq("FLT"), w / nf, -w / ng)
        alloc = (alloc[:, None, :] + D[None, :, :]).reshape(-1, len(meta))
    return observed, alloc

def rank_rows(M):  # Spearman needs ranks of each allocation's effect vector
    return np.apply_along_axis(rankdata, 1, M)

def spearman_rows(x_rank, Y):  # rho of fixed ranks against each row of Y
    yr = rank_rows(Y); xr = x_rank - x_rank.mean(); yr = yr - yr.mean(1, keepdims=True)
    return (yr @ xr) / (np.sqrt((yr ** 2).sum(1)) * np.sqrt((xr ** 2).sum()))

def holm(p):
    p = np.asarray(p); o = np.argsort(p); m = len(p); adj = np.empty(m); run = 0
    for k, i in enumerate(o): run = max(run, (m - k) * p[i]); adj[i] = min(1, run)
    return adj

tr = pd.read_csv(CT / "RR23_only_frozen_directions.csv.gz")
loc = pd.read_csv(MS / "local_exact_results_complete.csv.gz")
pub = pd.read_csv(CT / "all3692_local_directional_tests.csv.gz")
rows, effects = [], {}
for tag in ["RR10", "RRRM2"]:
    z = np.load(CT / f"{tag}_frozen_measurements.npz", allow_pickle=True)
    meta = pd.read_csv(CT / f"{tag}_included_animals.csv").set_index("sample").loc[z["samples"]]
    if tag == "RR10": meta["exposure"] = meta.group; meta["stratum"] = "all"
    observed, alloc = weights(meta)
    off, gate, ids = z["offsets"], z["gate"], z["cluster_ids"]
    for direction, contrast in [("FLT_minus_GC", "FLT_vs_GC"), ("FLT_minus_VIV", "FLT_vs_VIV")]:
        d = tr[tr.direction == direction]; assert list(d.cluster_id.drop_duplicates()) == list(ids)
        u = d.RR23_unit_vector.to_numpy(); base = d.RR23_comparator_mean.to_numpy()
        scores = np.add.reduceat(u[:, None] * (z["proportions"] - base[:, None]), off[:-1], axis=0)
        eff = scores @ observed
        # reproduce the published per-cluster effects exactly
        ref = pub[(pub.cohort == tag) & (pub.direction == direction)].set_index("cluster_id").loc[ids]
        assert np.allclose(eff[gate], ref.effect.to_numpy()[gate], atol=1e-12), (tag, direction)
        ev = -np.log10(loc[loc.contrast == contrast].set_index("cluster_id").loc[ids].exact_p.to_numpy())
        g = gate; xr = rankdata(ev[g])
        null_eff = alloc @ scores[g].T                      # allocations x clusters
        rho_null = spearman_rows(xr, null_eff); rho_obs = spearmanr(ev[g], eff[g]).statistic
        p = float((rho_null >= rho_obs - TOL).mean())
        top = ev[g] >= np.quantile(ev[g], 0.9)
        diff_obs = eff[g][top].mean() - eff[g][~top].mean(); diff_null = null_eff[:, top].mean(1) - null_eff[:, ~top].mean(1)
        rows.append(dict(cohort=tag, direction=direction, n_clusters=int(g.sum()), n_allocations=len(alloc),
                         spearman_rho=float(rho_obs), exact_p_one_sided=p,
                         null_rho_q975=float(np.quantile(rho_null, .975)), null_rho_q025=float(np.quantile(rho_null, .025)),
                         top_decile_n=int(top.sum()), top_minus_rest_pp=100 * float(diff_obs),
                         top_minus_rest_exact_p=float((diff_null >= diff_obs - TOL).mean())))
        effects[(tag, direction)] = pd.Series(np.where(gate, eff, np.nan), index=ids)
res = pd.DataFrame(rows); res["Holm_p_four"] = holm(res.exact_p_one_sided)
res.to_csv(HERE / "graded_transport_results.csv", index=False)

# external-external concordance (descriptive)
a, b = effects[("RR10", "FLT_minus_GC")], effects[("RRRM2", "FLT_minus_GC")]
both = a.notna() & b.notna()
ee = spearmanr(a[both], b[both]).statistic
per = pd.DataFrame({f"{c}_{d}": s for (c, d), s in effects.items()})
per["RR23_FLT_vs_GC_exact_p"] = loc[loc.contrast == "FLT_vs_GC"].set_index("cluster_id").exact_p
per.to_csv(HERE / "per_cluster_external_effects.csv")
summary = dict(results=rows, holm=list(res.Holm_p_four), external_external_spearman_RR10_RRRM2=float(ee),
               n_common=int(both.sum()), note="retrospective; exact allocation nulls; see PLAN_2026-10-15.md")
(HERE / "summary.json").write_text(json.dumps(summary, indent=1))
print(res.to_string()); print("RR10 vs RRRM2 Spearman on", int(both.sum()), "common clusters:", round(ee, 4))
