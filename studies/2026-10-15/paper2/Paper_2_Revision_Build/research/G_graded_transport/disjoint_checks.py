#!/usr/bin/env python3
"""Disjoint-animal check for panels c and e of Figure 10 (written after seeing the shared-animal correlations).

Contrasts that share a group share that group's sampling noise, which inflates their correlation across
clusters. Splitting the shared group into disjoint 4/5 halves (126 splits of 9 animals) removes the shared
noise. Two checks are run: c) flight is shared, so FLT_A-GC is compared with FLT_B-VIV; e) habitat control is
shared, so FLT-GC_A is compared with GC_B-VIV.
Per-cluster effect = native mean total variation, the quantity plotted in panel c.
"""
import itertools, json, pathlib
import numpy as np, pandas as pd
from scipy.stats import spearmanr

HERE = pathlib.Path(__file__).resolve().parent
T = HERE.parents[2] / "Paper_2_Complete_Submission_Package/extensions/multiway_splicing/tables"
comp = pd.read_csv(T / "animal_alternative_compositions.csv.gz")
grp = pd.read_csv(T / "analysis_animals.csv").set_index("sample").group
P = comp.pivot_table(index=["cluster_id", "alternative_index"], columns="sample", values="proportion").sort_index()
clusters = P.index.get_level_values(0).to_numpy(); uniq, inv = np.unique(clusters, return_inverse=True)
G = {g: [s for s in P.columns if grp[s] == g] for g in ["FLT", "GC", "VIV"]}
X = P.to_numpy()

def tv(a, b):  # per-cluster native total variation between two animal sets
    d = np.abs(X[:, [P.columns.get_loc(s) for s in a]].mean(1) - X[:, [P.columns.get_loc(s) for s in b]].mean(1))
    return 0.5 * np.bincount(inv, weights=d)

out = {}
shared_c = spearmanr(tv(G["FLT"], G["GC"]), tv(G["FLT"], G["VIV"])).statistic
shared_e = spearmanr(tv(G["FLT"], G["GC"]), tv(G["GC"], G["VIV"])).statistic
dc, de = [], []
for A in itertools.combinations(G["FLT"], 4):
    B = [s for s in G["FLT"] if s not in A]; dc.append(spearmanr(tv(list(A), G["GC"]), tv(B, G["VIV"])).statistic)
for A in itertools.combinations(G["GC"], 4):
    B = [s for s in G["GC"] if s not in A]; de.append(spearmanr(tv(G["FLT"], list(A)), tv(B, G["VIV"])).statistic)
# Precision control: noisier (low-depth) clusters show larger TV in every contrast. Partial Spearman on log depth
# and number of alternatives (rank residuals) separates shared precision from shared response.
loc = pd.read_csv(T / "local_exact_results_complete.csv.gz"); L = loc[loc.contrast == "FLT_vs_GC"].set_index("cluster_id").loc[uniq]
from scipy.stats import rankdata
Z = np.column_stack([np.ones(len(uniq)), rankdata(np.log(L.min_selected_reads)), rankdata(L.n_alternatives)])
def resid(v): r = rankdata(v); return r - Z @ np.linalg.lstsq(Z, r, rcond=None)[0]
def pcor(a, b): ra, rb = resid(a), resid(b); return float(ra @ rb / np.sqrt((ra @ ra) * (rb @ rb)))
pc = [pcor(tv(list(A), G["GC"]), tv([s for s in G["FLT"] if s not in A], G["VIV"])) for A in itertools.combinations(G["FLT"], 4)]
pe = [pcor(tv(G["FLT"], list(A)), tv([s for s in G["GC"] if s not in A], G["VIV"])) for A in itertools.combinations(G["GC"], 4)]
# noise-only reference: TV between two disjoint halves of the SAME group carries no group effect at all
pn = [pcor(tv(list(A), [s for s in G["GC"] if s not in A]), tv(G["FLT"], G["VIV"])) for A in itertools.combinations(G["GC"], 4)]
out_extra = dict(c_disjoint_partial_median=float(np.median(pc)), c_disjoint_partial_range=[float(min(pc)), float(max(pc))],
                 e_disjoint_partial_median=float(np.median(pe)), e_disjoint_partial_range=[float(min(pe)), float(max(pe))],
                 noise_only_GCsplit_vs_FLTVIV_partial_median=float(np.median(pn)), noise_only_range=[float(min(pn)), float(max(pn))])
out = dict(n_clusters=len(uniq), n_splits=len(dc),
           c_shared_FLT_rho=float(shared_c), c_disjoint_FLT_median_rho=float(np.median(dc)), c_disjoint_range=[float(min(dc)), float(max(dc))],
           e_shared_GC_rho_TV=float(shared_e), e_disjoint_GC_median_rho_TV=float(np.median(de)), e_disjoint_range=[float(min(de)), float(max(de))],
           **out_extra, note="Descriptive. Disjoint halves reuse the same animals across splits; no inferential claim.")
(HERE / "disjoint_animal_checks.json").write_text(json.dumps(out, indent=1)); print(json.dumps(out, indent=1))
