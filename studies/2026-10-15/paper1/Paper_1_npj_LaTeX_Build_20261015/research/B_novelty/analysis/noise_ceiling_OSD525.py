"""Noise-ceiling check for RR-23 (OSD-525) ground animals: HGC vs VGC contrast.

Package ground-holdout design: 200 splits, train 6 HGC + 5 VGC, test 3 HGC + 3 VGC (17 ground animals).
Technical-replicate libraries are summed per animal before normalisation (package collapses by animal).
Expected train/test cosine under ONE shared loading: S / sqrt((S + N_tr)(S + N_te)).
"""
import json
import os

import numpy as np
import pandas as pd

P = (""
     "Paper_1_Complete_Submission_Package/Paper1_Strengthened_Review_20261015")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "noise_ceiling_OSD525.json")
rng = np.random.default_rng(52520261015)

man = pd.read_csv(P + "/results/sample_manifest.csv")
m = man[(man.accession == "OSD-525") & (man.cohort.isin(["HGC", "VGC"]))]
counts = pd.read_csv(P + "/inputs/OSD-525_counts.csv.gz", index_col=0)[m["sample"].tolist()]
counts = counts.T.groupby(m.set_index("sample")["animal"]).sum().T  # animal-level
coh = m.drop_duplicates("animal").set_index("animal")["cohort"]
counts = counts[coh.index]
counts = counts.loc[counts.index.astype(str).str.match(r"^ENSMUSG\d+$")]
c = counts.to_numpy(float)
keep = c.mean(axis=1) >= 10
c = c[keep]
pos = (c > 0).all(axis=1)
lg = np.log(c[pos])
sf = np.exp(np.median(lg - lg.mean(axis=1, keepdims=True), axis=0))
Y = np.log2(c / sf + 1)
A = np.flatnonzero(coh.values == "HGC")
B = np.flatnonzero(coh.values == "VGC")
print("animals HGC", len(A), "VGC", len(B), "genes", Y.shape[0])

resid = np.concatenate([Y[:, A] - Y[:, A].mean(1, keepdims=True), Y[:, B] - Y[:, B].mean(1, keepdims=True)], axis=1)
s2 = (resid ** 2).sum(1) / (len(A) + len(B) - 2)
d = Y[:, A].mean(1) - Y[:, B].mean(1)
S = float((d ** 2 - s2 * (1 / len(A) + 1 / len(B))).sum())
N_tr = float((s2 * (1 / 6 + 1 / 5)).sum())
N_te = float((s2 * (1 / 3 + 1 / 3)).sum())
expected = S / np.sqrt((S + N_tr) * (S + N_te))


def cosv(a, b):
    return float(a @ b / np.sqrt((a @ a) * (b @ b)))


obs = []
for _ in range(200):
    ta = rng.choice(A, 3, replace=False)
    tb = rng.choice(B, 3, replace=False)
    ra = np.setdiff1d(A, ta)
    rb = rng.choice(np.setdiff1d(B, tb), 5, replace=False) if len(B) - 3 > 5 else np.setdiff1d(B, tb)
    dtr = Y[:, ra].mean(1) - Y[:, rb].mean(1)
    dte = Y[:, ta].mean(1) - Y[:, tb].mean(1)
    obs.append(cosv(dtr, dte))
res = {"n_HGC": int(len(A)), "n_VGC": int(len(B)), "n_genes": int(Y.shape[0]), "signal_energy_S": S,
       "noise_energy_train": N_tr, "noise_energy_test": N_te, "signal_to_noise_test": S / N_te,
       "expected_cosine_if_identical_loading_analytic": float(expected),
       "observed_random_split_cosine_mean": float(np.mean(obs)),
       "observed_q025_q975": [float(np.quantile(obs, .025)), float(np.quantile(obs, .975))],
       "note": "whole-genome HGC minus VGC contrast; splits drawn here, not the package's exact 200"}
print(json.dumps(res, indent=1))
json.dump(res, open(OUT, "w"), indent=1)
