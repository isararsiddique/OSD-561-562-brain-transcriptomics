"""Worker B check: is the low train/test loading cosine explained by sampling noise alone?

For each RRRM-2 ground deposit (OSD-561 cerebellum, OSD-562 hippocampus; 12 ground animals,
3 per age x collection-cohort cell, one library each):
  1. normalise as the package does (mean count >= 10, median-of-ratios, log2(c/sf + 1));
  2. reproduce the rPRG ground protocol effect (age-balanced ISS-T minus LAR, mean over panel genes);
  3. enumerate all 81 leave-one-animal-per-cell splits (train 8 / test 4) and compute the cosine
     between the training and test whole-genome cohort contrasts, and the same on Sensor31 genes;
  4. compute the cosine EXPECTED if training and test animals shared exactly the same true loading
     (perfect transfer), given per-gene residual variance: E[cos] ~= S / sqrt((S + N_tr)(S + N_te)),
     with S an unbiased estimate of true signal energy, N = sum of sampling variances of each contrast;
  5. confirm (4) by a parametric simulation that draws train and test data from ONE common loading.
Uses only ground animals; no flight data. Gene filter is applied on all 12 ground animals (the package
filters on training animals only), so values are close to, not identical with, the package's.
"""
import itertools
import json
import os

import numpy as np
import pandas as pd

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
P = (""
     "Paper_1_Complete_Submission_Package/Paper1_Strengthened_Review_20261015")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "noise_ceiling_results.json")
rng = np.random.default_rng(20261015)

man = pd.read_csv(P + "/results/sample_manifest.csv")
pan = pd.read_csv(P + "/results/panel_definitions.csv")
ann = pd.read_csv(P + "/inputs/gene_annotation.csv.gz")
sym2ens = ann.dropna(subset=["gene_symbol"]).groupby("gene_symbol")["ensembl_id"].apply(list).to_dict()


def panel_ids(name, index):
    ids = set()
    for s in pan.loc[pan.panel == name, "symbol"]:
        ids.update(sym2ens.get(s, []))
    return [g for g in index if g in ids]


def normalize(counts, min_mean=10):
    counts = counts.loc[counts.index.astype(str).str.match(r"^ENSMUSG\d+$")]
    c = counts.to_numpy(float)
    keep = c.mean(axis=1) >= min_mean
    c, genes = c[keep], counts.index[keep]
    pos = (c > 0).all(axis=1)
    lg = np.log(c[pos])
    sf = np.exp(np.median(lg - lg.mean(axis=1, keepdims=True), axis=0))
    return pd.DataFrame(np.log2(c / sf + 1), index=genes, columns=counts.columns)


def contrast(Y, cells, members):
    """Age-balanced ISS-T minus LAR contrast from the listed animals (columns)."""
    d = 0.0
    for age in ("12", "29"):
        a = [s for s in members if cells[s] == ("ISS_T", age)]
        b = [s for s in members if cells[s] == ("LAR", age)]
        d = d + (Y[a].mean(axis=1) - Y[b].mean(axis=1))
    return d / 2.0


def cos(a, b):
    return float(np.dot(a, b) / np.sqrt(np.dot(a, a) * np.dot(b, b)))


res = {}
for acc in ("OSD-561", "OSD-562"):
    m = man[(man.accession == acc) & (man.flight == "GC")].copy()
    m["age"] = m["age"].astype(str)
    counts = pd.read_csv(P + f"/inputs/{acc}_counts.csv.gz", index_col=0)
    counts = counts[m["sample"].tolist()]
    Y = normalize(counts)
    cells = dict(zip(m["sample"], zip(m["cohort"], m["age"])))
    samples = m["sample"].tolist()
    groups = {k: [s for s in samples if cells[s] == k] for k in
              [("ISS_T", "12"), ("ISS_T", "29"), ("LAR", "12"), ("LAR", "29")]}
    assert all(len(v) == 3 for v in groups.values()), groups

    # residual variance within cells (df = 12 - 4 = 8)
    R = pd.concat([Y[v].sub(Y[v].mean(axis=1), axis=0) for v in groups.values()], axis=1)
    s2 = (R ** 2).sum(axis=1) / (len(samples) - 4)

    def var_contrast(n_per_cell):
        # Var(d_g) = s2 * (1/4) * sum_age (1/n + 1/n) = s2 / n for equal n per cell
        return s2 / n_per_cell

    d_full = contrast(Y, cells, samples)
    rprg = panel_ids("rPRG", Y.index)
    sens = panel_ids("Sensor31", Y.index)

    out = {"n_genes": int(Y.shape[0]), "n_rPRG_genes_found": len(rprg), "n_Sensor31_genes_found": len(sens),
           "rPRG_ground_protocol_effect_mean_log2": float(d_full[rprg].mean())}

    for label, idx in (("whole_genome", Y.index), ("Sensor31", pd.Index(sens))):
        S = float((d_full[idx] ** 2 - var_contrast(3)[idx]).sum())
        N_tr, N_te = float(var_contrast(2)[idx].sum()), float(var_contrast(1)[idx].sum())
        expected = S / np.sqrt((S + N_tr) * (S + N_te)) if S > 0 else 0.0
        # empirical 81 splits
        obs = []
        for pick in itertools.product(*[range(3)] * 4):
            test = [groups[k][i] for k, i in zip(groups, pick)]
            train = [s for s in samples if s not in test]
            obs.append(cos(contrast(Y, cells, train)[idx].to_numpy(), contrast(Y, cells, test)[idx].to_numpy()))
        # parametric perfect-transfer simulation: one common loading delta (shrunk unbiasedly), Gaussian noise
        d0 = d_full[idx].to_numpy()
        shrink = np.sqrt(max(S, 0) / float((d0 ** 2).sum()))
        delta = d0 * shrink  # energy-matched (||delta||^2 = S) estimate of the common true loading
        sd = np.sqrt(s2[idx].to_numpy())
        sims = []
        for _ in range(400):
            tr = delta + rng.normal(0, sd * np.sqrt(1 / 2))
            te = delta + rng.normal(0, sd * np.sqrt(1 / 1))
            sims.append(cos(tr, te))
        out[label] = {"signal_energy_S": S, "noise_energy_train": N_tr, "noise_energy_test": N_te,
                      "signal_to_noise_test": S / N_te,
                      "observed_split_cosine_mean": float(np.mean(obs)),
                      "observed_split_cosine_q025_q975": [float(np.quantile(obs, .025)), float(np.quantile(obs, .975))],
                      "expected_cosine_if_identical_loading_analytic": float(expected),
                      "expected_cosine_if_identical_loading_simulated_mean": float(np.mean(sims)),
                      "expected_cosine_sim_q025_q975": [float(np.quantile(sims, .025)), float(np.quantile(sims, .975))]}
    res[acc] = out
    print(acc, json.dumps(out, indent=1))

with open(OUT, "w") as fh:
    json.dump(res, fh, indent=1)
print("wrote", OUT)
