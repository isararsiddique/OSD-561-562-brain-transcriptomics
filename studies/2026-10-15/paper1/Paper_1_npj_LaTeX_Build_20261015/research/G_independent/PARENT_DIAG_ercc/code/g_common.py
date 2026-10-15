"""Worker G: independent re-implementation helpers (written from methods.tex only).

No code from the paper's package was read or imported.
"""
import itertools
import json
import os
import re

import numpy as np
import pandas as pd

PKG = (""
       "Paper_1_Complete_Submission_Package/Paper1_Strengthened_Review_20261015")
BUILD = (""
         "Paper_1_npj_LaTeX_Build_20261015")
OUT = "Paper_1_npj_LaTeX_Build_20261015/research/G_independent/PARENT_DIAG_ercc"  # PARENT DIAGNOSTIC: outputs redirected
PANELS = ["rPRG", "dPRG", "SRG", "HSP", "Glia_Complement", "BBB_Endothelium", "Circadian", "Sensor31"]


# ---------------------------------------------------------------- loading
def load_deposit(acc):
    """Return (animal-level raw counts genes x animals, animal metadata) for one deposit.

    Technical libraries are summed within Source Name (= animal) before anything else.
    """
    meta = pd.read_csv(f"{PKG}/inputs/{acc}_metadata.csv.gz")
    cnt = pd.read_csv(f"{PKG}/inputs/{acc}_counts.csv.gz", index_col=0)
    cnt = cnt.loc[cnt.index.astype(str).str.match(r"^ENSMUSG\d+$")]  # PARENT DIAGNOSTIC: drop ERCC spike-ins
    meta = meta[meta["sample"].isin(cnt.columns)].copy()
    lib2animal = dict(zip(meta["sample"], meta["Source Name"]))
    cnt = cnt[list(lib2animal)]
    cnt = cnt.T.groupby(pd.Series(lib2animal)).sum().T  # sum technical libraries per animal
    rows = []
    for animal, sub in meta.groupby("Source Name"):
        r = {"animal": animal}
        if acc in ("OSD-561", "OSD-562"):
            # flight / cohort decoded from the Source Name, e.g. RRRM-2_FL_ISS_12, RRRM-2_HGC_LAR_05
            m = re.match(r"RRRM-2_(FL|HGC)_(ISS|LAR)_\d+$", animal)
            r["flight"] = "FLT" if m.group(1) == "FL" else "GC"
            r["cohort"] = "ISS_T" if m.group(2) == "ISS" else "LAR"
            r["age"] = int(sub["Factor Value[Age]"].iloc[0])
        else:
            sf = sub["Factor Value[Spaceflight]"].iloc[0]
            r["flight"] = {"Space Flight": "FLT", "Ground Control": "HGC", "Vivarium Control": "VGC"}[sf]
            r["cohort"] = "LAR"
            r["age"] = 0
        r["n_libraries"] = len(sub)
        rows.append(r)
    am = pd.DataFrame(rows).set_index("animal").loc[cnt.columns]
    return cnt.astype(float), am


def filter_and_normalise(cnt, ref_cnt=None):
    """Mean raw count >= 10 filter, median-ratio size factors on genes positive in every sample,
    log2(count/sf + 1). If ref_cnt is given, the filter and the geometric-mean reference are
    computed from ref_cnt (training animals) and applied to cnt (test animals)."""
    base = cnt if ref_cnt is None else ref_cnt
    keep = base.mean(axis=1) >= 10
    b = base.loc[keep]
    pos = (b > 0).all(axis=1)
    logref = np.log(b.loc[pos]).mean(axis=1)  # log geometric mean
    c = cnt.loc[keep]
    ratios = np.log(c.loc[pos]).sub(logref, axis=0)
    sf = np.exp(ratios.median(axis=0))  # zeros in test give -inf; median is robust to a few
    expr = np.log2(c.div(sf, axis=1) + 1.0)
    return expr, sf


# ---------------------------------------------------------------- annotation / panels
def symbol_map():
    ann = pd.read_csv(f"{PKG}/inputs/gene_annotation.csv.gz")
    d = {}
    for e, s in zip(ann["ensembl_id"], ann["gene_symbol"]):
        if isinstance(s, str):
            d.setdefault(s, []).append(e)
    return d


def panel_ids():
    p = pd.read_csv(f"{PKG}/results/panel_definitions.csv")
    sm = symbol_map()
    out = {}
    for name, sub in p.groupby("panel"):
        ids = []
        for s in sub["symbol"]:
            ids += sm.get(s, [])
        out[name] = ids
    return out


# ---------------------------------------------------------------- estimators
def resid_maker(Z):
    """Return function v -> (I - H_Z) v."""
    Q, _ = np.linalg.qr(Z)
    return lambda v: v - Q @ (Q.T @ v)


def build_Z(am, pooled):
    cols = [np.ones(len(am))]
    if am["age"].nunique() > 1:
        cols.append((am["age"].values == am["age"].max()).astype(float))
    if pooled:
        cols.append((am["tissue"].values == "HPC").astype(float))
    return np.column_stack(cols)


def exact_panel_test(Y, am, g, panels, pooled):
    """Y: genes x samples (DataFrame), am: sample meta with 'animal_id','age'; g: 0/1 per sample
    (constant within animal). Effects for each panel = mean OLS coefficient over panel genes.
    Exact two-sided p over whole-animal relabelings preserving group sizes within age blocks."""
    Z = build_Z(am, pooled)
    R = resid_maker(Z)
    an = am.groupby("animal_id").agg(age=("age", "first"))
    an["g"] = pd.Series(g, index=am.index).groupby(am["animal_id"]).first()
    assert (pd.Series(g, index=am.index).groupby(am["animal_id"]).nunique() == 1).all()
    # enumerate assignments
    blocks = []
    for age, sub in an.groupby("age"):
        ids = list(sub.index)
        k = int(sub["g"].sum())
        blocks.append([set(c) for c in itertools.combinations(ids, k)])
    assigns = []
    for combo in itertools.product(*blocks):
        pos = set().union(*combo)
        assigns.append(np.array([1.0 if a in pos else 0.0 for a in am["animal_id"]]))
    G = np.column_stack(assigns)  # samples x A
    Gt = R(G)
    denom = (Gt * Gt).sum(axis=0)
    res = {}
    g_obs = np.asarray(g, float)
    gt_obs = R(g_obs)
    for name, ids in panels.items():
        ids = [i for i in ids if i in Y.index]
        ybar = Y.loc[ids].values.mean(axis=0)
        null = (ybar @ Gt) / denom
        obs = float(ybar @ gt_obs / (gt_obs @ gt_obs))
        p = float(np.mean(np.abs(null) >= abs(obs) - 1e-12))
        res[name] = {"effect": obs, "p_exact": p, "n_assign": len(assigns), "n_genes": len(ids)}
    return res


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(q, 1)
    return out


def dump(obj, name):
    with open(os.path.join(OUT, name), "w") as f:
        json.dump(obj, f, indent=1, default=float)
