"""Part 2: identical-loading noise ceiling on ground-only holdouts (methods.tex eq. 3)."""
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from g_common import PKG, dump, load_deposit, panel_ids

SENSOR = set(panel_ids()["Sensor31"])
splits = pd.read_csv(f"{PKG}/extensions/ground_holdout/split_animal_roles.csv")


def norm(cnt, ref=None):
    """Filter (mean>=10) and median-ratio reference from `ref` (default cnt itself); apply to cnt.
    Ratios taken in linear space so a zero test count gives ratio 0 rather than -inf."""
    base = cnt if ref is None else ref
    keep = base.mean(axis=1) >= 10
    b = base.loc[keep]
    pos = (b > 0).all(axis=1)
    geo = np.exp(np.log(b.loc[pos]).mean(axis=1))
    c = cnt.loc[keep]
    sf = c.loc[pos].div(geo, axis=0).median(axis=0)
    return np.log2(c.div(sf, axis=1) + 1.0)


def coef(Y, age, g):
    """Per-gene OLS coefficient of g with intercept (+ age if it varies); returns coef, sigma2, c=g~'g~."""
    cols = [np.ones(len(g))]
    if len(np.unique(age)) > 1:
        cols.append((age == age.max()).astype(float))
    Z = np.column_stack(cols)
    Q, _ = np.linalg.qr(Z)
    gt = g - Q @ (Q.T @ g)
    c = gt @ gt
    Yv = Y.values
    Yt = Yv - (Yv @ Q) @ Q.T
    b = Yt @ gt / c
    X = np.column_stack([Z, g])
    df = len(g) - X.shape[1]
    resid = Yt - np.outer(b, gt)  # FWL: residual of full model
    s2 = (resid ** 2).sum(axis=1) / df if df > 0 else np.full(len(b), np.nan)
    return pd.Series(b, Y.index), pd.Series(s2, Y.index), c


def cos(a, b):
    return float(a @ b / np.sqrt((a @ a) * (b @ b)))


settings = {"CB": "OSD-561", "HPC": "OSD-562", "RR23": "OSD-525"}
out = {}
rng = np.random.default_rng(15102026)
for name, acc in settings.items():
    cnt, am = load_deposit(acc)
    sp = splits[splits.accession == acc]
    ground_animals = sorted(sp.animal.unique())
    lab = sp.drop_duplicates("animal").set_index("animal")["label"].astype(float)
    if acc == "OSD-525":
        # label 1 = habitat ground control (HGC), 0 = vivarium (VGC); check against metadata
        assert all((am.loc[a, "flight"] == "HGC") == (lab[a] == 1) for a in ground_animals)
    else:
        assert all((am.loc[a, "cohort"] == "ISS_T") == (lab[a] == 1) for a in ground_animals)
        assert all(am.loc[a, "flight"] == "GC" for a in ground_animals)
    gc = cnt[ground_animals]
    age = am.loc[ground_animals, "age"].values
    # ---- observed per-split cosines
    cw, cs = [], []
    for sid, sub in sp.groupby("split_id"):
        tr = list(sub[sub.role == "train"].animal)
        te = list(sub[sub.role == "test"].animal)
        Ytr = norm(cnt[tr])
        Yte = norm(cnt[te], ref=cnt[tr])
        btr, _, _ = coef(Ytr, am.loc[tr, "age"].values, lab[tr].values)
        bte, _, _ = coef(Yte, am.loc[te, "age"].values, lab[te].values)
        cw.append(cos(btr.values, bte.values))
        sg = [i for i in btr.index if i in SENSOR]
        cs.append(cos(btr[sg].values, bte[sg].values))
    # group-size constants of the split designs (constant across splits)
    sub0 = sp[sp.split_id == sp.split_id.min()]
    tr0, te0 = list(sub0[sub0.role == "train"].animal), list(sub0[sub0.role == "test"].animal)

    def cval(ids):
        _, _, c = coef(norm(cnt[ids]).iloc[:5], am.loc[ids, "age"].values, lab[ids].values)
        return c

    c_tr, c_te = cval(tr0), cval(te0)
    # ---- expected cosine from all ground animals of the deposit
    Yall = norm(gc)
    b, s2, c_all = coef(Yall, age, lab[ground_animals].values)
    res = {"n_ground_animals": len(ground_animals), "n_splits": int(sp.split_id.nunique()),
           "c_all": c_all, "c_tr": c_tr, "c_te": c_te,
           "observed_cos_whole_mean": float(np.mean(cw)), "observed_cos_whole_sd": float(np.std(cw)),
           "observed_cos_sensor_mean": float(np.mean(cs)), "observed_cos_sensor_sd": float(np.std(cs))}
    for gset, idx in (("whole", b.index), ("sensor", [i for i in b.index if i in SENSOR])):
        bb, ss = b[idx].values, s2[idx].values
        S = float(bb @ bb - ss.sum() / c_all)  # unbiased energy of the true contrast
        Ntr, Nte = float(ss.sum() / c_tr), float(ss.sum() / c_te)
        exp_cos = S / np.sqrt((S + Ntr) * (S + Nte))
        # 400 Gaussian simulations: one common loading with energy S, independent noise
        if S > 0:
            delta = bb * np.sqrt(S / (bb @ bb))
        else:
            delta = np.zeros_like(bb)
        sims = []
        for _ in range(400):
            dt = delta + rng.normal(0, np.sqrt(ss / c_tr))
            de = delta + rng.normal(0, np.sqrt(ss / c_te))
            sims.append(cos(dt, de))
        res[gset] = {"n_genes": len(idx), "S": S, "N_tr": Ntr, "N_te": Nte,
                     "expected_cos": float(exp_cos), "SNR_tr": S / Ntr,
                     "sim400_mean_cos": float(np.mean(sims))}
    out[name] = res
    print(name, res)

dump(out, "G_part2.json")
