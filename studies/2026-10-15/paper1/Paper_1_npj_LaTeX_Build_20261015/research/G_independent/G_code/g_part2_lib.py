"""Shared helpers for the Part 2 diagnostic (same functions as g_part2.py, without the run)."""
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from g_common import OUT, PKG, load_deposit, panel_ids  # noqa: E402,F401

SENSOR = set(panel_ids()["Sensor31"])
SPLITS = pd.read_csv(f"{PKG}/extensions/ground_holdout/split_animal_roles.csv")
SETTINGS = {"CB": "OSD-561", "HPC": "OSD-562", "RR23": "OSD-525"}


def norm_nofilter(cnt, ref=None):
    base = cnt if ref is None else ref
    pos = (base > 0).all(axis=1)
    geo = np.exp(np.log(base.loc[pos]).mean(axis=1))
    sf = cnt.loc[pos].div(geo, axis=0).median(axis=0)
    return np.log2(cnt.div(sf, axis=1) + 1.0)


def norm(cnt, ref=None):
    base = cnt if ref is None else ref
    keep = base.mean(axis=1) >= 10
    return norm_nofilter(cnt.loc[keep], None if ref is None else ref.loc[keep])


def coef(Y, age, g):
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
    df = len(g) - Z.shape[1] - 1
    resid = Yt - np.outer(b, gt)
    s2 = (resid ** 2).sum(axis=1) / df if df > 0 else np.full(len(b), np.nan)
    return pd.Series(b, Y.index), pd.Series(s2, Y.index), c


def cos(a, b):
    return float(a @ b / np.sqrt((a @ a) * (b @ b)))


def setup(acc):
    cnt, am = load_deposit(acc)
    sp = SPLITS[SPLITS.accession == acc]
    lab = sp.drop_duplicates("animal").set_index("animal")["label"].astype(float)
    return cnt, am, sp, lab
