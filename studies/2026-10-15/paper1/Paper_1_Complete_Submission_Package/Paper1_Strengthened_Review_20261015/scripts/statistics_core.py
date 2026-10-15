"""Transparent exploratory linear models and animal-restricted permutations.

Permutation p-values test label exchangeability conditional on listed blocks;
observational collection cohorts do not establish a randomized causal effect.
"""
import itertools
import math

import numpy as np
import pandas as pd
from scipy import stats


def bh(p):
    p = np.asarray(p, dtype=float)
    out = np.full(p.shape, np.nan)
    valid = np.flatnonzero(np.isfinite(p))
    order = valid[np.argsort(p[valid])]
    out[order] = np.minimum(1, np.minimum.accumulate(
        (p[order] * len(order) / np.arange(1, len(order) + 1))[::-1])[::-1])
    return out


def normalize(counts, min_mean=10):
    counts = counts.loc[counts.index.astype(str).str.match(r"^ENSMUSG\d+$")]
    if not counts.index.is_unique or not counts.columns.is_unique:
        raise ValueError("Duplicate gene or sample IDs")
    c = counts.to_numpy(float)
    if not np.isfinite(c).all() or (c < 0).any():
        raise ValueError("Counts must be finite and nonnegative")
    c = c[c.mean(axis=1) >= min_mean]
    genes = counts.index[counts.mean(axis=1) >= min_mean]
    positive = (c > 0).all(axis=1)
    if positive.sum() < 10:
        raise ValueError("Insufficient positive genes for median-ratio normalization")
    log = np.log(c[positive])
    sf = np.exp(np.median(log - log.mean(axis=1, keepdims=True), axis=0))
    return pd.DataFrame(np.log2(c / sf + 1), index=genes, columns=counts.columns), sf


class Contrast:
    def __init__(self, y, z=None):
        self.y = np.asarray(y, float)
        if self.y.ndim != 2 or not np.isfinite(self.y).all():
            raise ValueError("Expression must be a finite gene-by-sample matrix")
        n = self.y.shape[1]
        self.x0 = np.ones((n, 1)) if z is None else np.column_stack([np.ones(n), z])
        if not np.isfinite(self.x0).all():
            raise ValueError("Nuisance design must be finite")
        u, s, _ = np.linalg.svd(self.x0, full_matrices=False)
        tol = max(self.x0.shape) * np.finfo(float).eps * s[0]
        self.q = u[:, s > tol]
        self.df = n - self.q.shape[1] - 1
        if self.df <= 0:
            raise ValueError("No residual degrees of freedom")
        self.y0 = self.y - (self.y @ self.q) @ self.q.T
        self.syy = np.sum(self.y0**2, axis=1)

    def fit(self, group):
        g = np.asarray(group, float)
        if not np.isfinite(g).all() or len(g) != self.y.shape[1]:
            raise ValueError("Invalid contrast labels")
        gt = g - self.q @ (self.q.T @ g)
        sxx = gt @ gt
        if sxx < 1e-10:
            raise ValueError("Group contrast is not identifiable after nuisance adjustment")
        beta = (self.y0 @ gt) / sxx
        rss = np.maximum(self.syy - beta**2 * sxx, 0)
        se = np.sqrt(rss / self.df / sxx)
        t = np.divide(beta, se, out=np.where(beta == 0, 0., np.copysign(np.inf, beta)), where=se > 0)
        return beta, t

    def cluster_fit(self, group, animals):
        beta, _ = self.fit(group)
        g = np.asarray(group, float)
        gt = g - self.q @ (self.q.T @ g)
        rss_y = self.y0 - beta[:, None] * gt
        clusters, inv = np.unique(animals, return_inverse=True)
        n, k, m = self.y.shape[1], self.q.shape[1] + 1, len(clusters)
        if m < 3:
            raise ValueError("Fewer than three animals")
        meat = np.zeros(len(beta))
        for i in range(m):
            score = rss_y[:, inv == i] @ gt[inv == i]
            meat += score**2
        cr1 = (m / (m - 1)) * ((n - 1) / (n - k))
        se = np.sqrt(cr1 * meat) / (gt @ gt)
        t = np.divide(beta, se, out=np.where(beta == 0, 0., np.copysign(np.inf, beta)), where=se > 0)
        p = 2 * stats.t.sf(abs(t), m - 1)
        return beta, se, t, p, m - 1

    def permutation_betas(self, labels):
        g = np.asarray(labels, float)
        gt = g - self.q @ (self.q.T @ g)
        sxx = np.sum(gt**2, axis=0)
        if (sxx < 1e-10).any():
            raise ValueError("Permutation contains an unidentifiable contrast")
        return (self.y0 @ gt) / sxx


def animal_permutations(meta, group, blocks=("age",), max_exact=100000, seed=561562):
    """Enumerate all assignments including observed, or sample plus-one MC null.

    Validate labels and block membership agree across every region of an animal.
    Never drop duplicates before checking their consistency.
    """
    frame = meta[list(blocks) + ["animal"]].copy()
    frame["label"] = np.asarray(group, int)
    if not frame.label.isin([0, 1]).all():
        raise ValueError("Labels must be binary")
    for col in ["label", *blocks]:
        if (frame.groupby("animal")[col].nunique(dropna=False) > 1).any():
            raise ValueError(f"Inconsistent {col} within an animal")
    per = frame.drop_duplicates("animal").set_index("animal")
    strata = []
    grouper = list(blocks)
    groups = per.groupby(grouper, dropna=False, sort=True) if grouper else [(None, per)]
    for _, block in groups:
        names = sorted(block.index)
        strata.append((names, int(block.label.sum())))
    total = math.prod(math.comb(len(a), n) for a, n in strata)
    animals = frame.animal.to_numpy()
    if total <= max_exact:
        choices = [list(itertools.combinations(a, n)) for a, n in strata]
        selected = [set(itertools.chain.from_iterable(c)) for c in itertools.product(*choices)]
        method = "complete enumeration including observed"
    else:
        rng = np.random.default_rng(seed)
        selected = [set(itertools.chain.from_iterable(rng.choice(a, n, replace=False).tolist()
                                                    for a, n in strata)) for _ in range(max_exact)]
        method = "Monte Carlo plus-one"
    labels = np.column_stack([np.isin(animals, list(s)).astype(float) for s in selected])
    return labels, method, total


def permutation_p(observed, null, method):
    exceed = np.sum(np.abs(null) >= abs(observed) - 1e-12, axis=-1)
    n = null.shape[-1]
    return exceed / n if method.startswith("complete") else (exceed + 1) / (n + 1)


def learn_anchor(y, group, z, sensor_index):
    c = Contrast(y, z)
    beta, _ = c.fit(group)
    if np.linalg.norm(beta) < 1e-10:
        raise ValueError("Zero ground-anchor contrast")
    v = beta / np.linalg.norm(beta)
    w = np.zeros(len(beta))
    w[sensor_index] = beta[sensor_index]
    if np.linalg.norm(w) < 1e-10:
        raise ValueError("Zero sensor weights")
    w /= np.linalg.norm(w)
    score = w @ c.y0
    if score @ score < 1e-10:
        raise ValueError("Zero ground score variance")
    gamma = (c.y0 @ score) / (score @ score)
    return v, w, gamma


def bounded_gamma(w, gamma, max_gain=2.0):
    """Shrink the oblique component to cap spectral gain; preserve unit sensor removal."""
    if max_gain < 1 or not np.isclose(w @ w, 1):
        raise ValueError("Bound requires max_gain >= 1 and unit sensor weights")
    h = gamma - (w @ gamma) * w
    hnorm = np.linalg.norm(h)
    strength = min(1., np.sqrt(max_gain**2 - 1) / hnorm) if hnorm > 0 else 1.
    return w + strength * h


def transform(y, anchor, mode):
    v, w, gamma = anchor
    if mode == "projection":
        return y - np.outer(v, v @ y)
    if mode == "score":
        return y - np.outer(gamma, w @ y)
    if mode == "bounded_score":
        return y - np.outer(bounded_gamma(w, gamma), w @ y)
    if mode == "raw":
        return y.copy()
    raise ValueError(mode)


def lambda_diagnostic(t, df):
    # This is a dispersion diagnostic, never a biological truth criterion.
    p = np.clip(2 * stats.t.sf(abs(t), df), 1e-300, 1)
    return np.median(stats.chi2.isf(p, 1)) / stats.chi2.ppf(.5, 1)
