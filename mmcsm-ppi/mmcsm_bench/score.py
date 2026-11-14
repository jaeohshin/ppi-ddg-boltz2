"""Scoring, with the interval that decides whether a win is a win.

Correlations are quoted with a bootstrap over whole complexes, not over
mutations: 43% of the training rows sit on ten complexes, and resampling
mutations inside a fixed complex set understates the spread badly.
"""
from __future__ import annotations

import numpy as np
from scipy import stats


def metrics(y, p) -> dict:
    y, p = np.asarray(y, float), np.asarray(p, float)
    k = np.isfinite(y) & np.isfinite(p)
    y, p = y[k], p[k]
    return {"n": int(len(y)),
            "pearson": float(stats.pearsonr(y, p)[0]),
            "spearman": float(stats.spearmanr(y, p)[0]),
            "kendall": float(stats.kendalltau(y, p)[0]),
            "rmse": float(np.sqrt(((y - p) ** 2).mean())),
            "mae": float(np.abs(y - p).mean())}


def bootstrap_complexes(y, p, groups, n_boot=4000, seed=0, stat="pearson") -> tuple[float, float]:
    """95% interval for one method's score, resampling complexes with replacement."""
    y, p, g = np.asarray(y, float), np.asarray(p, float), np.asarray(groups)
    rng = np.random.default_rng(seed)
    uniq = np.unique(g)
    idx = {q: np.where(g == q)[0] for q in uniq}
    f = stats.pearsonr if stat == "pearson" else stats.spearmanr
    vals = []
    for _ in range(n_boot):
        m = np.concatenate([idx[q] for q in rng.choice(uniq, len(uniq))])
        if len(np.unique(y[m])) > 2:
            vals.append(f(y[m], p[m])[0])
    return (float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))) if vals else (np.nan, np.nan)


def paired_delta(y, p_a, p_b, groups, n_boot=4000, seed=0, stat="pearson") -> dict:
    """Is a better than b on these rows? Both are scored on the same resampled complexes."""
    y = np.asarray(y, float); a = np.asarray(p_a, float); b = np.asarray(p_b, float)
    g = np.asarray(groups)
    k = np.isfinite(y) & np.isfinite(a) & np.isfinite(b)
    y, a, b, g = y[k], a[k], b[k], g[k]
    f = stats.pearsonr if stat == "pearson" else stats.spearmanr
    rng = np.random.default_rng(seed)
    uniq = np.unique(g)
    idx = {q: np.where(g == q)[0] for q in uniq}
    ds = []
    for _ in range(n_boot):
        m = np.concatenate([idx[q] for q in rng.choice(uniq, len(uniq))])
        if len(np.unique(y[m])) > 2:
            ds.append(f(y[m], a[m])[0] - f(y[m], b[m])[0])
    lo, hi = (np.percentile(ds, [2.5, 97.5]) if ds else (np.nan, np.nan))
    return {"delta": float(f(y, a)[0] - f(y, b)[0]), "lo": float(lo), "hi": float(hi),
            "n_complexes": int(len(uniq)), "beats": bool(lo > 0)}
