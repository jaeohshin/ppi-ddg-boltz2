"""Scoring, and the error bars that go with it.

The headline metric is per-structure Pearson correlation: a correlation inside
each complex, averaged over complexes.  It only means anything on the full
cohort — the macro-average is over the 116 complexes that carry at least ten
variants, and a subset ablation does not have them.  Report pcc / srcc / mae on
subsets instead.

Every difference between two models is quoted with a paired bootstrap over
complexes, not over variants.  Resampling variants holds the complex set fixed
and understates the spread; in this project every point estimate believed
without a complex-level interval was later overturned.
"""
from __future__ import annotations

import numpy as np
from scipy import stats

MIN_GROUP = 10


def per_structure_pcc(pred: np.ndarray, y: np.ndarray, groups: np.ndarray,
                      min_group: int = MIN_GROUP) -> tuple[float, int]:
    """Mean within-complex Pearson r over complexes with >= `min_group` variants."""
    vals = []
    for g in np.unique(groups):
        m = groups == g
        if m.sum() >= min_group:
            r = stats.pearsonr(pred[m], y[m])[0]
            if np.isfinite(r):
                vals.append(r)
    return (float(np.mean(vals)) if vals else float("nan")), len(vals)


def per_structure_srcc(pred: np.ndarray, y: np.ndarray, groups: np.ndarray,
                       min_group: int = MIN_GROUP) -> tuple[float, int]:
    vals = []
    for g in np.unique(groups):
        m = groups == g
        if m.sum() >= min_group:
            r = stats.spearmanr(pred[m], y[m])[0]
            if np.isfinite(r):
                vals.append(r)
    return (float(np.mean(vals)) if vals else float("nan")), len(vals)


def paired_bootstrap(pred_a: np.ndarray, pred_b: np.ndarray, y: np.ndarray,
                     groups: np.ndarray, n_boot: int = 10_000, seed: int = 0,
                     min_group: int = MIN_GROUP) -> dict:
    """95% interval for (score of a) - (score of b), resampling whole complexes.

    Both models are scored on the same resampled complexes in each draw, so the
    interval is for the difference and not for either score on its own.
    """
    rng = np.random.default_rng(seed)
    per_g = {}
    for g in np.unique(groups):
        m = groups == g
        if m.sum() >= min_group:
            ra = stats.pearsonr(pred_a[m], y[m])[0]
            rb = stats.pearsonr(pred_b[m], y[m])[0]
            if np.isfinite(ra) and np.isfinite(rb):
                per_g[g] = (ra, rb)
    if not per_g:
        return {"delta": float("nan"), "lo": float("nan"), "hi": float("nan"), "n_groups": 0}
    A = np.array([v[0] for v in per_g.values()])
    B = np.array([v[1] for v in per_g.values()])
    n = len(A)
    idx = rng.integers(0, n, size=(n_boot, n))
    deltas = A[idx].mean(axis=1) - B[idx].mean(axis=1)
    lo, hi = np.percentile(deltas, [2.5, 97.5])
    return {"delta": float(A.mean() - B.mean()), "lo": float(lo), "hi": float(hi),
            "n_groups": n, "a": float(A.mean()), "b": float(B.mean())}
