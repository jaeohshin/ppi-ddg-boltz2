"""The reader's arithmetic, checked where it can be checked without the data.

Three things matter and are testable on synthetic blocks: that the dual solve
agrees with the textbook primal ridge, that the block normalisation actually
equalises blocks of different widths, and that test rows never touch the
standardisation.
"""
import numpy as np
import pytest

from prefold_probe.evaluate import paired_bootstrap, per_structure_pcc
from prefold_probe.reader import BlockRidge, _block_gram, kernel


def _blocks(n=60, seed=0):
    rng = np.random.default_rng(seed)
    return {"wide": rng.normal(size=(n, 200)).astype(np.float32),
            "narrow": rng.normal(size=(n, 4)).astype(np.float32)}


def test_dual_matches_primal():
    """Kernel ridge on the standardised design equals primal ridge on it."""
    rng = np.random.default_rng(1)
    n, p = 40, 12
    X = rng.normal(size=(n, p))
    y = rng.normal(size=n)
    train = np.zeros(n, dtype=bool)
    train[:30] = True
    alpha = 2.5

    mu, sd = X[train].mean(0), X[train].std(0)
    Z = (X - mu) / sd / np.sqrt(p)          # the reader's scaling, in primal form
    K = Z @ Z.T
    dual = BlockRidge().fit(K[np.ix_(train, train)], y[train]).predict(
        K[np.ix_(~train, train)], alpha)

    ym = y[train].mean()
    A = Z[train].T @ Z[train] + alpha * np.eye(p)
    w = np.linalg.solve(A, Z[train].T @ (y[train] - ym))
    primal = Z[~train] @ w + ym
    assert np.allclose(dual, primal, atol=1e-8)


def test_block_normalisation_equalises_width():
    """A block's contribution must not grow with its number of columns."""
    rng = np.random.default_rng(2)
    n = 50
    train = np.ones(n, dtype=bool)
    narrow = rng.normal(size=(n, 10))
    wide = np.repeat(narrow, 20, axis=1)      # same information, 20x the columns
    g_n = _block_gram(narrow.astype(np.float32), train)
    g_w = _block_gram(wide.astype(np.float32), train)
    assert np.allclose(g_n, g_w, rtol=1e-4, atol=1e-4)


def test_blocks_accumulate_independently():
    """Each block is divided by its own width, not by a later block's."""
    b = _blocks()
    train = np.ones(b["wide"].shape[0], dtype=bool)
    K = kernel(b, train)
    expected = _block_gram(b["wide"], train) + _block_gram(b["narrow"], train)
    assert np.allclose(K, expected)


def test_standardisation_uses_training_rows_only():
    """Shifting a test row must not change what the training rows see."""
    b = _blocks()
    n = b["wide"].shape[0]
    train = np.zeros(n, dtype=bool)
    train[: n // 2] = True
    K1 = kernel(b, train)[np.ix_(train, train)]
    b2 = {k: v.copy() for k, v in b.items()}
    b2["wide"][-1] += 1000.0
    K2 = kernel(b2, train)[np.ix_(train, train)]
    assert np.allclose(K1, K2)


def test_grid_shares_one_decomposition():
    b = _blocks()
    n = b["wide"].shape[0]
    train = np.zeros(n, dtype=bool)
    train[:40] = True
    rng = np.random.default_rng(3)
    y = rng.normal(size=n)
    K = kernel(b, train)
    r = BlockRidge(alphas=(0.1, 1.0, 10.0)).fit(K[np.ix_(train, train)], y[train])
    grid = r.predict_grid(K[np.ix_(~train, train)])
    assert set(grid) == {0.1, 1.0, 10.0}
    for a, p in grid.items():
        assert np.allclose(p, r.predict(K[np.ix_(~train, train)], a))


def test_per_structure_pcc_skips_small_groups():
    rng = np.random.default_rng(4)
    y = rng.normal(size=30)
    pred = y + rng.normal(scale=0.1, size=30)
    groups = np.array(["big"] * 20 + ["small"] * 10)
    score, n = per_structure_pcc(pred, y, groups, min_group=15)
    assert n == 1 and score > 0.9


def test_paired_bootstrap_is_about_the_difference():
    rng = np.random.default_rng(5)
    y = rng.normal(size=120)
    groups = np.repeat([f"g{i}" for i in range(6)], 20)
    good = y + rng.normal(scale=0.3, size=120)
    bad = y + rng.normal(scale=3.0, size=120)
    out = paired_bootstrap(good, bad, y, groups, n_boot=2000, min_group=10)
    assert out["n_groups"] == 6
    assert out["lo"] > 0          # the better model is better, and the interval says so
    assert out["delta"] == pytest.approx(out["a"] - out["b"])
