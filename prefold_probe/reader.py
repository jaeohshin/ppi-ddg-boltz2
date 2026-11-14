"""The reader: kernel ridge with one penalty scale per feature block.

A PreFold-dG feature vector is three blocks glued together — `wsi` and `ws`,
147,456 numbers each, and `wz`, 128.  Fitting one ridge penalty to the
concatenation lets the two large matrices decide the solution and buries `wz`,
which is where most of the signal is.  Scaling each block's Gram matrix by its
own width fixes that, and it is the whole difference between a reader that
scores about 0.53 and one that scores 0.59 to 0.62.

The fit is closed form.  There is no training loop, no seed, no epoch count and
no early stopping, so a given input and penalty produce one answer.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

CHUNK = 16_384  # columns per pass, to keep the standardisation out of core


def _block_gram(X: np.ndarray, train: np.ndarray) -> np.ndarray:
    """This block's Gram alone, standardised on training rows and width-normalised.

    Dividing by the block's own width is what puts a 128-wide block and a
    147,456-wide one on the same footing; it must happen to this block's Gram
    only, before anything else is added to it.  Accumulated in float64: the
    column means of these representations are large relative to their spread,
    and a float32 accumulation moves the third decimal of the score.
    """
    n, p = X.shape
    G = np.zeros((n, n), dtype=np.float64)
    for c0 in range(0, p, CHUNK):
        blk = X[:, c0:c0 + CHUNK].astype(np.float64, copy=False)
        mu = blk[train].mean(axis=0, keepdims=True)
        sd = blk[train].std(axis=0, ddof=0, keepdims=True)
        np.maximum(sd, 1e-8, out=sd)
        z = (blk - mu) / sd
        G += z @ z.T
    G /= p
    return G


def kernel(blocks: dict[str, np.ndarray], train: np.ndarray) -> np.ndarray:
    """Sum of per-block, width-normalised Gram matrices.

    `train` is a boolean mask over rows; only those rows set the standardisation,
    so the test rows never inform the fit.
    """
    n = next(iter(blocks.values())).shape[0]
    K = np.zeros((n, n), dtype=np.float64)
    for name, X in blocks.items():
        if X.shape[0] != n:
            raise ValueError(f"block {name!r} has {X.shape[0]} rows, expected {n}")
        K += _block_gram(X, train)
    return K


@dataclass
class BlockRidge:
    """Dual-form ridge over a precomputed kernel, solved once for many penalties.

    One eigendecomposition of the training kernel serves the whole penalty grid,
    which is what makes a 13-point sweep cost the same as a single fit.
    """

    alphas: tuple[float, ...] = tuple(np.logspace(-3, 3, 13))
    _eig: tuple[np.ndarray, np.ndarray] | None = field(default=None, repr=False)
    _vty: np.ndarray | None = field(default=None, repr=False)
    _ymean: float = field(default=0.0, repr=False)

    def fit(self, K_train: np.ndarray, y_train: np.ndarray) -> "BlockRidge":
        self._ymean = float(np.mean(y_train))
        w, V = np.linalg.eigh(K_train)
        self._eig = (w, V)
        self._vty = V.T @ (y_train - self._ymean)
        return self

    def predict(self, K_test_train: np.ndarray, alpha: float) -> np.ndarray:
        if self._eig is None:
            raise RuntimeError("fit() first")
        w, V = self._eig
        return K_test_train @ (V @ (self._vty / (w + alpha))) + self._ymean

    def predict_grid(self, K_test_train: np.ndarray) -> dict[float, np.ndarray]:
        """Predictions for every penalty in `alphas`, from the one decomposition."""
        return {a: self.predict(K_test_train, a) for a in self.alphas}
