"""Building the reader's input from the released PreFold-dG pickles.

The three blocks are the published ones.  Each is a weighted mean over the
interface pairs Boltz-2's own distogram puts within max(12 A, the 1st percentile
of interchain distances), weighted by 1/d^2:

    wsi  (384, 384)  outer product of the input single representation, s_inputs
    ws   (384, 384)  the same pooling applied to the trunk single representation
    wz   (128,)      the pair representation itself

The reader is fitted on the mutant-minus-wild-type difference of each block.
That difference is the only place the mutation enters: both complexes are folded
independently and pooled the same way.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

P_MAT = 384 * 384
P_WZ = 128
BLOCKS = ("wsi", "ws", "wz")


def load_pickles(mut_path: Path, wt_path: Path) -> tuple[dict, dict]:
    """Index the two released pickles by `subpath`.

    Each file is `(records, labels)`; every record is a dict carrying `subpath`,
    `wsi`, `ws`, `wz` and Boltz confidence scalars.
    """
    import joblib

    mut, _ = joblib.load(mut_path)
    wt, _ = joblib.load(wt_path)
    return {r["subpath"]: r for r in mut}, {r["subpath"]: r for r in wt}


def difference_blocks(mut: dict, wt: dict, mut_ids: list[str], wt_ids: list[str],
                      dtype=np.float32) -> dict[str, np.ndarray]:
    """Stack mutant-minus-wild-type differences, one array per block.

    Memory: `wsi` and `ws` are 147,456 wide, so the pair costs about 6.8 GB in
    float32 for the full 5,731-variant cohort.
    """
    n = len(mut_ids)
    if len(wt_ids) != n:
        raise ValueError("mut_ids and wt_ids must be the same length")
    out = {"wsi": np.empty((n, P_MAT), dtype=dtype),
           "ws": np.empty((n, P_MAT), dtype=dtype),
           "wz": np.empty((n, P_WZ), dtype=dtype)}
    missing = []
    for i, (s, w) in enumerate(zip(mut_ids, wt_ids)):
        if s not in mut or w not in wt:
            missing.append((s, w))
            continue
        for b in ("wsi", "ws"):
            out[b][i] = (np.asarray(mut[s][b], dtype=np.float32)
                         - np.asarray(wt[w][b], dtype=np.float32)).ravel()
        out["wz"][i] = (np.asarray(mut[s]["wz"], dtype=np.float32)
                        - np.asarray(wt[w]["wz"], dtype=np.float32)).ravel()
    if missing:
        raise KeyError(f"{len(missing)} variants absent from the pickles, first {missing[:3]}")
    return out
