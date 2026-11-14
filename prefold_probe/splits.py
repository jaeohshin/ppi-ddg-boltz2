"""Fold assignments.

Four splits ship with this package, all as `PDB,fold` over three folds.

  rdenet        what the field uses: grouped by PDB id (RDE-PPI / PreFold-dG).
                Bushuiev et al. measured 56% of its test interfaces as leaked.
  ppiformer     PPIRef's similarity-aware split.
  hom99         ours: GroupKFold over 99%-identity chain clusters. 611 distinct
                wild-type chain sequences over 333 PDB ids are clustered by
                single linkage, two PDB ids are linked when they share a chain
                cluster, and the resulting 127 homology groups replace the PDB
                id as the grouping variable.
  holdoutprot   ours: GroupKFold over SKEMPI 2.0's own `Hold_out_proteins`
                column, 159 groups. The curators' grouping, which PPIformer also
                split on, and the one measured at 0% leakage.

`hom99` and `holdoutprot` were built on unrelated principles and agree to within
0.013 on the reader-minus-head difference, which is why neither is quoted alone.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

SPLIT_DIR = Path(__file__).resolve().parent.parent / "data" / "split"
AVAILABLE = ("rdenet", "ppiformer", "hom99", "holdoutprot")


def load_split(name: str, pdbs: np.ndarray) -> np.ndarray:
    """Fold index per row; -1 where the PDB id is not in the split file."""
    path = SPLIT_DIR / f"{name}.csv"
    if not path.exists():
        raise FileNotFoundError(f"no split {name!r}; have {AVAILABLE}")
    m = pd.read_csv(path).drop_duplicates("PDB").set_index("PDB")["fold"].to_dict()
    return np.array([m.get(str(p), -1) for p in pdbs])
