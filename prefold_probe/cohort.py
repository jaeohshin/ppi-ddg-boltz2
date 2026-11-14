"""Assembling the cohort the published loader assembles.

The upstream loader groups a record by `subpath.split('_')[2]`, the PDB id, and
pairs every mutant of a PDB id with the first wild-type record carrying that id.
This package reproduces that pairing exactly, including its one known defect:
2C5D, 3SE3 and 3SE4 appear with more than one wild-type entry, so 35 variants
are paired against the wrong one.  The defect is kept rather than fixed, because
correcting it would make these numbers incomparable with the published ones.
"""
from __future__ import annotations

import numpy as np

RT_LOG10 = 1.364  # kcal/mol per unit of the log10(Kd, nM) label, at 25 C


def pdb_of(subpath: str) -> str:
    return subpath.split("_")[2]


def build_cohort(mut_records: list, mut_labels, wt_records: list, wt_labels):
    """Return (mut_ids, wt_ids, y, pdbs) in the loader's own order.

    `y` is the mutant-minus-wild-type difference of log10(Kd in nM), which is the
    label the published model is trained and scored on.  It is not kcal/mol;
    multiply by 1.364 for that.
    """
    wt_first: dict[str, int] = {}
    for j, rec in enumerate(wt_records):
        wt_first.setdefault(pdb_of(rec["subpath"]), j)

    mut_ids, wt_ids, y, pdbs = [], [], [], []
    for i, rec in enumerate(mut_records):
        g = pdb_of(rec["subpath"])
        j = wt_first.get(g)
        if j is None:
            continue
        mut_ids.append(rec["subpath"])
        wt_ids.append(wt_records[j]["subpath"])
        y.append(float(mut_labels[i]) - float(wt_labels[j]))
        pdbs.append(g)
    return mut_ids, wt_ids, np.asarray(y, dtype=np.float64), np.asarray(pdbs)
