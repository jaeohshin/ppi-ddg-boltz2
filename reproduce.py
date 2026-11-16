#!/usr/bin/env python3
"""Fit the reader on the released PreFold-dG pickles and score it.

    python reproduce.py --pickles <dir> --split rdenet

Reads `skempi_mut.pkl` and `skempi_wt.pkl`, builds the mutant-minus-wild-type
blocks, and runs the published three-fold cross-validation, reporting the
per-structure PCC for every penalty in the grid.  Nothing is selected here: the
curve is the output, because the inner selection surface is flat and its argmax
is not stable (see README, "Quote the gap, not a single score").

About 25 GB of RAM and some tens of minutes on the full cohort.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from prefold_probe.cohort import build_cohort
from prefold_probe.inputs import verify as verify_inputs
from prefold_probe.evaluate import per_structure_pcc, per_structure_srcc
from prefold_probe.features import difference_blocks
from prefold_probe.reader import BlockRidge, kernel
from prefold_probe.splits import load_split

T0 = time.time()


def log(msg: str) -> None:
    print(f"[{time.time() - T0:6.0f}s] {msg}", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pickles", required=True, type=Path,
                    help="directory holding skempi_mut.pkl and skempi_wt.pkl")
    ap.add_argument("--split", default="rdenet")
    ap.add_argument("--out", type=Path, default=Path("reader_scores.json"))
    ap.add_argument("--verify-digest", action="store_true",
                    help="hash the pickles against docs/PROVENANCE.md as well as "
                         "checking their size; about half a minute for 6.8 GB")
    args = ap.parse_args()

    log("checking the inputs against docs/PROVENANCE.md")
    verify_inputs(args.pickles, digest=args.verify_digest, log=log)

    import joblib
    log("loading pickles")
    mut_records, mut_labels = joblib.load(args.pickles / "skempi_mut.pkl")
    wt_records, wt_labels = joblib.load(args.pickles / "skempi_wt.pkl")
    mut_ids, wt_ids, y, pdbs = build_cohort(mut_records, mut_labels, wt_records, wt_labels)
    log(f"{len(y)} variants over {len(set(pdbs))} structures")

    mut = {r["subpath"]: r for r in mut_records}
    wt = {r["subpath"]: r for r in wt_records}
    del mut_records, wt_records
    log("building difference blocks")
    blocks = difference_blocks(mut, wt, mut_ids, wt_ids)
    del mut, wt
    log(f"blocks: " + ", ".join(f"{k}{v.shape}" for k, v in blocks.items()))

    folds = load_split(args.split, pdbs)
    if (folds < 0).any():
        log(f"{(folds < 0).sum()} variants have no fold in {args.split}; dropped")
    keep = folds >= 0
    blocks = {k: v[keep] for k, v in blocks.items()}
    y, pdbs, folds = y[keep], pdbs[keep], folds[keep]

    reader = BlockRidge()
    oof = {a: np.zeros_like(y) for a in reader.alphas}
    for f in sorted(set(folds)):
        te = folds == f
        tr = ~te
        log(f"fold {f}: {tr.sum()} train / {te.sum()} test — kernel")
        K = kernel(blocks, tr)
        fitted = BlockRidge().fit(K[np.ix_(tr, tr)], y[tr])
        preds = fitted.predict_grid(K[np.ix_(te, tr)])
        for a, p in preds.items():
            oof[a][te] = p
        del K
        log(f"fold {f} done")

    rows = []
    for a in reader.alphas:
        pcc, n = per_structure_pcc(oof[a], y, pdbs)
        srcc, _ = per_structure_srcc(oof[a], y, pdbs)
        rows.append({"alpha": float(a), "per_structure_pcc": pcc,
                     "per_structure_srcc": srcc, "n_scored": n})
        print(f"  alpha {a:>10.4g}   PCC {pcc:.4f}   SRCC {srcc:.4f}   ({n} structures)")
    best = max(rows, key=lambda r: r["per_structure_pcc"])
    print(f"\nbest on this grid: alpha {best['alpha']:g}, "
          f"PCC {best['per_structure_pcc']:.4f}, SRCC {best['per_structure_srcc']:.4f}")
    args.out.write_text(json.dumps(
        {"split": args.split, "n_variants": int(len(y)),
         "n_structures": int(len(set(pdbs))), "curve": rows}, indent=1))
    log(f"wrote {args.out}")


if __name__ == "__main__":
    main()
