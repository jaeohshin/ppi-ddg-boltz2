#!/usr/bin/env python3
"""Extract the 128-wide pair block once, then score it on every split.

    python pair_block.py --pickles <dir>          # extracts and scores
    python pair_block.py --cache pair_block.npz   # scores from the cache

`block_ablation.py` showed that on the homology-controlled splits the two
147,456-wide sequence blocks contribute nothing measurable: the pair block alone
is within 0.005 per-structure PCC of all three blocks together, while on the
56%-leaked `rdenet` split they are worth +0.038.  That makes the 128-wide block
the model to report for anything that claims generalisation, and it is small
enough that the whole cohort fits in 3 MB.

Two things the wide ablation could not do are done here.  The penalty grid is
extended far upward -- 128 features want a much heavier ridge than 295,040 do,
and the shared grid ran out at alpha 1000 with the curve still rising -- and the
block is cached, so any further question about these 128 numbers is instant.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from prefold_probe.evaluate import per_structure_pcc, per_structure_srcc
from prefold_probe.reader import BlockRidge, _block_gram
from prefold_probe.splits import AVAILABLE, load_split

T0 = time.time()
ALPHAS = np.logspace(-2, 7, 37)  # 128 features need the top end; see the docstring


def log(m: str) -> None:
    print(f"[{time.time() - T0:6.0f}s] {m}", flush=True)


def extract(pickles: Path, cache: Path) -> None:
    import joblib
    from prefold_probe.cohort import build_cohort
    log("loading pickles")
    mut_records, mut_labels = joblib.load(pickles / "skempi_mut.pkl")
    wt_records, wt_labels = joblib.load(pickles / "skempi_wt.pkl")
    mut_ids, wt_ids, y, pdbs = build_cohort(mut_records, mut_labels, wt_records, wt_labels)
    mut = {r["subpath"]: r for r in mut_records}
    wt = {r["subpath"]: r for r in wt_records}
    del mut_records, wt_records
    wz = np.stack([np.asarray(mut[s]["wz"], dtype=np.float32).ravel()
                   - np.asarray(wt[w]["wz"], dtype=np.float32).ravel()
                   for s, w in zip(mut_ids, wt_ids)])
    np.savez_compressed(cache, wz=wz, y=y, pdbs=np.asarray(pdbs, dtype=object),
                        mut_ids=np.asarray(mut_ids, dtype=object))
    log(f"wrote {cache}: wz{wz.shape}, {cache.stat().st_size/2**20:.1f} MiB")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pickles", type=Path)
    ap.add_argument("--cache", type=Path, default=Path("pair_block.npz"))
    ap.add_argument("--splits", nargs="+", default=list(AVAILABLE))
    ap.add_argument("--out", type=Path, default=Path("pair_block_scores.json"))
    args = ap.parse_args()

    if not args.cache.exists():
        if args.pickles is None:
            ap.error(f"{args.cache} missing and no --pickles given")
        extract(args.pickles, args.cache)

    d = np.load(args.cache, allow_pickle=True)
    wz, y_all, pdbs_all = d["wz"], d["y"], d["pdbs"]
    log(f"{wz.shape[0]} variants, {wz.shape[1]} features, "
        f"{len(set(pdbs_all.tolist()))} structures")

    result = {}
    for split in args.splits:
        folds_all = load_split(split, pdbs_all)
        keep = folds_all >= 0
        X, y, pdbs, folds = wz[keep], y_all[keep], pdbs_all[keep], folds_all[keep]
        oof = {a: np.zeros_like(y) for a in ALPHAS}
        for f in sorted(set(folds)):
            te = folds == f
            tr = ~te
            K = _block_gram(X, tr)
            fitted = BlockRidge(alphas=tuple(ALPHAS)).fit(K[np.ix_(tr, tr)], y[tr])
            for a, p in fitted.predict_grid(K[np.ix_(te, tr)]).items():
                oof[a][te] = p
        rows = []
        for a in ALPHAS:
            pcc, n = per_structure_pcc(oof[a], y, pdbs)
            srcc, _ = per_structure_srcc(oof[a], y, pdbs)
            rows.append({"alpha": float(a), "per_structure_pcc": pcc,
                         "per_structure_srcc": srcc, "n_scored": n})
        best = max(rows, key=lambda r: r["per_structure_pcc"])
        interior = rows[0]["per_structure_pcc"] < best["per_structure_pcc"] > rows[-1]["per_structure_pcc"]
        result[split] = {"n_variants": int(len(y)), "n_structures": int(len(set(pdbs.tolist()))),
                         "curve": rows, "best": best, "argmax_interior": bool(interior)}
        log(f"{split:12s} best PCC {best['per_structure_pcc']:.4f} "
            f"SRCC {best['per_structure_srcc']:.4f} at alpha {best['alpha']:.4g}"
            f"{'' if interior else '   <-- STILL AT A GRID EDGE'}")
    args.out.write_text(json.dumps(result, indent=1))
    log(f"wrote {args.out}")


if __name__ == "__main__":
    main()
