#!/usr/bin/env python3
"""Which of the three blocks actually generalises?

    python block_ablation.py --pickles <dir> --splits rdenet hom99 holdoutprot

The leave-one-family-out run on the mmCSM cohort found that dropping both
147,456-wide sequence blocks and keeping only the 128-wide pair block cost
0.0011 Pearson once the protein family was new (docs/FAMILY_SPLIT.md).  That was
one cohort and one kind of split.  This repeats the test here, on the published
SKEMPI cohort and all three of this repository's splits, with the per-structure
metric the rest of the results use.

Per fold the three per-block Gram matrices are computed once -- standardised on
that fold's training rows, each divided by its own width -- and every arm is a
sum of a subset of them, so the whole ablation costs one pass over the features
rather than one per arm.  Nothing is selected: the output is the penalty curve
for each arm, as elsewhere in this repository.

About 25 GB of RAM; roughly five minutes per split after the pickles are loaded.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from prefold_probe.cohort import build_cohort
from prefold_probe.evaluate import per_structure_pcc, per_structure_srcc
from prefold_probe.features import difference_blocks
from prefold_probe.reader import BlockRidge, _block_gram
from prefold_probe.splits import load_split

T0 = time.time()

# Each arm names the blocks it keeps.  "all" is the published reader.
ARMS = {
    "all": ("wsi", "ws", "wz"),
    "wz": ("wz",),
    "seq": ("wsi", "ws"),
    "wsi": ("wsi",),
    "ws": ("ws",),
}


def log(msg: str) -> None:
    print(f"[{time.time() - T0:6.0f}s] {msg}", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pickles", required=True, type=Path)
    ap.add_argument("--splits", nargs="+", default=["rdenet", "hom99", "holdoutprot"])
    ap.add_argument("--out", type=Path, default=Path("block_ablation.json"))
    ap.add_argument("--arms", nargs="+", default=list(ARMS),
                    help="subset of arms to run; the Gram cost is shared anyway")
    ap.add_argument("--save-oof", type=Path,
                    help="write the out-of-fold predictions, for paired bootstraps")
    args = ap.parse_args()

    import joblib
    log("loading pickles")
    mut_records, mut_labels = joblib.load(args.pickles / "skempi_mut.pkl")
    wt_records, wt_labels = joblib.load(args.pickles / "skempi_wt.pkl")
    mut_ids, wt_ids, y_all, pdbs_all = build_cohort(mut_records, mut_labels,
                                                    wt_records, wt_labels)
    mut = {r["subpath"]: r for r in mut_records}
    wt = {r["subpath"]: r for r in wt_records}
    del mut_records, wt_records
    log("building difference blocks")
    blocks_all = difference_blocks(mut, wt, mut_ids, wt_ids)
    del mut, wt
    log(f"{len(y_all)} variants over {len(set(pdbs_all))} structures")

    alphas = BlockRidge().alphas
    arms = {k: v for k, v in ARMS.items() if k in args.arms}
    if not arms:
        ap.error(f"no arms selected; have {list(ARMS)}")
    result = {}
    saved = {}
    for split in args.splits:
        folds_all = load_split(split, pdbs_all)
        keep = folds_all >= 0
        blocks = {k: v[keep] for k, v in blocks_all.items()}
        y, pdbs, folds = y_all[keep], pdbs_all[keep], folds_all[keep]
        log(f"{split}: {keep.sum()} variants, {len(set(pdbs))} structures, "
            f"{(~keep).sum()} without a fold")

        oof = {arm: {a: np.zeros_like(y) for a in alphas} for arm in arms}
        for f in sorted(set(folds)):
            te = folds == f
            tr = ~te
            log(f"{split} fold {f}: {tr.sum()} train / {te.sum()} test")
            grams = {}
            for name, X in blocks.items():
                grams[name] = _block_gram(X, tr)
                log(f"  {name} gram done")
            for arm, names in arms.items():
                K = sum(grams[n] for n in names)
                fitted = BlockRidge().fit(K[np.ix_(tr, tr)], y[tr])
                for a, p in fitted.predict_grid(K[np.ix_(te, tr)]).items():
                    oof[arm][a][te] = p
                del K
            del grams

        if args.save_oof is not None:
            saved[split] = {"y": y, "pdbs": np.asarray(pdbs, dtype=object),
                            "alphas": np.asarray(alphas),
                            **{f"oof_{arm}": np.stack([oof[arm][a] for a in alphas])
                               for arm in arms}}
        curves = {}
        for arm in arms:
            rows = []
            for a in alphas:
                pcc, n = per_structure_pcc(oof[arm][a], y, pdbs)
                srcc, _ = per_structure_srcc(oof[arm][a], y, pdbs)
                rows.append({"alpha": float(a), "per_structure_pcc": pcc,
                             "per_structure_srcc": srcc, "n_scored": n})
            curves[arm] = rows
            best = max(rows, key=lambda r: r["per_structure_pcc"])
            log(f"{split:12s} {arm:4s} ({len(arms[arm])} blocks)  "
                f"best PCC {best['per_structure_pcc']:.4f} at alpha {best['alpha']:g}  "
                f"SRCC {best['per_structure_srcc']:.4f}")
        result[split] = {"n_variants": int(len(y)),
                         "n_structures": int(len(set(pdbs))),
                         "arms": {k: {"blocks": list(arms[k]), "curve": v}
                                  for k, v in curves.items()}}
        args.out.write_text(json.dumps(result, indent=1))
        if args.save_oof is not None:
            np.savez_compressed(args.save_oof, **{f"{s}__{k}": v
                                                  for s, d in saved.items()
                                                  for k, v in d.items()})

    print("\nbest per-structure PCC on each split's penalty curve")
    print("arm   width      " + "  ".join(f"{s:>12s}" for s in args.splits))
    widths = {"all": 295040, "wz": 128, "seq": 294912, "wsi": 147456, "ws": 147456}
    for arm in arms:
        cells = []
        for s in args.splits:
            b = max(result[s]["arms"][arm]["curve"], key=lambda r: r["per_structure_pcc"])
            cells.append(f"{b['per_structure_pcc']:12.4f}")
        print(f"{arm:5s} {widths[arm]:>7d}    " + "  ".join(cells))
    log(f"wrote {args.out}")


if __name__ == "__main__":
    main()
