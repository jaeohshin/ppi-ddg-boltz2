#!/usr/bin/env python3
"""Recompute the README's Results table: reader, head, and the difference.

    python compare.py --pickles <dir> --split holdoutprot

Until this existed, the `head` column of that table was prose. The reader column
could be recomputed with `reproduce.py`, but the 58.4M-parameter head is not in this
repository and never will be, so the number it is compared against could not be
checked and the paired bootstrap behind each interval could not be rerun.

What makes it checkable is `head_predictions/`: the head's out-of-fold predictions
for the same cohort under the same three splits, one file per seed, as the retraining
runs wrote them. They are predictions, not scores, so everything downstream — the
per-structure correlation, the difference, the interval over complexes — is computed
here rather than quoted.

**The head column is a three-seed ensemble**: the seeds' predictions are averaged and
the average is scored once. Scoring each seed and averaging the scores is a different
and lower number, by 0.006 on `rdenet`, 0.002 on `hom99` and 0.008 on `holdoutprot`.
Both are printed, because the distinction is not visible in a single figure and the
published table reports the first.

The reader is fitted here exactly as `reproduce.py` fits it, and the difference is
reported **at every penalty in the grid**, not at one. That is deliberate and it is the
README's own rule: the inner selection surface is flat near its top, its argmax moves
between nominally identical runs, and a single reader score therefore carries about 0.02
of noise that the difference does not. A table of differences also checks a claim the
README makes and nothing else tested - that on `holdoutprot` the gap is positive with
its interval clearing zero at every grid point from alpha 1 to 1000.

It also explains two numbers that a single-alpha run would appear to contradict. The
README pairs `rdenet` conservatively, quoting the difference against the *lower* of two
reader runs (0.604, not 0.626), and its `hom99` reader figure of 0.573 is the best point
on that split's own curve rather than the value at any fixed alpha.

About 25 GB of RAM and some tens of minutes, the same as `reproduce.py`.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import time
from pathlib import Path

import numpy as np

from prefold_probe.cohort import build_cohort
from prefold_probe.evaluate import paired_bootstrap, per_structure_pcc, per_structure_srcc
from prefold_probe.features import difference_blocks
from prefold_probe.reader import BlockRidge, kernel
from prefold_probe.splits import load_split

HEAD_DIR = Path(__file__).resolve().parent / "head_predictions"
SEEDS = (42, 1337, 2024)
T0 = time.time()


def log(msg: str) -> None:
    print(f"[{time.time() - T0:6.0f}s] {msg}", flush=True)


def _verify(path: Path, body: bytes) -> None:
    sums = HEAD_DIR / "SHA256SUMS"
    if not sums.exists():
        return
    want = {n: d for d, n in
            (ln.split() for ln in sums.read_text().splitlines() if ln.strip())}
    expect = want.get(path.name.removesuffix(".gz"))
    if expect is None:
        return
    got = hashlib.sha256(body).hexdigest()
    if got != expect:
        raise SystemExit(
            f"{path} does not match {sums}.\n  expected {expect}\n  found    {got}\n"
            f"  These predictions are what the published head column was scored from.")


def load_head(split: str) -> dict[tuple[str, str], tuple[float, list[float]]]:
    """{(pdb, mut_id): (true_ddg, [pred per seed])} for one split."""
    per_seed = []
    for seed in SEEDS:
        f = HEAD_DIR / f"head_{split}_seed{seed}.csv.gz"
        if not f.exists():
            raise SystemExit(
                f"{f} is missing. `head_predictions/` holds the head's out-of-fold "
                f"predictions for rdenet, hom99 and holdoutprot; this split has none, "
                f"so the head column cannot be recomputed for it.")
        body = gzip.decompress(f.read_bytes())
        _verify(f, body)
        rows = {}
        for r in csv.DictReader(body.decode().splitlines()):
            rows[(r["pdb"], r["mut_id"])] = (float(r["true_ddg"]), float(r["pred_ddg"]))
        per_seed.append(rows)
    keys = set.intersection(*(set(d) for d in per_seed))
    return {k: (per_seed[0][k][0], [d[k][1] for d in per_seed]) for k in keys}


def fit_reader(pickles: Path, split: str):
    import joblib
    log("loading pickles")
    mut_records, mut_labels = joblib.load(pickles / "skempi_mut.pkl")
    wt_records, wt_labels = joblib.load(pickles / "skempi_wt.pkl")
    mut_ids, wt_ids, y, pdbs = build_cohort(mut_records, mut_labels, wt_records, wt_labels)
    log(f"{len(y)} variants over {len(set(pdbs))} structures")

    mut = {r["subpath"]: r for r in mut_records}
    wt = {r["subpath"]: r for r in wt_records}
    del mut_records, wt_records
    blocks = difference_blocks(mut, wt, mut_ids, wt_ids)
    del mut, wt

    folds = load_split(split, pdbs)
    keep = folds >= 0
    if (~keep).any():
        log(f"{(~keep).sum()} variants have no fold in {split}; dropped")
    blocks = {k: v[keep] for k, v in blocks.items()}
    y, pdbs, folds = y[keep], pdbs[keep], folds[keep]
    mut_ids = [m for m, k in zip(mut_ids, keep) if k]

    reader = BlockRidge()
    oof = {al: np.zeros_like(y) for al in reader.alphas}
    for f in sorted(set(folds)):
        te = folds == f
        tr = ~te
        log(f"fold {f}: {tr.sum()} train / {te.sum()} test")
        K = kernel(blocks, tr)
        fitted = BlockRidge().fit(K[np.ix_(tr, tr)], y[tr])
        for al, pr in fitted.predict_grid(K[np.ix_(te, tr)]).items():
            oof[al][te] = pr
        del K
    return oof, y, pdbs, mut_ids


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--pickles", required=True, type=Path,
                    help="directory holding skempi_mut.pkl and skempi_wt.pkl")
    ap.add_argument("--split", default="holdoutprot",
                    choices=("rdenet", "hom99", "holdoutprot"))
    ap.add_argument("--alpha", type=float, default=None,
                    help="report only this penalty instead of the whole grid")
    ap.add_argument("--out", type=Path, default=Path("comparison.json"))
    a = ap.parse_args()

    head = load_head(a.split)
    log(f"head predictions: {len(head)} variants, {len(SEEDS)} seeds")

    oof_by_alpha, y, pdbs, mut_ids = fit_reader(a.pickles, a.split)
    alphas = ([a.alpha] if a.alpha is not None
              else sorted(oof_by_alpha, key=float))

    # Score both on exactly the rows both cover, so no difference comes from the cohort.
    idx = [i for i, (p, m) in enumerate(zip(pdbs, mut_ids)) if (p, m) in head]
    if not idx:
        raise SystemExit(
            "no row is present in both the fitted cohort and the head predictions; the "
            "identifiers do not line up, so nothing here would be comparable.")
    idx = np.asarray(idx)
    dropped = len(y) - len(idx)
    if dropped:
        log(f"{dropped} fitted rows are absent from the head predictions; dropped")
    yv = y[idx]
    gv = np.asarray(pdbs)[idx]
    keys = [(pdbs[i], mut_ids[i]) for i in idx]

    head_label = np.array([head[k][0] for k in keys])
    if not np.allclose(head_label, yv, atol=1e-3):
        log(f"WARNING: labels differ from the head file by up to "
            f"{np.abs(head_label - yv).max():.4f}; the cohorts may not be the same rows")

    per_seed = np.array([[head[k][1][s] for k in keys] for s in range(len(SEEDS))])
    ens = per_seed.mean(axis=0)

    ens_pcc, n = per_structure_pcc(ens, yv, gv)
    ens_srcc, _ = per_structure_srcc(ens, yv, gv)
    seed_pccs = [per_structure_pcc(p, yv, gv)[0] for p in per_seed]

    print(f"\nsplit {a.split}   {len(idx)} variants, "
          f"{n} complexes scored (>= 10 variants each)\n")
    print(f"  head, 3-seed ensemble       PCC {ens_pcc:.4f}   SRCC {ens_srcc:.4f}"
          f"   <- the README's head column")
    print(f"  head, mean of seed scores   PCC {np.mean(seed_pccs):.4f}"
          f"   (sd {np.std(seed_pccs):.4f}, seeds "
          f"{', '.join(f'{v:.4f}' for v in seed_pccs)})")
    print(f"\n  {'alpha':>10s}  {'reader':>7s}  {'reader - head':>13s}  "
          f"{'95% interval':>20s}   clears 0")
    curve = []
    for al in alphas:
        pr = oof_by_alpha[al][idx]
        rp, _ = per_structure_pcc(pr, yv, gv)
        bt = paired_bootstrap(pr, ens, yv, gv)
        clears = "yes" if bt["lo"] > 0 else ("no" if bt["hi"] < 0 else "no (spans 0)")
        print(f"  {al:10.4g}  {rp:7.4f}  {bt['delta']:+13.4f}  "
              f"[{bt['lo']:+8.4f}, {bt['hi']:+8.4f}]   {clears}")
        curve.append({"alpha": float(al), "reader_per_structure_pcc": rp,
                      "delta": bt["delta"], "ci95": [bt["lo"], bt["hi"]],
                      "interval_excludes_zero": bool(bt["lo"] > 0)})
    best = max(curve, key=lambda r: r["reader_per_structure_pcc"])
    clear = [r["alpha"] for r in curve if r["interval_excludes_zero"]]
    n_clear = len(clear)
    print(f"\n  best reader on this grid: {best['reader_per_structure_pcc']:.4f} "
          f"at alpha {best['alpha']:g}, difference {best['delta']:+.4f} "
          f"[{best['ci95'][0]:+.4f}, {best['ci95'][1]:+.4f}]")
    if not clear:
        print(f"  the interval spans zero at every one of the {len(curve)} penalties")
    else:
        # The penalties that fail are the under-regularised end, so the useful
        # statement is the range over which the sign holds, not a bare count.
        tail = [r["alpha"] for r in curve
                if r["alpha"] >= min(clear)]
        contiguous = all(r["interval_excludes_zero"] for r in curve
                         if r["alpha"] >= min(clear))
        print(f"  the interval excludes zero at {n_clear} of {len(curve)} penalties"
              + (f", every one from alpha {min(clear):g} upward" if contiguous
                 else f" (not a contiguous range: {', '.join(f'{v:g}' for v in clear)})"))
    print(f"  the difference is what to quote; a single reader score carries about "
          f"0.02 of selection noise")

    report = {
        "split": a.split, "n_variants": int(len(idx)),
        "n_complexes_scored": int(n), "n_dropped_from_fit": int(dropped),
        "head_ensemble": {"per_structure_pcc": ens_pcc, "per_structure_srcc": ens_srcc},
        "head_per_seed_pcc": {str(s): float(v) for s, v in zip(SEEDS, seed_pccs)},
        "head_seed_mean_pcc": float(np.mean(seed_pccs)),
        "head_seed_sd_pcc": float(np.std(seed_pccs)),
        "penalty_curve": curve,
        "best_reader": best,
        "n_penalties_with_interval_excluding_zero": n_clear,
        "alphas_with_interval_excluding_zero": clear,
    }
    a.out.write_text(json.dumps(report, indent=1))
    print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
