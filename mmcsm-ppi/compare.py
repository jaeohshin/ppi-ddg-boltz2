#!/usr/bin/env python3
"""Score your predictions against mmCSM-PPI on its own benchmark.

    python compare.py --release data/release --predictions mine.csv --name "ours"

`mine.csv` has columns pdb, mutations, prediction. Rows are matched to the
benchmark by (PDB id, set of mutations), so order does not matter and any row
you do not predict is reported as a miss rather than silently dropped.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from mmcsm_bench.data import load_blind, load_predictions, load_zemu
from mmcsm_bench.score import bootstrap_complexes, metrics, paired_delta


def report(rows, mine, name, title):
    print(f"\n{title}  ({len(rows)} rows, "
          f"{len({r.pdb for r in rows})} complexes)")
    y = np.array([r.ddg for r in rows])
    g = np.array([r.pdb for r in rows])
    cols = {k: np.array([r.predictions.get(k, np.nan) for r in rows])
            for k in rows[0].predictions if not k.startswith("_")}
    if mine is not None:
        got = np.array([mine.get(r.key, np.nan) for r in rows])
        n_miss = int(np.isnan(got).sum())
        if n_miss:
            print(f"  ({n_miss} of {len(rows)} rows not present in your file)")
        cols[name] = got
    print(f"  {'method':20s} {'n':>5s} {'r':>8s} {'95% CI':>18s} {'rho':>7s} {'RMSE':>7s} {'MAE':>7s}")
    for k, p in cols.items():
        m = metrics(y, p)
        kk = np.isfinite(y) & np.isfinite(p)
        lo, hi = bootstrap_complexes(y[kk], p[kk], g[kk])
        print(f"  {k:20s} {m['n']:5d} {m['pearson']:+8.3f} "
              f"[{lo:+.3f},{hi:+.3f}] {m['spearman']:+7.3f} {m['rmse']:7.3f} {m['mae']:7.3f}")
    if mine is not None and "mmCSM-PPI" in cols:
        for stat in ("pearson", "spearman"):
            d = paired_delta(y, cols[name], cols["mmCSM-PPI"], g, stat=stat)
            verdict = "beats mmCSM-PPI" if d["beats"] else "not separated from mmCSM-PPI"
            print(f"  {name} - mmCSM-PPI ({stat}): {d['delta']:+.3f} "
                  f"[{d['lo']:+.3f}, {d['hi']:+.3f}] over {d['n_complexes']} complexes - {verdict}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--release", type=Path, default=Path("data/release"))
    ap.add_argument("--predictions", type=Path, default=None)
    ap.add_argument("--name", default="ours")
    a = ap.parse_args()

    mine = load_predictions(a.predictions, a.name) if a.predictions else None
    report(load_blind(a.release / "blindtest.csv"), mine, a.name,
           "BLIND TEST, four or more point mutations")
    report(load_zemu(a.release / "blindtest_zemu.csv"), mine, a.name,
           "ZEMu COMPARISON SET")


if __name__ == "__main__":
    main()
