#!/usr/bin/env python3
"""Why 118 of the 1,935 benchmark rows had no feature row, and which are recoverable.

The ridge run scored 586 of the 595 blind mutations and trained on 1,231 of the
1,340 released training rows. This script classifies every dropped row and writes
the list of structures that would have to be folded to close the gap.

    python diagnose_missing.py --experiments /data/work/pythia_ppi/experiments \
                               --mmcsm /data/work/deltadeltaG/mmCSM/data \
                               --skempi /data/work/PreFold/data/skempi_v2.csv

Paths are local to the machine the ridge was fitted on; nothing here is shipped.
"""
from __future__ import annotations

import argparse, collections, csv, json, math
from pathlib import Path

RT = 8.314 / 4184 * (273.15 + 25)          # kcal/mol at 25 C


def read(p, d=","):
    return list(csv.DictReader(open(p), delimiter=d))


def key(pdb, muts):
    return pdb.upper(), ",".join(sorted(muts.upper().split(",")))


def flip(muts):
    return ",".join(t[-1] + t[1:-1] + t[0] for t in muts.split(","))


def ddg_mmcsm(s):
    """mmCSM's sign convention is ln(K_wt) - ln(K_mt), the negative of SKEMPI's."""
    try:
        return -RT * math.log(float(s["Affinity_mut_parsed"]) / float(s["Affinity_wt_parsed"]))
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--experiments", type=Path, required=True)
    ap.add_argument("--mmcsm", type=Path, required=True)
    ap.add_argument("--skempi", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=Path("data/to_fold.csv"))
    a = ap.parse_args()

    raw = read(a.skempi, ";")
    mapped = read(a.experiments / "mini_full_compare_20261103/mapped_rows.csv")
    ix, src = collections.defaultdict(dict), collections.defaultdict(list)
    for r in mapped:
        for j in json.loads(r["source_rows"]):
            ix[key(r["pdb"], raw[j - 1]["Mutation(s)_PDB"])][r["mut_id"]] = r
            src[r["mut_id"]].append(raw[j - 1])
    folded = {r["pdb"].upper() for r in mapped}

    train = list(csv.reader(open(a.mmcsm / "train_set.csv")))[1:]
    miss = read(a.experiments / "mmcsm_ridge_blind_20261106/missing.csv")

    buckets, to_fold = collections.Counter(), []
    detail = collections.defaultdict(collections.Counter)
    for m in miss:
        pdb, direction = m["pdb"].upper(), m["direction"]
        lookup = flip(m["mutation"]) if direction == "reverse" else m["mutation"]
        cands = ix[key(pdb, lookup)]

        if pdb not in folded:
            b = "complex never folded"
        elif not cands:
            b = "complex folded, this mutant was not"
        elif len(cands) > 1 and m["set"] == "train":
            want = float(train[int(m["source_row"]) - 1][-2])
            sign = -1 if direction == "reverse" else 1
            hits = {mid for mid in cands for s in src[mid]
                    if (d := ddg_mmcsm(s)) is not None and abs(sign * d - want) < 0.01}
            b = ("several interfaces, resolved by label" if len(hits) == 1
                 else "several interfaces, unresolved")
        else:
            b = "other"
        buckets[b] += 1
        detail[b][f"{pdb} ({m['set']})"] += 1
        if b in ("complex never folded", "complex folded, this mutant was not"):
            to_fold.append({"pdb": pdb, "mutations": m["mutation"],
                            "set": m["set"], "reason": b})

    print("WHY 118 ROWS WERE DROPPED\n")
    for b, n in buckets.most_common():
        print(f"  {n:3d}  {b}")
        for k, v in detail[b].most_common():
            print(f"         {v:3d}  {k}")
        print()

    per = collections.defaultdict(set)
    for r in to_fold:
        per[r["pdb"]].add(r["mutations"])
    n_struct = sum(len(v) + (1 if p not in folded else 0) for p, v in per.items())
    print("STRUCTURES TO FOLD TO CLOSE THE GAP")
    for p, v in sorted(per.items()):
        wt = " + WT" if p not in folded else " (WT already folded)"
        print(f"  {p}: {len(v)} mutants{wt}")
    print(f"  total: {n_struct} structures")

    a.out.parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["pdb", "mutations", "set", "reason"])
        w.writeheader(); w.writerows(sorted(to_fold, key=lambda r: (r["pdb"], r["mutations"])))
    print(f"\n  wrote {a.out}")


if __name__ == "__main__":
    main()
