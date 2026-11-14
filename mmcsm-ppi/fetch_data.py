#!/usr/bin/env python3
"""Download the three mmCSM-PPI release files.

They are not vendored here: the site states no redistribution terms, so this
repository ships code and derived facts only.

    python fetch_data.py --out data/release
"""
import argparse
from pathlib import Path
import urllib.request

from mmcsm_bench.data import FILES

ap = argparse.ArgumentParser()
ap.add_argument("--out", type=Path, default=Path("data/release"))
a = ap.parse_args()
a.out.mkdir(parents=True, exist_ok=True)
for name, url in FILES.items():
    dest = a.out / f"{name}.csv"
    if dest.exists() and dest.stat().st_size:
        print(f"  {dest} already here"); continue
    urllib.request.urlretrieve(url, dest)
    print(f"  {dest}  {dest.stat().st_size} bytes")
