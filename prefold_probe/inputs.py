"""The pinned inputs, and a check that what is on disk is still them.

`docs/PROVENANCE.md` records a sha256 and a byte count for each released pickle. A
digest nobody verifies is half a protection: the numbers in this repository are only
about these files, and a silently different release would change every one of them
while everything still ran.

Hashing 6.8 GB takes about half a minute, so `verify()` checks the byte count always
and the digest only when asked. The sizes alone catch a different release; the digest
is there for when a result is being quoted.
"""
from __future__ import annotations

import hashlib
from pathlib import Path

# From docs/PROVENANCE.md, "The pinned inputs". These are the authors' Zenodo release
# at zenodo.org/records/20236663, whose own published byte counts agree with these.
PINNED = {
    "skempi_mut.pkl": (
        "f89501f2ae2bca19e3725c8e1fc050002f02ac7094f7a4b1b5fe96de3ef90e91", 6_443_150_681),
    "skempi_wt.pkl": (
        "a106fb130d91ffbca2826e685f7cf7356680a9fb2fe81d85c208c883d88bb0b2", 375_498_913),
}


def sha256(path: Path, chunk: int = 1 << 22) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def verify(pickles: Path, digest: bool = False, log=print) -> None:
    """Check the pickles in `pickles` against PINNED. Raises on a mismatch."""
    for name, (want_sha, want_size) in PINNED.items():
        f = Path(pickles) / name
        if not f.exists():
            raise SystemExit(
                f"{f} not found.\n"
                f"  The two released pickles are at zenodo.org/records/20236663 "
                f"({want_size / 1e9:.2f} GB for {name}).")
        size = f.stat().st_size
        if size != want_size:
            raise SystemExit(
                f"{f} is {size:,} bytes; the pinned release is {want_size:,}.\n"
                f"  Every number in this repository is about the pinned file, so this "
                f"is a different input rather than a tolerable difference. "
                f"docs/PROVENANCE.md records both.")
        if digest:
            got = sha256(f)
            if got != want_sha:
                raise SystemExit(
                    f"{f} has sha256 {got},\n  but the pinned release is {want_sha}.")
            log(f"  {name}: {size:,} bytes, sha256 verified")
        else:
            log(f"  {name}: {size:,} bytes as pinned (pass --verify-digest to hash it)")
