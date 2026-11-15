"""Recovering which complex a released mmCSM-PPI row refers to.

mmCSM-PPI's released files name a complex by its four-letter PDB id alone. A
structure deposited with two distinct interfaces -- 3SE3_B_A and 3SE3_B_C --
therefore collapses to one key that carries two different measured affinities,
and nothing in the row says which one it is.

The row's own experimental value settles it. That value is a measurement, not a
prediction, and it reproduces the SKEMPI affinities behind each candidate
exactly. An *exact* match is required and exactly one candidate must match, so
this recovers a row's identity rather than choosing the closest fit; the
distinction matters because the discriminant is the label.
"""
import math

RT = 8.314 / 4184 * (273.15 + 25)  # kcal/mol at 25 C


def ddg(kd_wt: float, kd_mut: float) -> float:
    """mmCSM-PPI's convention: ln(K_wt) - ln(K_mt).

    That is the negative of SKEMPI's RT*ln(Kd_mut/Kd_wt), which is why a naive
    comparison against SKEMPI misses every row.
    """
    return -RT * math.log(kd_mut / kd_wt)


def resolve(candidates, y, direction, tol=1e-4):
    """The single candidate whose measured ddG is exactly this row's label.

    `candidates` maps a candidate id to the ddG values of the SKEMPI rows behind
    it -- a candidate may stand for more than one row. `direction` is 'forward'
    or 'reverse' as the released file gives it; a reverse row carries the
    negated label. Returns None when no candidate matches or more than one does,
    which leaves the row dropped rather than guessed.
    """
    if direction not in ("forward", "reverse"):
        raise ValueError(f"direction must be forward or reverse, not {direction!r}")
    want = y if direction == "forward" else -y
    hits = [k for k, ds in candidates.items() if any(abs(d - want) < tol for d in ds)]
    return hits[0] if len(hits) == 1 else None
