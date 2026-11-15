# Inside the 128 numbers

Once `block_ablation.py` showed that only the 128-wide pair block generalises
(RESULTS.md §3.1b), the feature set became small enough to examine directly —
which 295,040 dimensions never were. This is what is in it.

All numbers are per-structure PCC over the 116 complexes with at least ten
variants, three-fold, on the two homology-controlled splits. `pair_block.py`
caches the block; everything here runs in a couple of minutes from that cache.

## The signal is about four dimensions wide

Dimensions ranked by their own within-complex PCC **recomputed inside each
training fold**, then the top k refitted. In-fold ranking matters: ranking on the
full cohort and refitting on it inflated k=1 from 0.508 to 0.524 by selection
over 128 candidates, which is how the first pass of this got it wrong.

| k | `hom99` | `holdoutprot` |
|---|---|---|
| 1 | 0.498 | 0.508 |
| 2 | 0.542 | 0.528 |
| **4** | **0.567** | **0.552** |
| 8 | 0.564 | 0.555 |
| 16 | 0.554 | 0.549 |
| 64 | 0.570 | 0.565 |
| 128 | 0.578 | 0.576 |

Paired over complexes, 10,000 draws:

| | `hom99` | `holdoutprot` |
|---|---|---|
| k=4 − k=128 | −0.012 [−0.034, +0.014] | −0.024 [−0.062, +0.004] |
| k=1 − k=128 | **−0.080 [−0.113, −0.050]** | **−0.068 [−0.106, −0.032]** |

**Four dimensions are indistinguishable from all 128; one is significantly
worse.** For scale, the 58.4M-parameter head scores 0.505 and 0.524 on these two
splits, so four of these numbers, read linearly, are above it and one is roughly
level with it on `hom99` and below it on `holdoutprot`.

k=16 scores *below* k=4 and k=8 on both splits. The ranking is a poor selector in
the middle of its range, and the heavy-ridge fit on all 128 is still the best
arm. So the practical recommendation is unchanged — **fit all 128 with a heavy
penalty** — and this section is a statement about the intrinsic dimension of the
signal, not a recommendation to select features.

## The strong dimensions are copies of each other

The eight strongest single dimensions on `holdoutprot`, standardised, correlate
with each other up to |r| = 0.91:

|  | 97 | 78 | 99 | 105 | 42 | 46 | 18 | 107 |
|---|---|---|---|---|---|---|---|---|
| 97 | +1.00 | +0.69 | +0.59 | −0.91 | −0.81 | +0.72 | −0.88 | +0.60 |
| 105 | −0.91 | −0.69 | −0.56 | +1.00 | +0.87 | −0.80 | +0.91 | −0.64 |
| 42 | −0.81 | −0.74 | −0.66 | +0.87 | +1.00 | −0.75 | +0.82 | −0.47 |
| 18 | −0.88 | −0.62 | −0.43 | +0.91 | +0.82 | −0.83 | +1.00 | −0.59 |

**Participation ratio of those eight: 1.90.** They are not eight directions, they
are about two, written down repeatedly with sign flips. The whole 128-column
block has a participation ratio of 7.9, with three principal components covering
50% of the column variance and 22 covering 90%.

So there are three different "dimensions" in play and they should not be
confused: 128 columns, about 8 directions of variance, and a predictive subspace
of about 2 that four columns are enough to span.

## No dimension should be named

The top-k sets differ across the three folds — two or three distinct sets at
every k below 128. **The identity of the strongest dimension is not stable; only
the subspace is.** Any claim of the form "dimension 97 encodes binding" is not
supported by this, and the earlier framing of it has been removed rather than
softened.

## What this does not answer

The subspace has not been connected to anything physical. The obvious tests —
whether its projection tracks interface size, burial, the ΔΔG sign, the mutated
residue's class, or hydrophobicity — are all cheap now and none of them is done.
Until one is, "about two dimensions" is a statement about the geometry of the
representation and not about biophysics.

Nor does any of it bear on the sequence blocks beyond RESULTS.md §3.1b: they were
shown not to generalise, not to be empty.

Artefacts: `nested128.json`, `pair_block_scores.json`, `pairblock_bootstrap.json`.
