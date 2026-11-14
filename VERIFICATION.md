# Verification

`reproduce.py` was run on the released PreFold-dG pickles and compared against the
reference solver in the working repository (`prefold-dg-mechanism`,
`analysis/reproduction_gap/`), which was written independently and checked itself against
a primal LSQR solve.

Full penalty curve, `rdenet`, 5,731 variants over 334 structures, 116 scored:

| alpha | per-structure PCC | reference |
|---|---|---|
| 0.001 | 0.5225 | |
| 0.00316 | 0.5470 | |
| 0.01 | 0.5708 | |
| 0.0316 | 0.5902 | |
| 0.1 | 0.6111 | 0.6111 |
| 0.316 | 0.6244 | 0.6244 |
| **1** | **0.6264** | **0.6264** |
| 3.16 | 0.6222 | 0.6222 |
| 10 | 0.6141 | 0.6141 |
| 31.6 | 0.6031 | 0.6031 |
| 100 | 0.5926 | 0.5926 |
| 316 | 0.5865 | 0.5865 |
| 1000 | 0.5819 | |

`holdoutprot`: best 0.5788 at alpha 10, against the 0.579 reported from the reference.

The cohort matches too: 5,817 mutant records pair into 5,817 variants, of which 86 carry
no fold in `rdenet`, leaving **5,731** — the published cohort size.

## A bug this check caught

The first implementation divided the shared kernel accumulator by each block's width
inside the per-block loop, so blocks added earlier were divided again by every later
block's width. It still ran, still produced a plausible U-shaped penalty curve, and
peaked at 0.600 — close enough to the target to have been accepted without this
comparison. Dividing each block's own Gram before adding it gives 0.6264.

`tests/test_reader.py::test_blocks_accumulate_independently` guards against its return.
