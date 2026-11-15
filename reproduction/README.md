# Does the released result survive replacing the authors' features?

Every fit in this repository reads the PreFold-dG authors' Zenodo pickles. That
is deliberate — it makes the comparison against their published head a comparison
of readers rather than of pipelines — but it leaves the released numbers resting
on files produced elsewhere, by settings the release does not record.

This directory answers one question: **swap the whole feature producer for a
local Boltz extraction and refit. Does the score hold?**

## Answer

It holds, with the author-matched extraction. Same 440 of 586 blind rows, 59 of
66 complexes, both arms on identical rows and the authors' labels:

| features | Pearson | rho | vs authors' features | vs mmCSM-PPI |
|---|---|---|---|---|
| authors' Zenodo pickles | 0.8723 | 0.8185 | — | +0.117 [+0.026, +0.342] |
| local `scaled_v1_r10_tplfull` | **0.8751** | 0.8165 | **+0.003 [−0.032, +0.022]** | **+0.119 [+0.022, +0.344]** |
| local `scaled_v1_recycle3` | 0.8187 | 0.7347 | −0.054 [−0.131, −0.029] | +0.063 [−0.049, +0.279] |
| mmCSM-PPI (unchanged) | 0.7556 | 0.7190 | | |

The two local rows are the whole point: **local extractions are not
interchangeable.** `scaled_v1_r10_tplfull` is recycling 10 (the paper's setting)
with every wild-type complex's deposited structure as a template;
`scaled_v1_recycle3` is Boltz's default recycling with no template. Both were
run. Only the first reproduces the result.

`scaled_v1_recycle3` was the variant this test was first pre-registered on — see
[PLAN.md](PLAN.md) — on a misreading of `docs/PROVENANCE.md`, which recorded
recycling 3 as the project's configuration without saying it is unmatched. Both
numbers are reported; the recycle3 row is a wrong-variant result, not a property
of local folding.

## Why the extraction matters more than its raw fidelity suggests

The reader never sees Boltz output. It sees `mut − wt`, a small residual of two
nearly identical vectors, and subtraction amplifies any disagreement — hardest in
the block where the residual is the smallest fraction of the vector:

| block | ‖Δ‖ / ‖wt‖ | raw-vector cosine | difference cosine | angle amplification |
|---|---|---|---|---|
| `wsi` | 0.115 | 0.9867 | 0.9392 | 2.1× |
| `ws`  | 0.105 | 0.9925 | 0.8567 | 4.4× |
| `wz`  | **0.048** | 0.9929 | **0.6343** | **7.2×** |

(`scaled_v1_recycle3`, 400 blind-cohort rows.) `wz` — the only block that
generalises across protein families — carries a residual that is 4.8% of the
vector norm. Raw output matched to cosine 0.993 still leaves the quantity the
model reads disagreeing at 0.63. The matched extraction reaches raw 0.9998 and
difference 0.9764 on `wz`, and that is the whole gap between the two local rows.

**So raw-output fidelity is the wrong reproduction target.** Repeated attempts to
close the gap by tightening agreement with the authors' raw tensors did not
converge, and this is why.

What is not the cause, on evidence recorded in `prefold-dg-mechanism`: the
checkpoint (the author confirmed it is bit-identical), the MSA scheme (ColabFold
greedy default, confirmed), and run-to-run nondeterminism (one pipeline rerun
across 120 GPU tasks self-agrees at tensor cosine 1.00000–1.00005, though only
0.994–0.998 on the difference). What remains is the authors' undocumented
sequence repair — aligning the RCSB and SKEMPI sequences and filling or trimming
— volunteered in correspondence, absent from the paper and their notebook, never
replicated.

## Joining the two producers

This had to be solved before anything could be compared. The two sides key their
records differently:

| | positions | token order |
|---|---|---|
| authors' prepared cohort | SKEMPI `Mutation(s)_cleaned` | sorted |
| local extraction | SKEMPI `Mutation(s)_PDB` | SKEMPI's original |

SKEMPI prints both numbering columns in the same token order for every row, so
zipping them yields a per-(complex, chain) position dictionary — 4,008 entries,
no conflicts — which normalises local keys into the authors' convention and makes
the join a dictionary lookup. `build_join.py` joins 5,326 of 5,729 mutant records
and 335 of 336 wild-type records. Eight antibody complexes differ only in how the
chain groups are spelled (`1DAN_HL_TU` vs `1DAN_HL_UT`) and are matched by trying
the permutations.

The join is validated against the labels stored inside both pickles, independently
of the features: Pearson 0.988, median difference 0.005, the residual explained by
SKEMPI duplicate-affinity aggregation. Labels are never taken from the local side;
both arms use the cohort's, so label-convention differences cannot enter the
comparison.

## Running it

```
PY=/path/to/python ./run_variant.sh scaled_v1_r10_tplfull
```

Inputs are machine-specific and located by [`paths.py`](paths.py); override them
with `PREFOLD_EXPERIMENTS`, `PREFOLD_AUTHOR_PKLS` and `PREFOLD_LOCAL_PKLS`.
Nothing in this directory ships feature data — the pickles are the authors'
release and a local Boltz run, both several gigabytes.

| script | |
|---|---|
| `build_join.py` | the SKEMPI-mediated join; writes `join_map_<variant>.json` |
| `compare_sources.py` | both arms, identical rows, the blind fit verbatim from the locked experiment |
| `diagnose.py` | validates the join on labels; splits the disagreement into wild-type, mutant and difference, per block |
| `block_arms.py` | each block alone, both producers |
| `weight_swap.py` | the whole (block weight, alpha) surface for both producers, selecting nothing |

Results are in [`results/`](results/), one directory per variant.

## What this does not license

This is a producer-swap test: one producer for a whole fit. It does **not** show
that the two can be mixed inside one training set, which is what closing the
cohort with the unfolded structures would require. That gate is still the
held-out-features test in `docs/PROVENANCE.md`.
