# Where the features came from

Every number in this repository is fitted on **features the PreFold-dG authors
released**, not on Boltz output produced here. That is a deliberate choice — it
is what makes the comparison against their published head a comparison of
readers rather than of pipelines — but it constrains what may be added to the
cohort later, and the constraint is not obvious from the code.

The released fits still use the authors' features. What that choice costs is now
measured rather than assumed: swapping in the author-matched local extraction
moves the blind-set Pearson by +0.003 [−0.032, +0.022]. See
[The swap test](#the-swap-test-run-2026-11-15).

## The pinned inputs

From `prepared_v1/config.json` of the run the released fits descend from:

| role | file | sha256 | bytes |
|---|---|---|---|
| mutant representations | `skempi_mut.pkl` | `f89501f2ae2bca19e3725c8e1fc050002f02ac7094f7a4b1b5fe96de3ef90e91` | 6,443,150,681 |
| wild-type representations | `skempi_wt.pkl` | `a106fb130d91ffbca2826e685f7cf7356680a9fb2fe81d85c208c883d88bb0b2` | 375,498,913 |
| cohort | `oof_predictions.csv` | `2e8be68020a0439d959355a0b01b9b43f1a6ca610c0b21a2e114d3d613473ce8` | 1,073,111 |
| labels | `skempi_v2.csv` | `76f8e60d09eaa4f05fb265063886859a86895c0b2fb225e31fe82928a883972d` | 1,602,208 |

Field widths 147,456 / 147,456 / 128 (`wsi`, `ws`, `wz`). Labels are
`log10(Kd_MT/Kd_WT)` — multiply by 1.364 for kcal/mol.

The two pkls are the authors' Zenodo release. The Boltz version, checkpoint and
sampling settings behind them are **not recorded in that release**, and were not
recovered.

## Why that matters for closing the cohort

The 118 structures listed in `mmcsm-ppi/data/to_fold.csv` have no author
features, because the authors never folded them. Folding them here and
concatenating would produce a training set that is part author features and part
local features. The ridge would fit it without complaint.

Two measurements bound how different those halves would be.

**Matching the authors' settings matters a great deal at the feature level.**
Median cosine of the `mut − wt` difference against the released pkls, over 5,471
matched variants, for a local re-fold with and without the wild-type PDB as a
template:

| field | with template | template-free |
|---|---|---|
| `wsi` | 0.948 | 0.905 |
| `ws`  | 0.923 | 0.855 |
| `wz`  | **0.953** | **0.577** |

Template-free output is not in the same space as the authors'. `wz` is the field
the ridge leans on hardest.

**But cosine is not an acceptance test.** Twice now, moving local features
*closer* to the authors' failed to help, once moving the opposite way:

- the 2026-09-15 template ablation raised feature cosine 0.98 → 0.999 and gave
  no `group_pcc` gain;
- switching to the paper's MSA setting raised cosine on *every* field
  (+0.020 to +0.061) while moving `group_pcc` **−0.10**.

Feature agreement with the published pkls must not be used as a proxy target.
The acceptance test is a score on held-out complexes, and one has now been run —
see [The swap test](#the-swap-test-run-2026-11-15).
The likely mechanism: this project caches MSAs per sequence, so an unmutated
chain gets a byte-identical alignment in the WT and mutant runs and `mut − wt`
carries no MSA sampling noise. A pipeline that regenerates the MSA per complex
injects noise into exactly the difference the reader keys on. That is a
difference in the *noise structure of the difference vector*, which no cosine on
the vectors themselves will reveal.

One robust observation on top: scored with the authors' own head, every local
feature variant sits below the authors' features (pcc 0.80 vs 0.84, srcc
0.63–0.65 vs 0.69, mae 0.86 vs 0.79) — a gap larger than any difference between
local variants. At cosine 0.975 the head is still reading a measurably different
signal.

## The test to run before mixing anything

Not a cosine. A held-out-features test:

1. pick complexes the authors **did** fold,
2. re-fold them locally with the intended settings (wild-type PDB as template),
3. fit the ridge on the authors' features with those complexes held out,
4. score the held-out rows twice — once on the authors' features, once on the
   local re-fold.

If the two scores are indistinguishable under the paired complex-level
bootstrap, local features may be appended. Any other outcome means the cohort
stays as it is and the 110 rows stay missing.

If features are ever mixed, the row's producer belongs in the cohort file as a
column, so a later reader can split on it.

## Local extractions, and which one matches the authors

Several local extractions exist in `prefold-dg-mechanism/processed_data/`. They
are **not interchangeable**, and the difference between them is larger than the
difference between the matched one and the authors' release.

| extraction | recycling | template | status |
|---|---|---|---|
| `scaled_v1_r10_tplfull` | 10 (the paper's) | every WT complex's deposited structure | **the author-matched extraction. Use this one.** |
| `scaled_v1_recycle10{,_fixed}` | 10 | none | not matched |
| `scaled_v1_recycle3` | 3 (Boltz default) | none | not matched |
| `scaled_v1_{own,paper}msa574`, `notemplate600`, `template600`, `r10_tplA`, `r10_tplpmsa574` | various | various | MSA / template ablations |

Common to all of them: Boltz 2.2.1, checkpoint `boltz2_conf.ckpt` sha256
`090e82ac8c92f5e943fa1b39e7410a44027bea7243c0bbb3caa67a77fc1428e1`, seed 42,
per-sequence MSA cache. The author confirmed by correspondence that this
checkpoint is bit-identical to theirs (sha1 `ef15d2d4...`) and that MSA pairing is
the ColabFold greedy default, so no weight or MSA-scheme difference exists.
Recycling 3 is the Boltz default and recycling 10 is the paper's; both were run.

**Do not read "recycling 3" as the matched setting.** An earlier version of this
section recorded recycling 3 as the project's configuration without saying that
it is unmatched, and that wording caused a wrong-variant comparison to be
pre-registered and run on 2026-11-15.

## The swap test, run 2026-11-15

The whole feature producer was replaced and the mmCSM-PPI blind comparison refit,
on the 440 of 586 blind rows that join between the two producers (59 of 66
complexes), both arms using identical rows and the authors' labels:

| features | Pearson | rho | vs authors' features | vs mmCSM-PPI |
|---|---|---|---|---|
| authors' Zenodo pkls | 0.8723 | 0.8185 | — | +0.117 [+0.026, +0.342] |
| local `scaled_v1_r10_tplfull` | **0.8751** | 0.8165 | **+0.003 [−0.032, +0.022]** | **+0.119 [+0.022, +0.344]** |
| local `scaled_v1_recycle3` | 0.8187 | 0.7347 | −0.054 [−0.131, −0.029] | +0.063 [−0.049, +0.279] |

With the matched extraction the released headline is reproducible from local
Boltz output and the margin over mmCSM-PPI is unchanged. With an unmatched one it
is not. This is a producer-swap test, not the held-out-features test below, and it
does **not** license mixing the two inside one training set.

Run and scripts: `/data/work/local_refold_recycle3/` (`run_variant.sh`,
`build_join.py`, `compare_sources.py`, `diagnose.py`), summary in
`run/SUMMARY.md`.

## Why the extraction matters more than its raw fidelity suggests

The reader never sees Boltz output; it sees `mut − wt`, a small residual of two
nearly identical vectors. Subtraction amplifies any disagreement, and most of
all in the block where the residual is the smallest fraction of the vector
(`scaled_v1_recycle3` arm, 400 blind-cohort rows):

| block | ‖Δ‖ / ‖wt‖ | raw-vector cosine | difference cosine | angle amplification |
|---|---|---|---|---|
| `wsi` | 0.115 | 0.9867 | 0.9392 | 2.1× |
| `ws`  | 0.105 | 0.9925 | 0.8567 | 4.4× |
| `wz`  | **0.048** | 0.9929 | **0.6343** | **7.2×** |

`wz` — the only block that generalises across protein families — carries a
residual that is 4.8% of the vector norm, so it is amplified hardest. Raw output
matched to cosine 0.993 still leaves the quantity the model reads disagreeing at
0.63. The matched extraction reaches raw 0.9998 and difference 0.9764 on `wz`,
and that is the whole difference between the two local rows in the table above.

This is why repeated attempts to reproduce the authors' inference by tightening
raw-output agreement did not converge: raw fidelity is the wrong target.

`prefold-dg-mechanism` bounds the floor. Rerunning one pipeline unchanged across
120 GPU tasks (2026-10-06) agrees with itself at tensor-level cosine
1.00000–1.00005 on both mutant and wild-type tensors — deterministic to ~1e-5
across hardware — yet self-agreement on the *difference* vector is only
0.9939 / 0.9976 / 0.9963 (`wsi`/`ws`/`wz`). Author-vs-local on the same variants
is 0.9748 / 0.9859 / 0.9665, a median gap of −0.008 / −0.003 / −0.006: inside the
noise floor. No global setting is missing.

## What is still open

One documented difference has never been replicated. The author volunteered it
and it appears in neither the paper nor their `process_data.ipynb`: *"where the
RCSB PDB sequence and the SKEMPI DB sequence disagreed, I preprocessed by
aligning and then filling in or trimming."* A different input sequence changes
one target's features while leaving every global setting alone, which is the
shape of the per-target residual that survives on a quarter to a half of
variants. Asking for the alignment rule would close it.

## Joining the two producers

They key records differently, which has to be solved before anything can be
compared. The authors' cohort uses SKEMPI's `Mutation(s)_cleaned` positions with
the token list sorted; the local extractions use `Mutation(s)_PDB` positions in
SKEMPI's original order. SKEMPI prints both columns in the same token order, so
zipping them yields a per-(complex, chain) position dictionary — 4,008 entries,
no conflicts — which normalises local keys into the authors' convention. That
joins 5,326 of 5,729 mutant records and 335 of 336 wild-type records. Eight
antibody complexes differ only in chain-group spelling (`1DAN_HL_TU` vs
`1DAN_HL_UT`) and are matched by trying the permutations.

The join is validated against the labels stored inside both pkls independently of
the features: Pearson 0.988, median difference 0.005, the residual explained by
SKEMPI duplicate-affinity aggregation. Labels are never taken from the local
side; both arms use the cohort's.
