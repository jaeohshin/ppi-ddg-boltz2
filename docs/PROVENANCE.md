# Where the features came from

Every number in this repository is fitted on **features the PreFold-dG authors
released**, not on Boltz output produced here. That is a deliberate choice — it
is what makes the comparison against their published head a comparison of
readers rather than of pipelines — but it constrains what may be added to the
cohort later, and the constraint is not obvious from the code.

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

## Local Boltz settings, for the record

Any folding done on ada would currently use Boltz 2.2.1, checkpoint
`boltz2_conf.ckpt` sha256
`090e82ac8c92f5e943fa1b39e7410a44027bea7243c0bbb3caa67a77fc1428e1`, seed 42,
`recycling_steps` 3, per-sequence MSA cache. There is no evidence this matches
what produced the released pkls.
