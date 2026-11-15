# prefold-dg-linear-probe

A linear reader of frozen Boltz-2 representations for protein–protein ΔΔG, and a
re-examination of the split it is measured on.

**Claim.** The ΔΔG signal the PreFold-dG head exploits is linearly readable from the
pooled Boltz-2 representations it is given. A ridge with one penalty scale per feature
block matches the 58.4M-parameter head on the split the field uses, and significantly
outperforms it once homology is controlled — on identical inputs, the same cohort and
the same folds.

## Results

Per-structure Pearson correlation: a correlation inside each complex, averaged over the
116 complexes carrying at least ten variants. Intervals are paired bootstraps over
complexes.

| split | external leakage | reader | head | reader − head |
|---|---|---|---|---|
| `rdenet` (what the field uses, grouped by PDB id) | 56% | 0.604 – 0.626 | 0.589 | +0.015 [−0.023, +0.055] |
| `hom99` (ours, 99% chain identity, 127 groups) | not measured | 0.573 | 0.505 | **+0.068 [+0.030, +0.109]** |
| `holdoutprot` (SKEMPI curators' own groups, 159 groups) | 0% | 0.579 | 0.524 | **+0.055 [+0.021, +0.088]** |

Leakage percentages are Bushuiev et al.'s iAlign measurements of these splits, not ours.
The `rdenet` reader entry is a range, not a point, and the difference beside it is
measured at its lower end: two runs of the same procedure on the same inputs selected
different penalties and landed at 0.604 and 0.626, which is the selection spread described
below. They are two points on one curve — `reproduce.py` prints 0.6031 at alpha 31.6 and
0.6264 at alpha 1 — and the difference quoted pairs the 0.604 run against the head's 0.589
on the same complexes. The two controlled rows do not have this ambiguity; there each
difference is the subtraction of the two columns beside it.

**Read the controlled rows.** `rdenet` groups by PDB id, which is not enough: at 99%
chain identity the cohort's 50 antibody–antigen PDB ids collapse to 22 independent
clusters, one of them fourteen hen-lysozyme complexes. The same antibody sits on both
sides of the split while the model is being fitted. Under control both models fall — so
the leakage is a property of split strictness in general, not something specific to
antibodies — but the head falls twice as far (−0.086 against the reader's −0.053). That
is the finding: **the leakage was helping the larger model more.**

The two controlled splits were built on unrelated principles — ours by clustering chain
sequences, the other by the SKEMPI curators' own `Hold_out_proteins` column, which
PPIformer also used and which Bushuiev et al. measured at 0% leakage. They agree on the
gap to within 0.013. `holdoutprot` is also the *looser* grouping of the two (159 groups
against 127), so the gap is not an artefact of over-strict grouping either.

## Reproduction of the published numbers

Before comparing anything, the published result had to be reproduced from the released
pickles alone:

| | this repository | PreFold-dG, as printed |
|---|---|---|
| SKEMPI per-structure PCC | 0.622 | 0.611 |
| SKEMPI per-structure SRCC | 0.597 | 0.573 |
| SKEMPI overall PCC | 0.744 | 0.750 |
| HER2 external test, PCC | 0.546 [0.483, 0.606] | 0.524 |
| HER2 external test, SRCC | 0.542 | 0.529 |

The HER2 protocol excludes nine SKEMPI entries homologous to HER2; that list was
re-derived from the rule in the paper's Supplementary Fig. S1 (Smith-Waterman against the
HER2 extracellular domain, BLOSUM62, gap open −12, extend −2, score above 350) and
reproduces its Supplementary Table S4 exactly. The published values lie inside our
interval, so this is a match, not a win — and the reader never sees HER2 at all, while
the published protocol selects a stopping epoch on it.

## How the reader works

Three blocks arrive glued together: `wsi` and `ws` at 147,456 numbers each, and `wz` at
128. Fit one ridge penalty to that concatenation and the two large matrices decide the
solution, burying `wz` — which is where most of the signal is. That costs about 0.09:

| | per-structure PCC |
|---|---|
| one shared penalty over the concatenation | ~0.53 |
| one penalty scale per block | 0.59 – 0.62 |

The fix is three lines. Standardise each block on the training rows, form its Gram
matrix, divide by the block's width, and sum:

```python
K = np.zeros((n, n))
for X in (wsi, ws, wz):                 # mutant minus wild type, per block
    z = (X - X[train].mean(0)) / X[train].std(0)
    K += (z @ z.T) / X.shape[1]         # <- the whole trick
```

Then solve the dual once and read off every penalty from the same eigendecomposition:

```python
w, V = np.linalg.eigh(K[train][:, train])
pred = K[test][:, train] @ (V @ ((V.T @ (y[train] - y[train].mean())) / (w + alpha))) + y[train].mean()
```

That is the entire model. No seed, no epochs, no early stopping, no model selection: the
same input and penalty give the same answer every time.

An earlier reading of the block penalty — that it protects weak blocks from strong ones —
was withdrawn. A control with a deliberately weak extra block showed the 0.049 is the
cost of the single shared penalty itself.

### Quote the gap, not a single score

The inner penalty-selection surface is flat near its top: 0.006 across a 300-fold range
of alpha, while the outer consequence of the choice is ten times larger. Nested
cross-validation makes it worse, because its inner folds train on 53% of the data against
the outer 67%, biasing selection toward over-regularisation. Three nominally identical
runs selected three different penalty sets and landed between 0.604 and 0.622.

So the reader's absolute score carries about ±0.02 of selection noise, and this
repository reports the whole penalty curve rather than a selected point. The
reader-minus-head difference does not have that problem: on `holdoutprot` it is positive
with intervals clearing zero at every grid point from alpha 1 to 1000.

### Verified against the reference implementation

`reproduce.py` in this repository, run on the released pickles, reproduces the working
repository's reference solver across the whole penalty grid:

| split | this repository | reference | penalty |
|---|---|---|---|
| `rdenet` | 0.6264 | 0.6264 | alpha 1 |
| `holdoutprot` | 0.5788 | 0.579 | alpha 10 |

All thirteen grid points agree to four decimal places on `rdenet`. The run takes under two
minutes on 14 cores once the pickles are loaded.

## Install

```bash
pip install -e .
```

NumPy, SciPy, pandas, scikit-learn and joblib. No GPU. The reader itself is about a
hundred lines.

## Usage

```bash
python reproduce.py --pickles <dir with skempi_{mut,wt}.pkl> --split rdenet
python reproduce.py --pickles <dir> --split holdoutprot
python -m prefold_probe.spec          # what the model is, and the experiment behind each choice
```

About 25 GB of RAM: the two 147,456-wide blocks are 6.8 GB in float32 for the full
cohort, and the kernel is formed in float64 because these representations have column
means large relative to their spread, and a float32 accumulation moves the third decimal.

## What the reader is given

| | |
|---|---|
| Required | `skempi_mut.pkl` and `skempi_wt.pkl` as released by the PreFold-dG authors, and a `PDB,fold` split file |
| Used | the mutant-minus-wild-type difference of `wsi`, `ws` and `wz` |
| Never used | a sequence model, an inverse-folding model, any structure beyond what Boltz-2 already pooled, HER2 at any point during fitting |

The label is `log10(Kd_mutant in nM) − log10(Kd_wild-type in nM)`, not kcal/mol — the
upstream code never converts, so published MAE and RMSE are in label units. Multiply by
1.364 for energies.

This repository reproduces the published loader's pairing exactly, including its one
known defect: it pairs a mutant with the first wild-type record sharing its four-letter
PDB id, so 35 variants in 2C5D, 3SE3 and 3SE4 carry the wrong wild type. Keeping it is
deliberate — fixing it would make these numbers incomparable with the published ones.

## What did not work

Measured, not abandoned early. Per-structure PCC is the figure of merit, and every
difference carries a paired complex-level bootstrap.

| approach | outcome |
|---|---|
| ESM-2 650M as an extra block | 0.596 against 0.615 for the reader alone; not significant |
| Pythia-PPI predictions as a block | 0.587 against 0.593; output averaging 0.591; not significant |
| richer inputs from the saved Boltz traces | response-ranked interface edges, distance shells, earlier blocks: at most +0.01 over a reader that reads all three tensors properly |
| other Pairformer layers and recycles | 88 snapshots scanned; flat after recycle 3–4; best +0.005 [−0.019, +0.031] over the final one |
| a nonlinear head on the same input | MLP minus ridge −0.002 [−0.021, +0.017] on the full data |
| head capacity, 0.1M to 58.4M parameters | every variant scores below the reader. A finite tuning grid cannot prove the converse, so the head grid was deliberately left unrun |
| combining with Pythia-PPI | residual correlation 0.765 — the two families fail on the same variants — and a gate that picks the better model per variant is not learnable (accuracy 0.47 against a 0.50 base rate) |
| distillation from a teacher's soft labels | +0.011 [+0.004, +0.018]; real but a regularisation effect |

## Where this problem stands

| route | status |
|---|---|
| a larger or better-tuned head | closed: capacity is not the limit, and the reader is never below the head under any condition tested |
| richer features from the same representations | closed: ≤ +0.01, below the cohort's resolution |
| fusing another learned representation | closed on this task: ESM-2 and Pythia-PPI each add nothing |
| more experimental data | out of reach here: Hummer et al. estimate ≥90,000 measured ΔΔG values are needed for a test Pearson above 0.85. SKEMPI has 5,731 |

Which leaves the representation itself. Everything here is frozen Boltz-2; nothing in
this work tests what happens when the trunk is allowed to move.

One subset is worth naming separately. Antibody–antigen variants score 0.47 against 0.64
for the rest, a gap of +0.17 to +0.18 that is not a composition artefact. Our pipeline is
handed the wild-type crystal structure as a template, so the difficulty is not structural
uncertainty. The candidate explanation is the representation: antibody CDRs are
somatically hypermutated and carry little of the evolutionary covariation an MSA-consuming
model relies on. That is cheap to test — correlate per-structure performance against MSA
depth and Boltz-2 confidence — and has not been done.

## How to read the numbers here

Rules this project arrived at the hard way; every point estimate believed without them was
later overturned.

- Per-structure PCC is interpretable only on the full cohort, where 116 complexes carry
  ten or more variants. Subset ablations report pcc / srcc / mae instead.
- Every difference comes with a paired bootstrap **over complexes**. Resampling variants
  holds the complex set fixed and understates the spread.
- A score on `rdenet` includes 56% leakage. Generalisation claims belong to `hom99` or
  `holdoutprot`.
- A single reader score carries ±0.02 of penalty-selection noise. Quote the gap.

## Data

This repository ships code and split files, not data. The feature pickles are the
PreFold-dG authors' release and are hosted at
[zenodo.org/records/20236663](https://zenodo.org/records/20236663) (6.4 GB for the mutant
file); pass the directory holding them to `reproduce.py`.

The two controlled splits are the part of this repository most likely to be useful
elsewhere:

| file | what |
|---|---|
| `data/split/hom99.csv` | three folds by GroupKFold over 99%-identity chain clusters, 127 groups |
| `data/split/holdoutprot.csv` | three folds by GroupKFold over SKEMPI 2.0's `Hold_out_proteins`, 159 groups |
| `data/split/rdenet.csv`, `ppiformer.csv` | the published splits, for comparison |

## Repository layout

| path | contents |
|---|---|
| `prefold_probe/reader.py` | the block-penalised kernel ridge |
| `prefold_probe/features.py` | assembling the three difference blocks from the released pickles |
| `prefold_probe/cohort.py` | the published loader's pairing, defect included |
| `prefold_probe/evaluate.py` | per-structure PCC and the paired complex bootstrap |
| `prefold_probe/splits.py` | fold assignments |
| `prefold_probe/spec.py` | what the model is and which measurement settled each choice |
| `reproduce.py` | fit on the released pickles and print the penalty curve |
| `RESULTS.md` | the full record, with the experiment behind every number |
| `docs/PROVENANCE.md` | which files every number is fitted on, with hashes, and what may be added to the cohort |

## A second baseline: `mmcsm-ppi/`

`mmcsm-ppi/` is a separate, self-contained harness for a different published method —
mmCSM-PPI (Rodrigues, Pires & Ascher, *NAR* 2021), which predicts the effect of *multiple*
point mutations from structure-based signatures rather than from a foundation model. It
scores any predictor against mmCSM-PPI, FoldX, Discovery Studio and ZEMu on the same rows,
using the predictions shipped inside the authors' own released files, so nothing has to be
reimplemented. Its README carries the full comparison and what the released files actually contain. On
mmCSM-PPI's own 595-mutation blind set, the block-penalised ridge scores Pearson 0.868 and
RMSE 1.41 kcal/mol against mmCSM-PPI's 0.687 and 2.06 — a paired difference of +0.166
[+0.065, +0.356] over 66 complexes — and a trained network on the same features does not
beat it there either. The two halves of this repository share no code and no cohort.

## Citation

The method and the data this builds on:

- Park et al. PreFold-dG. *Bioinformatics* 42(Supplement_2), btag489 (2026).
- Jankauskaitė et al. SKEMPI 2.0. *Bioinformatics* 35, 462–469 (2019).
- Passaro et al. Boltz-2. (2025).
