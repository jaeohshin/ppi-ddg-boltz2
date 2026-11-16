# ppi-ddg-boltz2

A linear reader of frozen Boltz-2 representations for protein–protein ΔΔG, and a
re-examination of the split it is measured on.

**Claim.** The ΔΔG signal the PreFold-dG head exploits is linearly readable from the
pooled Boltz-2 representations it is given. A ridge with one penalty scale per feature
block matches the 58.4M-parameter head on the split the field uses, and significantly
outperforms it once homology is controlled — on identical inputs, the same cohort and the
same folds. Of the 295,040 numbers it is handed, 128 do the work.

## Where to read what

This README is an overview. It is not the primary record of anything below.

| | |
|---|---|
| `RESULTS.md` | the full record: every experiment in the order it was run, with the measurement that settled it, and a section for the ones withdrawn |
| `docs/FAMILY_SPLIT.md` | the mmCSM-PPI comparison with similarity controlled — what the margin was worth, and which 128 features did the work |
| `docs/PAIR_BLOCK.md` | inside the 128: four dimensions are enough, and the strong ones are copies of each other |
| `docs/SELECTION.md` | how the penalty is chosen, four rules compared, and the two things that are not settled |
| `docs/PROVENANCE.md` | which files every number is fitted on, with hashes, and what may be added to the cohort |
| `reproduction/` | refitting on our own Boltz-2 inference instead of the authors' pickles, and what that test does and does not license |
| `mmcsm-ppi/` | a separate, self-contained harness scoring any predictor against mmCSM-PPI, FoldX, Discovery Studio and ZEMu on the same rows |
| `VERIFICATION.md` | what is checked before a number here is believed |

## Results

Per-structure Pearson correlation: a correlation inside each complex, averaged over the
116 complexes carrying at least ten variants. Intervals are paired bootstraps over
complexes.

| split | external leakage | reader | head | reader − head | **pair block alone** |
|---|---|---|---|---|---|
| `rdenet` (what the field uses, grouped by PDB id) | 56% | 0.604 – 0.626 | 0.589 | +0.015 [−0.023, +0.055] | 0.589 |
| `hom99` (ours, 99% chain identity, 127 groups) | not measured | 0.573 | 0.505 | **+0.068 [+0.030, +0.109]** | **0.578** |
| `holdoutprot` (SKEMPI curators' own groups, 159 groups) | 0% | 0.579 | 0.524 | **+0.055 [+0.021, +0.088]** | **0.576** |
| `ppiformer` (PPIRef's similarity-aware split) | — | — | — | — | 0.573 |

The last column is 128 features; the `reader` column is 295,040. Subtracting it from the
`head` column, a closed-form ridge on 128 frozen numbers beats the 58.4M-parameter head by
+0.073 on `hom99` and +0.052 on `holdoutprot`, and lands exactly on it on the leaked split.
`RESULTS.md` §3.1b refits every block subset from scratch on its own penalty curve, and has
the paired intervals against the full reader: the 294,912 sequence dimensions are worth a
significant +0.038 on the leaked split and nothing distinguishable from zero on either
controlled one. Leakage percentages are Bushuiev et al.'s iAlign measurements of these
splits, not ours.

**Read the controlled rows.** `rdenet` groups by PDB id, which is not enough: at 99% chain
identity the cohort's 50 antibody–antigen PDB ids collapse to 22 independent clusters, one
of them fourteen hen-lysozyme complexes. Under control both models fall, so the leakage is
a property of split strictness in general rather than something specific to antibodies —
but the head falls twice as far (−0.086 against the reader's −0.053). That is the finding:
**the leakage was helping the larger model more.** The two controlled splits were built on
unrelated principles and agree on the gap to within 0.013.

The `rdenet` reader entry is a range because two runs of the same procedure on the same
inputs selected different penalties and landed at 0.604 and 0.626. A single reader score
carries about ±0.02 of penalty-selection noise; the difference does not. `docs/SELECTION.md`
has why, and `RESULTS.md` §4 has the controlled refit in full.

The published numbers were reproduced from the released pickles before anything was
compared to them — SKEMPI per-structure PCC 0.622 against 0.611 as printed, HER2 external
test 0.546 [0.483, 0.606] against 0.524, with the paper's nine-entry homology exclusion
re-derived from its Supplementary Fig. S1 rule. The published values lie inside our
intervals, so this is a match, not a win. `RESULTS.md` §2 and §3.3.

## The model

Three blocks arrive glued together: `wsi` and `ws` at 147,456 numbers each, and `wz` at
128. Fit one ridge penalty to that concatenation and the two large matrices decide the
solution, burying `wz` — which is where most of the signal is. That costs about 0.09. The
fix is three lines: standardise each block on the training rows, form its Gram matrix,
divide by the block's width, and sum.

```python
K = np.zeros((n, n))
for X in (wsi, ws, wz):                 # mutant minus wild type, per block
    z = (X - X[train].mean(0)) / X[train].std(0)
    K += (z @ z.T) / X.shape[1]         # <- the whole trick

w, V = np.linalg.eigh(K[train][:, train])
pred = K[test][:, train] @ (V @ ((V.T @ (y[train] - y[train].mean())) / (w + alpha))) + y[train].mean()
```

That is the entire model: one eigendecomposition serves the whole penalty grid. No seed,
no epochs, no early stopping — the same input and penalty give the same answer every time.
`reproduce.py` matches the working repository's reference solver to four decimal places at
all thirteen grid points on `rdenet`. `RESULTS.md` §3.2 has the block-penalty experiment,
including an earlier reading of it that was withdrawn.

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

About 25 GB of RAM: the two 147,456-wide blocks are 6.8 GB in float32 for the full cohort,
and the kernel is formed in float64 because these representations have column means large
relative to their spread, and a float32 accumulation moves the third decimal.

## What the reader is given

| | |
|---|---|
| Required | `skempi_mut.pkl` and `skempi_wt.pkl` as released by the PreFold-dG authors, and a `PDB,fold` split file |
| Used | the mutant-minus-wild-type difference of `wsi`, `ws` and `wz` |
| Never used | a sequence model, an inverse-folding model, any structure beyond what Boltz-2 already pooled, HER2 at any point during fitting |

Two things to know about the inputs. The label is
`log10(Kd_mutant in nM) − log10(Kd_wild-type in nM)`, not kcal/mol — the upstream code
never converts, so published MAE and RMSE are in label units; multiply by 1.364 for
energies. And this repository reproduces the published loader's pairing exactly, including
its one known defect: it pairs a mutant with the first wild-type record sharing its
four-letter PDB id, so 35 variants in 2C5D, 3SE3 and 3SE4 carry the wrong wild type.
Keeping it is deliberate — fixing it would make these numbers incomparable with the
published ones.

## How to read the numbers here

Rules this project arrived at the hard way; every point estimate believed without them was
later overturned.

- Per-structure PCC is interpretable only on the full cohort, where 116 complexes carry ten
  or more variants. Subset ablations report pcc / srcc / mae instead.
- Every difference comes with a paired bootstrap **over complexes**. Resampling variants
  holds the complex set fixed and understates the spread.
- A score on `rdenet` includes 56% leakage. Generalisation claims belong to `hom99` or
  `holdoutprot`.
- A single reader score carries ±0.02 of penalty-selection noise. Quote the gap.

## What is closed, and what is open

Eight alternatives were measured rather than abandoned early, each with a paired
complex-level bootstrap: ESM-2 650M as an extra block, Pythia-PPI as a block or by output
averaging, richer inputs from the saved Boltz traces, other Pairformer layers and recycles,
a nonlinear head, head capacity from 0.1M to 58.4M parameters, combining with Pythia-PPI
through a learned gate, and distillation from soft labels. None beats the reader by more
than this cohort can resolve. `RESULTS.md` §5 lists every one with its interval.

| route | status |
|---|---|
| a larger or better-tuned head | closed: capacity is not the limit, and the reader is never below the head under any condition tested |
| richer features from the same representations | closed: ≤ +0.01, below the cohort's resolution |
| fusing another learned representation | closed on this task: ESM-2 and Pythia-PPI each add nothing |
| more experimental data | out of reach here: Hummer et al. estimate ≥90,000 measured ΔΔG values are needed for a test Pearson above 0.85. SKEMPI has 5,731 |

Which leaves the representation itself. Everything here is frozen Boltz-2; nothing in this
work tests what happens when the trunk is allowed to move.

One subset is worth naming separately. Antibody–antigen variants score 0.47 against 0.64
for the rest, a gap that is not a composition artefact, and our pipeline is handed the
wild-type crystal structure as a template, so the difficulty is not structural uncertainty.
The candidate explanation is the representation: antibody CDRs are somatically hypermutated
and carry little of the evolutionary covariation an MSA-consuming model relies on. That is
cheap to test — per-structure performance against MSA depth and Boltz-2 confidence — and
has not been done. `RESULTS.md` §4.1 and §7.

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

## Code

| path | contents |
|---|---|
| `prefold_probe/reader.py` | the block-penalised kernel ridge |
| `prefold_probe/features.py` | assembling the three difference blocks from the released pickles |
| `prefold_probe/cohort.py` | the published loader's pairing, defect included |
| `prefold_probe/evaluate.py` | per-structure PCC and the paired complex bootstrap |
| `prefold_probe/splits.py` | fold assignments |
| `prefold_probe/spec.py` | what the model is and which measurement settled each choice |
| `reproduce.py` | fit on the released pickles and print the penalty curve |
| `block_ablation.py` | refit on block subsets; which of the three blocks generalises |
| `pair_block.py` | extract the 128-wide pair block once (2.6 MiB) and score it on all four splits |
| `pair_block_dimensions.py` | how many of the 128 carry the signal, with the ranking recomputed per fold |

The import path is `prefold_probe`, from the repository's earlier name.

## A second baseline: `mmcsm-ppi/`

`mmcsm-ppi/` scores any predictor against mmCSM-PPI (Rodrigues, Pires & Ascher, *NAR*
2021), which predicts the effect of *multiple* point mutations from structure-based
signatures rather than from a foundation model. It uses the predictions shipped inside the
authors' own released files, so nothing has to be reimplemented. On mmCSM-PPI's own
595-mutation blind set the ridge scores Pearson 0.865 against 0.687, a paired difference of
+0.163 [+0.062, +0.351] over 66 complexes.

That margin is a property of the split, and the split is lenient: 516 of the 586 scorable
rows sit on complexes the training release contains. Refitted leave-one-family-out, the
ridge falls to 0.555 against mmCSM-PPI's 0.702 — but mmCSM-PPI saw every held-out family
while the ridge did not, so the difference is undetermined rather than a loss. The same run
showed that only the pair block generalises: dropping both 147,456-wide sequence blocks
costs 0.0011 Pearson once the family is new. `docs/FAMILY_SPLIT.md` has it, and
`mmcsm-ppi/README.md` has the full comparison and what the released files actually contain.
The two halves of this repository share no code and no cohort.

## Citation

The method and the data this builds on:

- Park et al. PreFold-dG. *Bioinformatics* 42(Supplement_2), btag489 (2026).
- Jankauskaitė et al. SKEMPI 2.0. *Bioinformatics* 35, 462–469 (2019).
- Passaro et al. Boltz-2. (2025).
