# mmCSM-PPI benchmark

A harness for comparing a multiple-mutation ΔΔG predictor against mmCSM-PPI on
mmCSM-PPI's own published test sets, on the same rows, using the predictions its
authors released.

The point is a head-to-head number: **does a new model beat 0.687 Pearson and 2.06
kcal/mol RMSE on the 595-mutation blind set?** Nothing here reimplements mmCSM-PPI, and
nothing here needs its code — its predictions ship inside its own data files.

## The numbers to beat

Recomputed here from the prediction columns of the released files, not copied from the
paper. Intervals are bootstraps over whole complexes.

**Blind test — 595 mutations of four or more point mutations, 68 complexes**

| method | Pearson | 95% CI | Spearman | RMSE | MAE |
|---|---|---|---|---|---|
| **mmCSM-PPI** | **+0.687** | [+0.456, +0.804] | +0.638 | **2.059** | 1.507 |
| FoldX | +0.395 | [+0.234, +0.518] | +0.370 | 5.277 | 4.128 |
| Discovery Studio | −0.390 | [−0.515, −0.262] | −0.419 | 5.575 | 4.274 |

**ZEMu comparison set — 272 mutations, 24 complexes**

| method | Pearson | 95% CI | Spearman | RMSE | MAE |
|---|---|---|---|---|---|
| **mmCSM-PPI** | **+0.732** | [+0.592, +0.870] | +0.753 | **1.730** | 1.310 |
| ZEMu | +0.636 | [+0.494, +0.787] | +0.651 | 2.134 | 1.641 |

Against the paper (Rodrigues, Pires & Ascher, *Nucleic Acids Res* 49:W417–W424, 2021):

| | here | paper |
|---|---|---|
| blind, mmCSM-PPI Pearson | 0.687 | 0.70 |
| blind, mmCSM-PPI RMSE | 2.059 | 2.02 in the body, **2.06 in the abstract** |
| blind, FoldX | 0.395 / 5.277 | 0.39 / 5.27 |
| ZEMu set, mmCSM-PPI | 0.732 / 1.730 | 0.73 / 1.72 |
| ZEMu set, ZEMu | 0.636 / 2.134 | 0.64 / 2.11 |
| blind, Discovery Studio | **−0.390** / 5.575 | +0.39 / 3.07 |

Everything lines up except Discovery Studio, whose released predictions are
*anti*-correlated with experiment as given. Negating them reproduces the paper exactly
(+0.390, RMSE 3.073), so the published figures are computed on the negated values. Which
side carries the sign error is not something these files can settle.

## One thing to know before claiming a win

The blind set is separated from the training set **by mutation count** — two or three
point mutations train, four or more test. The paper calls it "a non-redundant blind test
at mutation level", which is exactly accurate and easy to read as more than it says. At
the complex level it is not independent:

| | |
|---|---|
| blind complexes also in the released training set | **51 of 68 (75%)** |
| blind rows sitting on a trained complex | **525 of 595 (88%)** |

Split on that boundary, mmCSM-PPI scores 0.714 / RMSE 1.978 on complexes it trained on
and 0.672 / RMSE 2.590 on the 70 rows it did not. The correlation barely moves; the error
grows by 31%.

So the strongest form of a win is to **train without those 68 complexes and beat
mmCSM-PPI anyway** — it saw 88% of the test rows' complexes, and you saw none.
`data/blind_complexes.txt` lists them, marking which are in the training release.

## Use

```bash
python fetch_data.py --out data/release        # the three files, from the authors' site
python compare.py --release data/release       # the baselines alone
python compare.py --release data/release --predictions mine.csv --name ours
```

`mine.csv` needs three columns — `pdb`, `mutations`, `prediction`. Mutations are matched
as a set, so order does not matter, and rows you did not predict are reported rather than
dropped. `compare.py` prints each method's score with a complex-level interval, then the
paired difference against mmCSM-PPI with its own interval, and says plainly whether the
interval clears zero.

## Two traps in the released files

Both cost an afternoon here; `tests/test_data.py` pins them.

1. **`train_set` has an unquoted mutation list.** A triple-mutant row has one more field
   than a double-mutant row, so the columns must be read from the right.
2. **`train_set`'s `pdb` column is not always a PDB id.** Modelled reverse entries are
   named `1a22-db319a_eb320a_kb321a`. Counting that column gives 344 complexes where
   there are 130 — which is what the paper's 147-complex figure refers to, across all
   1,721 multiple mutations.

## Data

Not vendored. The site states no redistribution terms, so `fetch_data.py` pulls the three
files and the repository ships only code and derived facts. What is derived and kept:
`data/blind_complexes.txt`.

## Still to do

Deliberately out of scope for now; recorded so the choice is visible rather than
forgotten.

| | why it was left |
|---|---|
| **Homology clustering of the 68 blind complexes** | a different PDB id is not a different protein. Of the 17 blind complexes absent from training, an unknown number are homologues of training complexes — barnase/barstar, β-lactamase/BLIP and lysozyme/antibody families all recur in SKEMPI under several ids. Until they are clustered by sequence identity, "0.672 on unseen complexes" is an upper bound. Sequences for the 147 real PDB entries are downloadable; the clustering is single-linkage over chains at 99% and 90% identity |
| **Scoring the training set's own cross-validation** | the paper reports CV1 0.87 (random 80/20), CV4 0.76 (leave-one-complex-out) and CV5 0.73 (leave-one-binding-site-out). The released `prediction_mmCSM-PPI` column on the training file is not identified as out-of-fold, so it cannot stand in for any of them |
| **The abstract's 10-fold CV figure** | the abstract reports 0.75 / RMSE 1.64 under 10-fold CV. No such number appears in the body, whose closest entry is CV4 at 0.76 / 1.59. Unresolved |
| **Double versus triple, and per-complex scores** | not reported in the paper's body; possibly in supplementary tables S8–S11, unread |
| **Which side owns the Discovery Studio sign** | needs a Discovery Studio run, or the authors |

## Citation

Rodrigues CHM, Pires DEV, Ascher DB. mmCSM-PPI: predicting the effects of multiple point
mutations on protein–protein interactions. *Nucleic Acids Research* 49(W1):W417–W424
(2021). doi:10.1093/nar/gkab273
