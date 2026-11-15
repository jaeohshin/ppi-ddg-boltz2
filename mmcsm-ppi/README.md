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

## Result: the block-penalised ridge beats it

The reader described in the parent directory — a ridge over frozen Boltz-2 blocks with one
penalty scale each, fitted on the 2-to-3-mutation training rows — scores on 586 of the 595
blind mutations (66 of 68 complexes; nine rows had no feature row available).

| method | Pearson | 95% CI | Spearman | RMSE | MAE |
|---|---|---|---|---|---|
| **Boltz block ridge** | **+0.865** | [+0.769, +0.907] | **+0.813** | **1.436** | 1.119 |
| mmCSM-PPI | +0.687 | [+0.456, +0.804] | +0.638 | 2.059 | 1.507 |
| FoldX | +0.395 | [+0.234, +0.518] | +0.370 | 5.277 | 4.128 |
| Discovery Studio | −0.390 | [−0.515, −0.262] | −0.419 | 5.575 | 4.274 |

Paired against mmCSM-PPI on the same rows, bootstrapping whole complexes:

| | difference | 95% CI | complexes |
|---|---|---|---|
| Pearson | **+0.163** | [+0.062, +0.351] | 66 |
| Spearman | **+0.163** | [+0.055, +0.329] | 66 |

Both intervals clear zero. The error falls by 30%, from 2.06 to 1.44 kcal/mol.

On the 36 ZEMu-set rows that overlap: ridge +0.867 / 1.108 against mmCSM-PPI's +0.732 /
1.730 and ZEMu's +0.636 / 2.134, though ten complexes make that interval wide.

Reproduce it with the shipped predictions:

```bash
python compare.py --release data/release \
    --predictions data/predictions_boltz_block_ridge_cohort1237.csv \
    --name "Boltz block ridge"
```

### The number moves more with the penalty than with the data

Two fits of the same model on the same features are shipped, differing only in six
training rows (`cohort1237` is the complete cohort; `cohort1231` dropped the ambiguous
3SE3 rows — see below). Each selects its penalty by the same pre-registered rule: the
argmax of pooled inner out-of-fold Spearman over block weight × alpha, with the blind set
untouched. They do not select the same penalty.

| | training rows | selected | Pearson | RMSE |
|---|---|---|---|---|
| `cohort1231` | 1,231 | w=0.5, α=1 | 0.868 | 1.408 |
| `cohort1237` | 1,237 | w=1, α=3.16 | **0.865** | **1.436** |

Reading that as "six more rows made it worse" would be wrong. Holding the penalty fixed
and changing only the data isolates the two effects:

| penalty held fixed at | 1,231 rows | 1,237 rows | data effect |
|---|---|---|---|
| w=0.5, α=1 | 0.8683 | 0.8697 | **+0.0014** |
| w=1, α=3.16 | 0.8634 | 0.8648 | **+0.0014** |

The six rows help, by the same small amount at either penalty. The whole apparent
regression is the argmax moving: −0.0049 from the penalty, +0.0014 from the data.

So **0.5% more training data is enough to move the selected penalty**, and that move costs
3.5 times what the data gained. Scoring the blind set across the entire grid as a
diagnostic — not to select anything — the reachable range is 0.579 to 0.875, with the best
cell (w=2, α=1) at 0.8748 for the complete cohort, 0.010 above what the rule picked.

The honest reading: the margin over mmCSM-PPI is robust — every cell of the grid above
α≈0.1 clears it, and the paired interval excludes zero in both fits — but the third decimal
place of the headline is penalty-selection noise, not a property of the model. Quote the
gap, not the digits. The full grid for both cohorts is in
`data/penalty_sensitivity.json`, and the rule itself is `mmcsm_bench/disambiguate.py`
with `tests_mmcsm/test_disambiguate.py` pinning it on the real 3SE3 affinities. A selection rule stable enough to pin them down is an open problem
here, and it is the same instability that shows up as the `rdenet` range in the parent
directory's results.

### A trained network does not beat the ridge here either

The same features, the same folds, the same training rows, with the published
architecture's shared-head difference and an MSE on ΔΔG. Both sides are the 1,231-row
cohort, which is what makes this like for like; the ridge row is therefore 0.868, not the
0.865 of the complete cohort above:

| model | Spearman | Pearson | RMSE |
|---|---|---|---|
| **Boltz block ridge** | **0.815** | **0.868** | **1.408** |
| neural head, mean of 3 seeds | 0.770 | 0.823 | 1.632 |
| mmCSM-PPI | 0.650 | 0.702 | 2.035 |

Network minus ridge: ΔSpearman −0.045 [−0.115, +0.006]. This is the third cohort in which
a closed-form reader of these representations is not beaten by a trained head on the same
input — after SKEMPI single-point mutations under both the published and the
homology-controlled splits.

### What this result does not claim

- 586 of 595 test rows and 1,237 of the 1,340 training rows were usable; nine test rows,
  95 forward and 8 reverse training rows had no feature row. Most of the shortfall is 1KBH.
- **516 of the test rows share a PDB entry with the training set.** That is the paper's own
  mutation-count split, and mmCSM-PPI trained under the same condition, so the comparison is
  like for like — but neither number is a test of generalisation to a new complex.
- Overlap between these complexes and Boltz-2's pretraining was not audited.
- This benchmark had been looked at before the run, so it is not a fully blind confirmation.

## How far the margin reaches

A PDB id is not a protein family. Clustering the 147 complexes in these files by chain
sequence identity — single linkage, two entries joined when they share a chain cluster —
turns the 130 training entries into 73 families at 99% identity, 67 at 95% and 60 at 90%.
Grading the blind set by what the training set actually contained
(`analyse_homology.py`):

| the blind rows where… | rows | complexes | mmCSM-PPI | ridge | difference |
|---|---|---|---|---|---|
| everything | 595 | 68 | +0.687 | +0.865 | **+0.163 [+0.062, +0.352]** |
| the same PDB entry is in training | 525 | 51 | +0.714 | +0.875 | **+0.142 [+0.039, +0.351]** |
| the PDB entry is new | 70 | 17 | +0.672 | +0.753 | +0.082 [−0.029, +0.322] |
| **the family is new** | **19** | **7** | **+0.306** | **−0.082** | −0.388 [−1.481, −0.004] |

Pearson; the last two rows are at 99% identity, and 95% and 90% give 18 rows over 6
complexes with mmCSM-PPI at −0.010 and the ridge at −0.225.

**The margin shrinks monotonically as the test gets stricter**: +0.163 over everything,
+0.142 where the complex was trained on, +0.082 on a new PDB entry with an interval that
now includes zero, and nothing favourable where the family is new.

Two things follow, and they pull in opposite directions.

The comparison above is still **fair**: mmCSM-PPI was trained under exactly the same
condition, on the same mutation-count split, so on the paper's own protocol the ridge wins
it, and that is what the protocol was built to measure.

But the win is **earned almost entirely on complexes the model had already seen**. Where
the family is new there are 19 rows over 7 complexes — far too few to conclude anything,
with an interval spanning 1.5 — and on that little evidence **both methods fail**, mmCSM-PPI
at +0.31 and the ridge below zero (−0.08). Neither model has been shown to transfer to a new
protein family, because neither has been tested on one. Any claim that stops at the first
row of that table is claiming more than these files can support.

## Why 112 rows are still dropped, and what it would take to get them back

The ridge scored 586 of the 595 blind mutations and trained on 1,237 of the 1,340 released
rows. `diagnose_missing.py` classifies every dropped row; it needs the local feature
cohort, so it is a record of the diagnosis rather than something a reader can re-run.

| | rows | what it is |
|---|---|---|
| complex never folded | 110 | 1KBH 86, 3VR6 14, 4GNK 4, 2NYY 3, 2NZ9 2, 4GXU 1 — no Boltz features exist for these complexes at all |
| several interfaces, resolved by label | 6 | 3SE3 — **recovered**, see below; these are in the 1,237-row cohort and absent from the 1,231-row one |
| complex folded, this mutant was not | 2 | 2C5D, blocked on a defect in SKEMPI, see below |

**The six 3SE3 rows have been recovered.** 3SE3 has two interfaces, `3SE3_B_A` and
`3SE3_B_C`, and a mutation on chain B appears under both with very different measured
effects — ΔΔG −2.458 against −0.191 kcal/mol for the same three substitutions. mmCSM-PPI's
released files name the complex by its four-letter id only, so nothing in the row says
which interface it is, and the preparation code, which required exactly one candidate
feature row, skipped all six.

The row's own experimental value settles it. That value is a measurement, not a
prediction, and it reproduces the SKEMPI affinities behind each candidate exactly, once
the sign convention is handled: mmCSM-PPI's ΔΔG is ln(K_wt) − ln(K_mt), the negative of
RT·ln(Kd_mut/Kd_wt).

| row | label | interface | residual |
|---|---|---|---|
| 1038 | −0.190819 | `3SE3_B_C` | 0.0000 |
| 1039 (reverse) | +0.190819 | `3SE3_B_C` | 0.0000 |
| 1042 | −0.173393 | `3SE3_B_C` | 0.0000 |
| 1043 (reverse) | +0.173393 | `3SE3_B_C` | 0.0000 |
| 1046 | −2.457730 | `3SE3_B_A` | 0.0000 |
| 1048 | −0.769763 | `3SE3_B_A` | 0.0000 |

An **exact** match is required, to 1e-4, and exactly one candidate must match; a row that
matches none or several stays dropped. That makes this recovery of a row's identity rather
than a choice of the closest fit, which matters because the discriminant is the label. The
test set is untouched at 586 rows, so the blind comparison is scored on exactly the same
rows as before and the only difference is six training rows.

**The two 2C5D rows are not a parsing problem.** SKEMPI 2.0's own `Mutation(s)_PDB` field
carries a duplicated token for them:

```
2C5D   RA310E,KA312E,RB310E,KA312E      KA312E twice
2C5D   EC56R,EC59R,ED56R,EC59R          EC59R twice
```

mmCSM-PPI copied the field faithfully, so the defect is upstream of both. Whether the
intended set is K312E on chain A twice or on chains A and B cannot be read off the file,
and a mutant cannot be folded until it is settled. These two rows are blocked on a
definition, not on compute.

**The remaining 110 need 118 structures folded** — the six absent complexes with their
wild types, plus the two 2C5D mutants if their mutation set is ever resolved.
`data/to_fold.csv` lists them. 1KBH alone is 86 of the 110, so folding one complex
recovers three quarters of the shortfall and takes the training set from 1,237 to 1,323.

Folding is not free, and not because of compute. Every fitted feature here comes from the
PreFold-dG authors' released pickles; locally folded structures would be a second producer
in one training set, and the ridge would fit the mixture without complaint.
[`docs/PROVENANCE.md`](../docs/PROVENANCE.md) gives the hashes, how far a local re-fold
agrees, why feature cosine is not a valid acceptance test, and the held-out-features
comparison that would have to pass first.

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
| **A split by complex, not by mutation count** | the comparison above is fair but neither model is tested on an unseen complex. Refitting both on a complex-level or family-level split is the experiment that would show whether the margin survives |
| **Closing the cohort** | diagnosed below. The six ambiguous rows are recovered; 110 still need 118 structures folded; two are blocked on a defect in SKEMPI itself. Folding is not a free operation: the fitted features are the PreFold-dG authors' released pickles, so new local features would be a second producer in one training set. [`docs/PROVENANCE.md`](../docs/PROVENANCE.md) gives the hashes, how far local re-folds agree, why cosine is not an acceptance test, and the held-out-features test to run first |
| **A benchmark with new families in it** | done and reported above: only 7 of the 68 blind complexes are a family the training set never saw, carrying 19 rows. No protocol over these files can measure out-of-family transfer, because the data for it is not there. Testing that needs a different cohort |
| **Scoring the training set's own cross-validation** | the paper reports CV1 0.87 (random 80/20), CV4 0.76 (leave-one-complex-out) and CV5 0.73 (leave-one-binding-site-out). The released `prediction_mmCSM-PPI` column on the training file is not identified as out-of-fold, so it cannot stand in for any of them |
| **The abstract's 10-fold CV figure** | the abstract reports 0.75 / RMSE 1.64 under 10-fold CV. No such number appears in the body, whose closest entry is CV4 at 0.76 / 1.59. Unresolved |
| **Double versus triple, and per-complex scores** | not reported in the paper's body; possibly in supplementary tables S8–S11, unread |
| **Which side owns the Discovery Studio sign** | needs a Discovery Studio run, or the authors |
| **A penalty-selection rule that does not move with 0.5% of the data** | measured above: six training rows change the selected penalty, and that change costs 3.5x what the rows gained. The margin survives it, the third decimal does not. Repeated-fold or bootstrap-averaged inner selection is the obvious next thing to try, and it would also close the `rdenet` range in the parent directory |

## Citation

Rodrigues CHM, Pires DEV, Ascher DB. mmCSM-PPI: predicting the effects of multiple point
mutations on protein–protein interactions. *Nucleic Acids Research* 49(W1):W417–W424
(2021). doi:10.1093/nar/gkab273
