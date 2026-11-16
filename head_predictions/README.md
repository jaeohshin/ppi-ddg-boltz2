# The head's predictions

The 58.4M-parameter PreFold-dG head is not in this repository and will not be: it is the
other paper's model. But the claim this repository makes is a *difference* — the reader
against that head, on identical inputs and folds — and a difference whose second term is
a number in prose cannot be checked.

So the head's out-of-fold predictions are here instead, as the retraining runs wrote
them. Nine files: three splits, three seeds each, 5,731 variants per file.

```
head_rdenet_seed{42,1337,2024}.csv.gz
head_hom99_seed{42,1337,2024}.csv.gz
head_holdoutprot_seed{42,1337,2024}.csv.gz
```

Columns: `pdb`, `mut_id`, `fold`, `true_ddg`, `pred_ddg`. Predictions, not scores, so
`compare.py` computes the per-structure correlation, the difference and the paired
interval over complexes rather than quoting any of them.

## The aggregation, which the table did not state

The `head` column of the README's Results table is a **three-seed ensemble**: the seeds'
predictions are averaged, and that average is scored once. Scoring each seed and
averaging the three scores is a different, lower number:

| split | ensemble (the published figure) | mean of seed scores | per-seed |
|---|---|---|---|
| `rdenet` | **0.5893** | 0.5832 | 0.5850, 0.5869, 0.5776 |
| `hom99` | **0.5049** | 0.5030 | 0.5046, 0.5039, 0.5005 |
| `holdoutprot` | **0.5235** | 0.5156 | 0.5103, 0.5235, 0.5131 |

The published table reports the first and did not say so, which makes the difference
against the reader smaller than the other reading would give — by 0.006, 0.002 and 0.008
on the three splits. `compare.py` prints both so the choice is visible rather than
inherited.

## Where they came from

These are the retraining runs described in RESULTS.md §4.2, not the authors' own
checkpoints: the head was refitted from scratch under each split, three seeds each, on
the same pickles and the same folds the reader uses. That is what makes the comparison
one of models rather than one of training protocols.

Their provenance outside this repository is the private experiment tree, under
`results/direct_delta_ablation/subtract/` for `rdenet` and
`analysis/reproduction_gap/results/head_{hom99,holdoutprot}/subtract/` for the other two.
Only the five columns above were carried over; the runs also recorded the separate
mutant and wild-type predictions, the checkpoint tag and the pair index, none of which
the comparison needs.

## Checking them

`SHA256SUMS` holds a digest of each **decompressed** file, and `compare.py` verifies the
ones it reads. To check by hand:

```bash
cd head_predictions
for f in *.csv.gz; do printf '%s  %s\n' "$(gzip -dc "$f" | sha256sum | cut -d' ' -f1)" "${f%.gz}"; done | diff - SHA256SUMS && echo OK
```

The scores in the table above are reproduced by `compare.py` from these files, so a
disagreement with it is a disagreement with something checkable.
