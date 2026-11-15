# Choosing the penalty, and why it was harder than it looked

The block ridge has two hyperparameters: one scale for the two sequence blocks
relative to the pair block, and one ridge penalty. Both are chosen by
cross-validation inside the training rows, with the blind set untouched. This
records what happened when that choice was examined properly, because the answer
changed the status of the headline number.

## What started it

Recovering six ambiguous training rows — 0.5% of the cohort — moved the selected
penalty from (weight 0.5, alpha 1) to (weight 1, alpha 3.16) and the blind
Pearson from 0.8697 to 0.8648. Holding the penalty fixed, those six rows were
worth **+0.0014 at either setting**. The entire apparent regression was the
argmax moving: −0.0049 against +0.0014.

So the selection was carrying more than three times the signal the data carried,
and the single 3-fold split that produced it decides everything on one partition
of 73 PDB groups.

## Four rules, all selecting on inner out-of-fold Spearman

Same cohort, features, grid and reverse-row handling throughout; only the inner
split differs. Each rule's blind score is at the penalty that rule selected,
scored once, after its selection was closed.

| inner split | selected | inner ρ | blind ρ | blind r | RMSE |
|---|---|---|---|---|---|
| single 3-fold, grouped by PDB (1,231 rows) | w=0.5, α=1 | — | 0.8197 | 0.8697 | 1.402 |
| single 3-fold, grouped by PDB (1,237 rows) | w=1, α=3.16 | — | 0.8128 | 0.8648 | 1.436 |
| 10×3 repeats, grouped by PDB | w=0.5, α=3.16 | 0.702 | 0.8014 | 0.8565 | 1.482 |
| 10×3 repeats, row-level | w=4, α=0.1 | 0.885 | 0.8213 | 0.8424 | 1.618 |
| **30 Monte-Carlo, mixture-matched** | **w=4, α=0.32** | 0.865 | **0.8301** | 0.8573 | 1.507 |
| *grid optimum on ρ, for reference* | *w=4, α=1* | — | *0.8352* | *0.8713* | *1.404* |

### Repeating the split did not help, and that was informative

Averaging 10 repeats of the grouped split gave a stable argmax and a **worse**
blind score. The per-repeat argmaxes had been spread over six distinct grid
cells, three of them at most — so no single grouped split could have pinned the
penalty — but the stabilised answer, α = 3.16, was further from the optimum than
either single split had landed.

That is a bias, not variance. The whole repeat-averaged inner surface prefers
α ≈ 3.16–10 at every weight, while the blind surface prefers α ≈ 1. More repeats
cannot remove an offset.

### The two splits miss in opposite directions

Grouping by PDB means inner validation rows sit on complexes the inner training
set never saw. The blind set is not like that: **516 of its 586 rows share a PDB
entry with the training release**, because mmCSM-PPI split by mutation count.
Grouped inner CV estimates new-complex generalisation; the blind set measures
mostly seen-complex accuracy. Harder task, more regularisation preferred.

A row-level split — complexes free to appear on both sides — was run to test
that. Its α fell to 0.1, confirming the direction, but its blind Pearson fell
too, which the prediction had got wrong. The reason: it also pushed the weight
from 0.5 to 4. Letting a complex appear on both sides makes the sequence blocks
look valuable for memorising it, so the rule cranks them up and releases the
penalty.

| | α | inner ρ | actual blind ρ | inner estimate is |
|---|---|---|---|---|
| grouped | 3.16 | 0.702 | 0.801 | **0.119 pessimistic** |
| row-level | 0.1 | 0.885 | 0.821 | **0.064 optimistic** |

### Matching the mixture

The blind set's condition is intermediate: 70 of 586 rows on complexes absent
from training, 516 on complexes present. The last rule reproduces that — each
Monte-Carlo split holds some complexes out of training entirely until their rows
reach 11.95% of the validation set, and fills the rest from complexes that stay
in training, with every shared complex keeping at least one training row.

That 11.95% is counted from the released files' composition. No label or score
enters it.

It is the best of the rules on the criterion all of them optimise: blind ρ
0.8301 against a grid optimum of 0.8352, and it recovers the right weight (4),
which the grouped rule had wrong by a factor of eight. Argmax concentration rose
from 3/10 to 9/30.

## Two things this does not settle

**The blind set has now been scored five times.** Everything after the first two
runs was designed while the blind surface
(`mmcsm-ppi/data/penalty_sensitivity.json`) was known. The mixture rule is
justified by a stated principle, and its target comes from public row counts
rather than labels, but a rule designed with the answer in view is not confirmed
by that answer. **The published headline is therefore left at the pre-registered
single-split run**, and the mixture rule will be judged where no penalty has been
chosen yet: the family-level split.

**The protocol selects on Spearman and reports Pearson.** These rank the penalty
differently — the mixture rule is first on ρ and fourth on r; the single-split
rule is the reverse. That is inherited from the original run script and it is a
defect, not a detail. Either the selection criterion should be the reported
metric, or both should be reported with the selected one named. Until that is
settled, the only safe reading of any of these numbers is the one the repository
already states: **quote the gap over mmCSM-PPI, not the digits.** The gap holds
across every rule here — +0.146 to +0.166 Pearson, every interval clear of zero.

## Files

`mmcsm-ppi/data/selection_study/` holds each run's pre-registration, its closed
selection and its results; `*_PLAN.md` was written before the corresponding fit.
`grouped_inner_scores.npz` is the full 10 × 5 × 13 inner surface.
`mmcsm-ppi/data/penalty_sensitivity.json` is the blind surface, a diagnostic that
selects nothing. The runs themselves are
`/data/work/pythia_ppi/experiments/mmcsm_ridge_{repeated,matched,mixture}_20261115`.
