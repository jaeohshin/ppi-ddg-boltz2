# Mixture-matched inner selection, pre-registered 2026-11-15

## Honesty note, written first and not to be removed
By this point the blind surface has been computed and looked at as a diagnostic
(data/penalty_sensitivity.json in the published repo). I know that the blind
optimum is near weight 2, alpha 1 at Pearson 0.8748. **Nothing from this run can
be reported as blind confirmation.** What it can honestly establish is whether a
selection rule justified on a stated principle, with its target set from the
public composition of the released files rather than from their labels, lands
nearer the optimum than the two schemes already tried. That is a weaker claim
and it is the only one this run is allowed to make.

Three runs have now scored this blind set. That is three looks, and the right
way to settle a selection rule is a cohort nobody has scored yet -- the
family-level split, which is the next item and where no penalty has been chosen.

## What the two previous runs showed
Blind Pearson at each rule's selected penalty, complete 1,237-row cohort:

| inner split | selected | inner score | blind Pearson |
|---|---|---|---|
| grouped by PDB, 10x3 | w=0.5, alpha=3.16 | 0.702 | 0.8565 |
| row-level, 10x3 | w=4, alpha=0.1 | 0.885 | 0.8424 |
| (grid optimum) | w=2, alpha=1 | -- | 0.8748 |

The grouped inner estimate is 0.119 pessimistic against the blind Spearman of
0.821; the row-level one is 0.064 optimistic. They miss the penalty in opposite
directions. The blind set's own condition is a mixture: 516 of 586 rows sit on
complexes present in the training release, 70 do not.

## The rule, fixed before any fit
Monte-Carlo inner validation, 30 random splits, matching that mixture.

Each split: shuffle complexes; set aside complexes until their rows reach
11.95% of the target validation size (70/586, the blind set's own proportion)
and hold those complexes out of training entirely; fill the rest of the
validation set with rows drawn from complexes that remain in training. Score
pooled Spearman over the forward validation rows. **Selected penalty: the argmax
of the mean over the 30 splits.** Same cohort, features, grid, standardisation
and reverse-row handling as the other two runs. 30 fits, so the compute matches.

The 11.95% comes from counting rows in the released files, not from any label or
score.

## Prediction
The selected alpha lands strictly between 0.1 and 3.16, and the selected weight
between 0.5 and 4. Whether the resulting blind score beats 0.8565 is recorded
but is not evidence for the rule, for the reason in the honesty note.
