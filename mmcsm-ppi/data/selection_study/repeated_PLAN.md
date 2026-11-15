# Repeated-fold penalty selection, pre-registered 2026-11-15

## The problem this addresses, measured before this run
On the 1,237-row mmCSM cohort, six extra training rows (0.5% of the data) moved
the inner-OOF argmax from (weight 0.5, alpha 1) to (weight 1, alpha 3.16) and
the blind Pearson from 0.8697 to 0.8648. Holding the penalty fixed, those rows
were worth +0.0014 at either setting, so the entire move was selection, not
data: -0.0049 against +0.0014. The single 3-fold GroupKFold split used by
mmcsm_ridge_blind_20261106/run.py decides the penalty on one partition of 73
PDB groups, and that is too little evidence to pin it.

## The rule, fixed before any fit
Features, cohort, labels, reverse-row handling, candidate grid and the
"validation scored on forward rows only" convention are all unchanged from
mmcsm_ridge_blind_20261115. The only change is how the penalty is chosen.

- 10 repeats. Repeat r shuffles the 73 training PDB groups with seed r and
  assigns them greedily to 3 folds, smallest fold first by row count, so folds
  are balanced in rows and no group is split.
- For each repeat, pooled out-of-fold Spearman over the forward training rows,
  for every (block weight, alpha) in the existing grid: weight in
  {0.25, 0.5, 1, 2, 4}, alpha in logspace(-3, 3, 13).
- **The reported penalty is the argmax of the mean over the 10 repeats.** This
  is the rule; it is chosen now and not revisited.
- Two further quantities are computed as diagnostics and will be reported
  whatever they say: the per-repeat argmax (how often the single-split rule
  would have landed somewhere else) and the one-standard-error choice (the
  largest alpha within 1 SE of the best mean). Neither selects the reported
  model.
- Refit once on all 1,237 training rows at the selected penalty. The blind set
  is scored exactly once, after selection is closed.

## What would count as success
Not a higher blind score -- a selection that does not move. Specifically: the
repeat-averaged surface should have a single clear argmax, and the per-repeat
argmaxes should concentrate on it. If the per-repeat argmaxes are still spread
across the grid, the instability is in the data, not in the splitting, and the
honest conclusion is that this cohort cannot pin the penalty at all. That
outcome gets reported as such.

Standardisation is fit on inner-training rows only, as before. Reverse rows
enter the inner training side only and never the validation side.
