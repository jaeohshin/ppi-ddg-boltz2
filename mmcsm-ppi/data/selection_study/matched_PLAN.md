# Matching the inner split to the test condition, pre-registered 2026-11-15

## The observation this tests
The repeated-fold run (mmcsm_ridge_repeated_20261115) stabilised the penalty
argmax and produced a WORSE blind score: 0.8565 against 0.8648 and 0.8697 from
two single-split runs, with the grid's best cell at 0.8748. The whole inner
surface, averaged over 10 repeats, peaks at alpha 3.16-10 for every block
weight. The blind surface peaks at alpha 1. That is a systematic offset, not
noise, and more repeats cannot remove it.

Candidate explanation: the inner folds are grouped by PDB, so inner validation
rows sit on complexes the inner training set never saw. The blind set is not
like that -- 516 of its 586 rows share a PDB entry with the training release,
because mmCSM-PPI split by mutation count, not by similarity. Grouped inner CV
therefore estimates new-complex generalisation while the blind set measures
mostly seen-complex accuracy. Harder task, more regularisation preferred.

## Prediction, fixed before the fit
An inner split that matches the test condition -- rows assigned at random,
complexes free to appear on both sides, which is exactly the blind set's
relationship to training -- will select a SMALLER alpha, nearer 1, and the model
it selects will score HIGHER on the blind set than the grouped rule's 0.8565.

If instead the matched split also lands at alpha 3.16 or above, the explanation
is wrong and the offset is something else; that gets reported as a refutation.

## Protocol
Identical to the repeated run -- same cohort, features, grid, 10 repeats,
3 folds, standardisation on inner-training rows only, reverse rows on the
training side only, validation scored on forward rows only -- with one change:
folds are formed by shuffling ROWS, not PDB groups.

Both rules are reported side by side. The blind set is scored once per rule,
after that rule's selection is closed. Neither selection looks at it.

## Interpretation fixed in advance
A matched split is the right estimator for THIS benchmark's condition, and
saying so is not an endorsement of the condition. The grouped rule remains the
right estimator for the question "will this work on a new complex", and the
homology grading already answers that one unfavourably. If the matched split
wins here, what it demonstrates is that the penalty must be chosen for the
deployment condition -- and that this benchmark's condition is a lenient one.
