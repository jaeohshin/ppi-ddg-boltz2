# Leave-one-family-out, pre-registered 2026-11-15

Nothing in this cohort has been scored. It is the arbiter for the selection
question docs/SELECTION.md leaves open, and it is the first test here of whether
the margin over mmCSM-PPI survives a similarity-controlled split.

## Design, and the asymmetry in it

mmCSM-PPI's predictions exist only for the 595 blind rows; the training file's
prediction column is not identified as out-of-fold and cannot stand in for one.
So the test side has to sit inside the blind set.

Families are the 99%-chain-identity single-linkage clusters already in
`mmcsm-ppi/data/families_0.99.json`: 79 families over the 147 complexes, 73 on
the training side, 36 on the blind side, 30 on both.

For each of the 36 blind families f: fit on every training row whose family is
not f, predict the blind rows whose family is f. Pool the 36 held-out blocks
into one prediction per blind row and compare against mmCSM-PPI's published
value on the same rows.

**The asymmetry must be stated with every number this produces.** mmCSM-PPI's
published blind predictions come from a model trained on the full training set,
which contains family f. Our model is denied it. The comparison therefore
handicaps us, by exactly the thing being measured. If the ridge still wins, that
is a strong result. If it loses, the result is confounded and says nothing about
which method generalises better -- only that one of them was allowed to see the
family. No number from this run may be reported without that sentence.

## The selection question this settles

docs/SELECTION.md ends with a principle -- choose the penalty by an inner split
that matches the deployment condition -- and an untested rule. The principle
makes a falsifiable prediction here, because the deployment condition has
changed. The blind set was 88% seen complexes, so matching it meant a mostly
ungrouped inner split. A family-disjoint test is 100% unseen families, so the
principle now demands a family-grouped inner split, which is the rule that
scored WORST on the blind set.

Two rules, both nested properly inside each outer fold, penalty chosen without
any sight of family f:

- **A, principle-matched:** inner folds grouped by family, 10 repeats x 3 folds,
  argmax of mean inner out-of-fold Spearman.
- **B, blind-set-matched:** the mixture rule from docs/SELECTION.md, 30
  Monte-Carlo splits at 11.95% held-out-complex rows, unchanged.

**Prediction, fixed now: A beats B on the pooled family-disjoint score.** If B
wins or they tie, the principle is wrong as stated and docs/SELECTION.md gets
corrected rather than defended.

A third arm, **C**, applies the fixed penalty (weight 2, alpha 1) that the blind
surface showed to be optimal under the lenient split. C is a diagnostic, not a
rule: it is the question "does the lenient optimum transfer", and it is reported
whatever it says.

## Everything else unchanged

Same features, same 1,237-row cohort, same grid, same standardisation on
inner-training rows only, reverse rows on training sides only, validation and
scoring on forward rows only. Paired bootstrap over blind families rather than
over complexes, since families are the unit of independence here.
