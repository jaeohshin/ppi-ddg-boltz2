# What survives a similarity-controlled split

The headline comparison in `mmcsm-ppi/` is fair but lenient: mmCSM-PPI separated
its blind set by mutation count, so 516 of its 586 scorable rows sit on complexes
the training release contains. This is the same comparison run with the
similarity controlled, and it changes two things — one about the benchmark, one
about the model.

## Design, and the asymmetry that qualifies every number here

mmCSM-PPI's predictions exist only for the blind rows, so the test side has to
stay inside the blind set. Families are the 99%-identity single-linkage clusters
in `mmcsm-ppi/data/families_0.99.json`. For each of the 34 families the cohort
covers: fit on every training row outside that family, predict that family's
blind rows, pool the 34 held-out blocks. The penalty is chosen inside each outer
fold, with the held-out family absent.

**mmCSM-PPI's published predictions come from a model trained on the full
training set, which contains the family being held out. Ours is denied it.** The
comparison handicaps the ridge by exactly the quantity under test. That sentence
belongs with every number below. If the training file's
`prediction_mmCSM-PPI` column were identified as out-of-fold the asymmetry could
be removed; it is not, so one question to the authors would settle this
experiment, and it is worth asking.

## Result 1: the margin over mmCSM-PPI does not survive

| | Pearson | ρ | RMSE | MAE |
|---|---|---|---|---|
| mmCSM-PPI *(saw every family)* | **0.702** | 0.650 | 2.035 | 1.489 |
| ridge, family-grouped selection | 0.555 | 0.620 | 2.380 | 1.800 |
| ridge, mixture-matched selection | 0.403 | 0.405 | 2.934 | 2.134 |
| ridge, fixed w=2 α=1 (the lenient optimum) | 0.451 | 0.460 | 2.695 | 1.982 |

Paired over families, best arm against mmCSM-PPI: Pearson −0.076
[−0.343, +0.370], ρ −0.010 [−0.193, +0.299]. Both intervals span zero widely.

**This is neither a win nor a loss — it is undetermined, and undetermined in the
ridge's disfavour by design.** On rank correlation the handicapped ridge is level
with an mmCSM-PPI that saw the family, which is worth noting and is not worth
more than noting, given the interval.

What is no longer in doubt: the repository's +0.163 Pearson margin is a property
of the lenient split. Scored family-disjoint, the ridge falls from 0.865 to
0.555. The parent README already said the win was earned on complexes the model
had seen; this is the number for it.

## Result 2: the penalty was never unstable — the validation split was wrong

[`docs/SELECTION.md`](SELECTION.md) left a principle untested: choose the penalty
by an inner split that matches the deployment condition. A family-disjoint test
is 100% unseen families, so the principle demands family-grouped inner folds —
the rule that had scored *worst* on the lenient blind set. Pre-registered
prediction: it beats the mixture rule here.

It does, by 0.152 Pearson, and the stronger result is in how it chose:

| inner split | penalty chosen, across 34 outer folds |
|---|---|
| **family-grouped** | **w=0.25 in 34 of 34**; α=3.16 in 28, α=10 in 6 |
| mixture-matched | w=4 in 30, scattered over four cells |

The instability that started all of this — six training rows, 0.5% of the
cohort, moving the argmax and costing 3.5× what the rows gained — **was not in
the data.** Matched to the task, the selection is pinned: one weight for every
fold. The earlier wandering was an inner split measuring one task while the
penalty was being chosen for another.

The lenient optimum (w=2, α=1) transfers badly: 0.451, worse than the matched
rule by 0.10. Tuning on a lenient condition is not merely uninformative for a
strict one, it is harmful.

## Result 3: only 128 of the 295,040 features generalise

w=0.25 was the bottom of the grid, so the grid was extended downward and to zero
— zero meaning the 128-wide pair block alone, both 147,456-wide sequence blocks
dropped. Family-disjoint scores (a diagnostic over the grid, selecting nothing;
`mmcsm-ppi/data/family_split/lofo_weight_surface.npz`):

| sequence-block weight | Pearson | ρ |
|---|---|---|
| **0 — pair block alone, 128 features** | **0.5579** | 0.6229 |
| 0.03125 | 0.5590 | 0.6279 |
| 0.125 | 0.5545 | **0.6341** |
| 0.25 *(what the rule chose)* | 0.5454 | 0.6322 |
| 1 | 0.5124 | 0.5850 |
| 4 *(what the lenient rules chose)* | 0.4487 | 0.4753 |

**Dropping both 147,456-wide blocks costs 0.0011 Pearson and 0.0112 ρ**, and the
score falls monotonically as their weight rises above about 0.1.

So the 294,912 sequence dimensions contribute nothing measurable once the family
is new. On the lenient split they were worth a great deal and the selection rules
pushed their weight to 4 — which is what memorising a complex looks like from the
inside. The honest model for a new protein family is a ridge on **128** frozen
pair-block features, and it is within 0.001 of the best this feature set can do
under that condition.

That is a smaller claim than beating mmCSM-PPI, and a more defensible one.

## What this says to do next

- Refit and report the pair-block-only model as the primary model for
  generalisation claims. It is two orders of magnitude smaller and loses nothing
  that can be measured here.
- Any future benchmark must be family-disjoint before a margin is quoted from it.
- Ask the mmCSM-PPI authors whether the released training-set prediction column
  is out-of-fold. One answer removes the only asymmetry in Result 1.
- The 128 dimensions are now small enough to look at directly. Which of them
  carry the signal, and whether they correspond to anything interpretable about
  the interface, is a question this feature set can finally be asked.

Runs: `/data/work/pythia_ppi/experiments/mmcsm_family_lofo_20261115`.
Artefacts, including the pre-registration written before the fits:
`mmcsm-ppi/data/family_split/`.
