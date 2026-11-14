# Results

The record behind every number in the README. Source experiments are in the working
repository `jaeohshin/prefold-dg-mechanism`, whose ledger (`docs/EXPERIMENTS.md`) carries
the commit, inputs, command, job id and artifact path for each; the date in parentheses
is the ledger entry.

Metric throughout: **per-structure PCC**, a Pearson correlation inside each complex
averaged over the complexes with at least ten variants (116 of 334 on the full cohort).
Labels are differences of log10(Kd in nM), not kcal/mol — multiply by 1.364. Intervals
are paired bootstraps over complexes.

---

## 1. How the question changed

The project asked three questions in sequence; each answer produced the next question.

| stage | period | question | outcome |
|---|---|---|---|
| 1 | 2026-09-02 – 09-08 | at which Boltz-2 layer and recycle does ΔΔG information appear? | on a 100-pair pilot, mutation-centred representations did not beat the mean baseline: MAE 1.508 against 1.500 kcal/mol, −0.008 [−0.114, +0.095]. Extended to the full 5,735-variant cohort |
| 2 | 09-12 – 10-07 | why do our own extracted features train 0.12 worse than the authors' released ones? | the wild-type template. Confirmed by the corresponding author |
| 2b | 09-20 – 10-06 | does changing the head, or combining with Pythia-PPI, help? | every route's interval covers zero |
| 3 | 10-07 – 10-08 | what is the head reading out of its input? | a linear reader matches it, and beats it once homology is controlled |

## 2. The reproduction gap, and what closed it

Running the published `main.py` on the authors' released pickles gives per-structure PCC
0.590 against the 0.611 in print. Running it on features extracted by this project's own
Boltz-2 pipeline gave 0.492. That 0.12 was traced over three weeks:

| date | test | result | reading |
|---|---|---|---|
| 09-13 | recycle 3, own trace | 0.492 against the authors' 0.590 | recycle count suspected |
| 09-14 | re-extract at recycle 10 | 0.468, worse | recycle count is not the cause |
| 09-17 | decompose the gap from existing artifacts | 131 of 328 wild-type complexes pool over the wrong interface; ipTM correlation with the authors only 0.485 | structure quality |
| 09-18 | freeze the authors' head, replace the features of the 12 worst complexes | 0.159 → 0.451 (authors 0.540) | 77% of the gap recovered; full re-extraction submitted |
| 09-19 | full template re-extraction, recycle 10 | mutant-minus-wild-type cosine on `wz` 0.577 → 0.953; ipTM 0.911 against the authors' 0.914; training 0.5289 ± 0.009 | half the gap closed |
| 10-06 | labels | wild type identical 5,471/5,471, mutant 99.3% | labels are not the cause |
| 10-07 | re-run our own pipeline against itself | self-reproduction cosine 1.00000 | no global setting remains; the residual is per-variant |
| 10-07 (f) | **the corresponding author replied** | both mutant and wild type were folded with the wild-type PDB as a template; checkpoint SHA-1 matches; MSA greedy pairing; sequence mismatches fixed by alignment | **closed** |

The authors' Boltz runs saw the deposited structure and ours did not. With templates, our
`wz` mutation response agrees with theirs at cosine 0.95 — better than two of our own
template-free runs agree with each other (0.92).

One difference remains untested: the authors corrected RCSB/SKEMPI sequence mismatches by
alignment and we did not.

## 3. The head reads something linear

### 3.1 Which blocks carry the signal

Branch ablation on a frozen head, replacing one branch by its mean (09-20):

| branch removed | cost |
|---|---|
| `wz` (128 numbers) | −0.082 / −0.087 |
| `ws` (147,456) | −0.046 |
| `wsi` (147,456) | −0.025 |

58.3M of the head's 58.4M parameters sit in the two flattening layers that read `wsi` and
`ws`. The 128 numbers of `wz` carry most of the signal.

### 3.2 The block penalty

Reading the three blocks under one shared ridge penalty scores about 0.53; the two large
matrices dominate the small one. Giving each block its own penalty — dividing each block's
Gram by its width — gives 0.59 to 0.62 (10-07 b, c).

| comparison | reader | head | difference |
|---|---|---|---|
| our template features, `rdenet` | 0.589 | 0.521 | +0.069 [+0.039, +0.099] |
| authors' features, `rdenet` | 0.607 | 0.580 | +0.027 [−0.009, +0.059] |
| full cohort, 5,735 | | | +0.051 [+0.023, +0.078] |
| Stage A split, variants shared with Pythia-PPI | 0.553 | 0.513 (Pythia 0.441) | +0.040 [+0.010, +0.072] |

A fair summary of these is **+0.04 to +0.05**. Six checks passed, including an independent
primal solver reproducing the kernel path to within 0.001, selection by the head's own
criterion (inner MAE) rather than the reported metric, and the full cohort with the
published loader's pairing. The reader is never below the head under any condition tested.

Reading the authors' features and ours with the same reader differs by −0.017 [−0.045,
+0.010] — indistinguishable. **The residual gap was the head, not the features.**

The GPU head-tuning grid was deliberately not run: no finite grid can establish that no
training procedure reaches the reader.

An earlier explanation of the block penalty — that it protects weak blocks from strong
ones — was **withdrawn** (10-08 b). A control adding a deliberately weak block showed the
0.049 is the cost of the single shared penalty itself.

### 3.3 Clean-room reproduction and the external test

From the released pickles alone (10-07 d, e):

| | here | printed |
|---|---|---|
| SKEMPI per-structure PCC | 0.622 | 0.611 |
| SKEMPI per-structure SRCC | 0.597 | 0.573 |
| SKEMPI overall PCC | 0.744 | 0.750 |
| SKEMPI AUROC | 0.798 | 0.790 |
| HER2, 419 variants, PCC | 0.546 [0.483, 0.606] | 0.524 |
| HER2, SRCC | 0.542 | 0.529 |

Higher on all five HER2 folds, but the published value is inside the interval: a match,
not a win. Two asymmetries favour the reader and should be stated — it never sees HER2,
while the published protocol selects a stopping epoch on it.

### 3.4 The absolute ΔG task comes free

The paper's Table 4 task — predict mutant and wild-type ΔG separately and subtract — is
also linear (10-08 a): overall PCC 0.723 / SRCC 0.725 against the paper's 0.711 / 0.712.
RMSE is worse (2.26 against 1.72 kcal/mol; a ranking penalty compresses between-complex
scale). The decomposition costs −0.001, because for a linear reader
`w·x_mut − w·x_wt = w·(x_mut − x_wt)` — no shared neural module is needed to do it.

## 4. Homology

### 4.1 The antibody subset, and why the split had to change

Antibody–antigen variants (933, 50 PDB ids, 25 scored) score 0.47–0.48 against 0.64–0.66
for the rest, a difference of +0.17 to +0.18 [+0.04, +0.31] (10-08 d). Not a composition
artefact: within-complex ΔΔG spread 0.476 against 0.488, 27 variants per complex against
24. This is substantially the same data as Hummer et al.'s `Experimental_ddG_608`.

Computing chain identity between those complexes from the Boltz input sequences (10-08 e):
at 99% identity the **50 PDB ids collapse to 22 clusters** and the 25 scored ones to 13,
with a single cluster of 14 hen-lysozyme complexes. A split by PDB id puts the same
antibody on both sides.

### 4.2 Refitting under control

611 distinct wild-type chain sequences over the 333 cohort PDB ids were clustered by
single linkage; two PDB ids are linked when they share a chain cluster (10-08 f, g, h).

| split | reader | head | reader − head |
|---|---|---|---|
| `rdenet` (PDB id) | 0.604 / 0.626 | 0.589 | +0.015 [−0.023, +0.055] |
| `hom99` (127 groups) | 0.573 | 0.505 | +0.068 [+0.030, +0.109] |
| `holdoutprot` (159 groups) | 0.579 | 0.524 | +0.055 [+0.021, +0.088] |

The reader falls −0.053 [−0.086, −0.021] from `rdenet` to `hom99`; the head falls −0.086,
and −0.162 on the antibody subset. **Both fall**, so the leakage is about split
strictness in general rather than antibodies specifically — and part of the headline 0.62
is homology recall. But the head falls twice as far.

Note against a "more room to fall" account: the reader starts higher (0.626 against 0.589)
*and* falls less.

The two `rdenet` reader values are the open penalty-selection problem of §7.1 in the
headline table: the `table2` out-of-fold run scores 0.6040 and the report run 0.6216, and
the +0.0147 difference quoted against the head's 0.5893 is measured on the first of them
(2026-10-08 d). The `hom99` and `holdoutprot` rows carry no such ambiguity — there the
difference is the subtraction of the two columns. This is the clearest reason the project
settled on quoting the gap rather than a reader score.

`holdoutprot` is the SKEMPI curators' own grouping, is what PPIformer split on, and is the
one Bushuiev et al. measured at 0% leakage. It is also the looser of the two controlled
splits (159 groups against 127), so the gap is not an artefact of over-strict grouping.
Two splits built on unrelated principles agree to within 0.013.

## 5. What did not work

| approach | outcome | date |
|---|---|---|
| mutation-centred representations, 100-pair pilot | MAE −0.008 [−0.114, +0.095] against a mean baseline | 09-08 |
| other Pairformer layers and recycles (88 snapshots) | flat after recycle 3–4; best +0.005 [−0.019, +0.031] | 09-20 |
| site-local features | pcc/MAE signal present, per-structure interval covers zero | 09-20 |
| an MLP on the same input | −0.002 [−0.021, +0.017] on the full data | 09-21 |
| distillation from teacher soft labels | +0.011 [+0.004, +0.018]; a regularisation effect | 09-20 |
| combining with Pythia-PPI | residual correlation 0.765; a per-variant gate is not learnable (0.47 against a 0.50 base rate) | 09-30, 10-06 |
| richer inputs from the saved traces | ≤ +0.01 over a reader that reads all three tensors | 10-07 a |
| ESM-2 650M as an extra block | 0.596 against 0.615; not significant | 10-07 i |
| Pythia-PPI predictions as a block | 0.587 against 0.593; averaging outputs 0.591 | 10-08 b |
| head capacity, 0.1M to 58.4M | all below the reader; the converse is not provable on a finite grid | 10-07 c |

Half the training complexes to all of them is worth +0.029 [+0.013, +0.047] to the reader
and +0.040 [+0.021, +0.058] to an MLP (09-21), so coverage contributes — but it is not the
only bottleneck, and Hummer et al. put the requirement for a test Pearson above 0.85 at
≥90,000 measured values against SKEMPI's 5,731.

## 6. Withdrawn and corrected

Kept because the corrections are part of the record.

| original claim | correction |
|---|---|
| templates do not improve training (09-15) | per-structure PCC on 600 variants averages six complexes — not measurable (09-16). On the full cohort, +0.06 (09-19) |
| the authors' Boltz settings are unpublished (09-16) | they are in Supplementary Table S2 (10-06) — though those settings alone do not reproduce the features |
| a "sample-size effect" explains the gap (09-15) | a metric artefact; withdrawn (09-16) |
| the pooled input is the limit, not the head (09-20) | half right. Head capacity is not the answer, but the cause was wrong: a single-penalty ridge was hiding the head's loss (10-07) |
| block penalties protect weak blocks from strong ones (10-07 i) | it is the cost of the single penalty itself; withdrawn (10-08 b) |
| reader and head are "comparable" (10-07) | only on `rdenet`; under homology control the reader leads significantly (10-08 g, h) |

## 7. Open

1. **Penalty selection.** Three nominally identical runs chose three different penalty
   sets, landing between 0.604 and 0.622. The inner surface is flat (0.006 across a
   300-fold range) while the outer consequence is 0.064, and nested CV's smaller inner
   folds bias toward over-regularisation. float32 and float64 behave identically. The
   conclusion is unaffected — the gap is positive from alpha 1 to 1000 — but a single
   reader score must not be quoted.
2. The paper's RMSE and MAE are higher than a retrained head's, and the retrained head
   scores 0.028 below the printed PCC. Unexplained.
3. Whether antibody–antigen is intrinsically harder or only more homologous. Under control
   both subsets fall, so this is unsettled. Testable cheaply: correlate per-structure
   performance against MSA depth and Boltz-2 confidence, on the hypothesis that somatically
   hypermutated CDRs carry little evolutionary covariation.
4. The effect of the sequence-mismatch alignment correction the authors applied.
5. Output form. A Letter addressed to one paper is neither planned nor excluded; the
   argument against has strengthened, since the authors were careful by the standards of
   this literature — they excluded nine HER2-homologous entries by Smith-Waterman before
   their external test, which others using the same test set do not — while a survey found
   that **no paper in this literature reports a linear or trivial baseline as a control**.
   That makes the finding a property of the field rather than of one paper.
