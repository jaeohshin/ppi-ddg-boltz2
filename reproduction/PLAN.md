# Does the blind performance survive my own Boltz features?

Written before the fits.

One question: replace the PreFold-dG authors' `skempi_{mut,wt}.pkl` with a local
re-fold (`scaled_v1_recycle3`: recycling 3, wild-type PDB as template, seed 42)
and refit the mmCSM-PPI blind comparison. Does the score stay the same?

Fixed in advance:

- One local variant, `scaled_v1_recycle3`, chosen because it is the setting
  recorded in `docs/PROVENANCE.md` as this project's documented configuration.
  Not chosen by score, and no other variant will be substituted afterwards.
- Rows the local re-fold does not cover are dropped from **both** arms, so the
  two numbers are computed on identical rows.
- Everything else is byte-identical to `mmcsm_ridge_blind_20261115/run.py`:
  same cohort file, same alpha grid, same five block weights, same 3-fold
  grouped inner CV, same selection on forward pooled inner OOF Spearman.
  The penalty is selected separately within each arm, as it would be if the
  local features were all anyone had.
- Reported either way. A lower local score is the answer to the question, not a
  reason to try another variant.

Reference numbers to beat/match (authors' features, full 586 test rows):
ridge Pearson 0.8565, mmCSM-PPI 0.7019 on the same rows.
