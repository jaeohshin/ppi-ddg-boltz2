"""Same cohort, same procedure, two feature producers.

Arm A: the PreFold-dG authors' released skempi_{mut,wt}.pkl, as every published
fit in this repository uses.
Arm B: a local Boltz re-fold (scaled_v1_recycle3, recycling 3, wild-type PDB as
template, seed 42).

Rows missing from either producer are dropped from BOTH arms, so the two scores
are computed on identical rows.  Nothing else differs: same cohort file, same
penalty grid, same inner grouped folds, same selection metric.
"""
from paths import (EXP, AUTHOR, LOCAL, VARIANT, BASE, PREP, DELTA, SKEMPI,
                   BLIND_COHORT, WIDTHS, D, OFFSETS, N_PREPARED_ROWS,
                   join_map, out_dir, missing)
VAR = VARIANT
import os, sys, json, time, csv
os.environ.update(OPENBLAS_NUM_THREADS='8', OMP_NUM_THREADS='8', MKL_NUM_THREADS='8')
from pathlib import Path

import numpy as np, pandas as pd, joblib
from scipy.stats import spearmanr, pearsonr, kendalltau
from sklearn.model_selection import GroupKFold

OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else out_dir()
OUT.mkdir(parents=True, exist_ok=True)
if missing():
    sys.exit('missing inputs:\n  ' + '\n  '.join(missing()))

WIDTHS = [147456, 147456, 128]
D = sum(WIDTHS)
T0 = time.time()
def log(m): print(f"[{time.time()-T0:6.0f}s] {m}", flush=True)

# ---- cohort -----------------------------------------------------------------
c = pd.read_csv(BLIND_COHORT)
prep = list(csv.DictReader(open(PREP)))
mut_id = np.array([r['mut_id'] for r in prep])
wt_id  = np.array([r['wt_id']  for r in prep])
log(f"mmCSM cohort {len(c)} rows; prepared cohort {len(prep)} feature rows")

def fv(sample):
    parts = []
    for k, shape in [('wsi', (384, 384)), ('ws', (384, 384)), ('wz', (128,))]:
        v = sample[k]
        if hasattr(v, 'detach'): v = v.detach().cpu().float().numpy()
        v = np.asarray(v, dtype=np.float32)
        if tuple(v.shape) != shape or not np.isfinite(v).all():
            raise ValueError('invalid features ' + sample['subpath'] + ' ' + k)
        parts.append(v.reshape(-1))
    return np.concatenate(parts)

# ---- which feature rows does the local re-fold cover? ------------------------
jm = json.load(open(join_map()))
MUTMAP, WTMAP = jm['mut'], jm['wt']
log(f"join map: {len(MUTMAP)} mutant records, {len(WTMAP)} wild-type records")

log("loading local wt pkl")
lws, _ = joblib.load(LOCAL / f'{VAR}_wt.pkl')
local_wt = {s['subpath']: fv(s) for s in lws}
del lws

need = np.unique(c.feature_row.to_numpy())
log("loading local mut pkl")
lms, _ = joblib.load(LOCAL / f'{VAR}_mut.pkl')
local_mut_ix = {s['subpath']: i for i, s in enumerate(lms)}
log(f"local: {len(local_wt)} wt, {len(lms)} mut records")

covered = np.array([r for r in need
                    if mut_id[r] in MUTMAP and MUTMAP[mut_id[r]] in local_mut_ix
                    and wt_id[r] in WTMAP and WTMAP[wt_id[r]] in local_wt])
dropped = sorted(set(need) - set(covered.tolist()))
log(f"feature rows needed {len(need)}; covered by local re-fold {len(covered)}; dropped {len(dropped)}")

keep = c.feature_row.isin(set(covered.tolist())).to_numpy()
c = c[keep].reset_index(drop=True)
log(f"cohort after dropping uncovered rows: {len(c)} "
    f"(train {(c['set']=='train').sum()}, test {(c['set']=='test').sum()})")

# ---- compact feature matrices, identical row order for both arms -------------
rows = np.unique(c.feature_row.to_numpy())
compact = {int(r): i for i, r in enumerate(rows)}
fi = np.array([compact[int(r)] for r in c.feature_row])
n = len(rows)
log(f"building two {n} x {D} float32 matrices ({n*D*4/2**30:.1f} GiB each)")

A = np.empty((n, D), dtype='float32')
src = np.memmap(DELTA, mode='r',
                dtype='float32', shape=(5729, D))
for i, r in enumerate(rows): A[i] = src[int(r)]
del src
log("arm A (authors) built from the existing delta_features.f32")

B = np.empty((n, D), dtype='float32')
for i, r in enumerate(rows):
    r = int(r)
    B[i] = fv(lms[local_mut_ix[MUTMAP[mut_id[r]]]]) - local_wt[WTMAP[wt_id[r]]]
del lms, local_mut_ix, local_wt
log("arm B (local recycle3) built")

cos = np.array([float(np.dot(A[i], B[i]) /
                      (np.linalg.norm(A[i]) * np.linalg.norm(B[i]) + 1e-12))
                for i in range(n)])
log(f"per-row cosine A vs B: median {np.median(cos):.4f}  "
    f"p05 {np.percentile(cos,5):.4f}  p95 {np.percentile(cos,95):.4f}")

# ---- the fit, verbatim from mmcsm_ridge_blind_20261115/run.py ----------------
tr = np.flatnonzero(c['set'] == 'train'); te = np.flatnonzero(c['set'] == 'test')
y = c.y.to_numpy(); sg = c.sign.to_numpy()
assert not set(fi[tr]) & set(fi[te])
alphas = np.logspace(-3, 3, 13); weights = [.25, .5, 1, 2, 4]

def kernels(x, a, b):
    ks, vs, off = [], [], 0
    for width in WIDTHS:
        k = np.zeros((len(a), len(a))); v = np.zeros((len(b), len(a)))
        for l in range(off, off + width, 2048):
            r = min(l + 2048, off + width)
            z = np.asarray(x[fi[a], l:r], dtype='float64') * sg[a, None]
            mu = z.mean(0); sd = np.maximum(z.std(0), 1e-8); z = (z - mu) / sd
            t = (np.asarray(x[fi[b], l:r], dtype='float64') * sg[b, None] - mu) / sd
            k += z @ z.T / width; v += t @ z.T / width
        ks.append(k); vs.append(v); off += width
    return ks, vs

def fit(x, label):
    ip = np.full((5, 13, len(c)), np.nan)
    for fold, (a, b) in enumerate(GroupKFold(3).split(tr, groups=c.iloc[tr].pdb)):
        a = tr[a]; b = tr[b]; b = b[sg[b] == 1]
        assert not set(c.iloc[a].pdb) & set(c.iloc[b].pdb)
        ks, vs = kernels(x, a, b)
        for wi, w in enumerate(weights):
            k = w * (ks[0] + ks[1]) + ks[2]; v = w * (vs[0] + vs[1]) + vs[2]
            ev, u = np.linalg.eigh(k); ev = np.maximum(ev, 0)
            ym = y[a].mean(); rhs = u.T @ (y[a] - ym)
            for ai, alpha in enumerate(alphas):
                ip[wi, ai, b] = v @ (u @ (rhs / (ev + alpha))) + ym
        log(f"  {label}: inner fold {fold} done")
    fwd = tr[sg[tr] == 1]
    sc = np.array([[spearmanr(y[fwd], p[fwd]).statistic for p in arm] for arm in ip])
    assert np.isfinite(sc).all()
    wi, ai = np.unravel_index(np.argmax(sc), sc.shape)
    w = weights[wi]; alpha = float(alphas[ai])
    ks, vs = kernels(x, tr, te)
    k = w * (ks[0] + ks[1]) + ks[2]; v = w * (vs[0] + vs[1]) + vs[2]
    ev, u = np.linalg.eigh(k); ev = np.maximum(ev, 0); ym = y[tr].mean()
    p = v @ (u @ ((u.T @ (y[tr] - ym)) / (ev + alpha))) + ym
    log(f"  {label}: selected w={w} alpha={alpha:.4g}")
    return dict(matrix_weight=w, alpha=alpha), p

def metrics(a, p):
    return dict(n=len(a), spearman=float(spearmanr(a, p).statistic),
                pearson=float(pearsonr(a, p).statistic),
                kendall=float(kendalltau(a, p).statistic),
                rmse=float(np.sqrt(np.mean((a - p) ** 2))),
                mae=float(np.mean(abs(a - p))))

selA, pA = fit(A, 'authors'); del A
selB, pB = fit(B, 'local');  del B
yt = y[te]; paper = c.iloc[te].paper.to_numpy()

# paired bootstrap over complexes: local minus authors, and each minus mmCSM
pdbs = c.iloc[te].pdb.to_numpy(); groups = [np.flatnonzero(pdbs == g) for g in pd.unique(pdbs)]
rng = np.random.default_rng(20261115); d_la, d_am, d_lm = [], [], []
for _ in range(4000):
    ix = np.concatenate([groups[j] for j in rng.integers(len(groups), size=len(groups))])
    if len(np.unique(yt[ix])) < 3: continue
    ra = pearsonr(yt[ix], pA[ix]).statistic
    rb = pearsonr(yt[ix], pB[ix]).statistic
    rm = pearsonr(yt[ix], paper[ix]).statistic
    d_la.append(rb - ra); d_am.append(ra - rm); d_lm.append(rb - rm)
def ci(d, point):
    if not d: return [point, None, None]
    return [point, float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))]
rA = pearsonr(yt, pA).statistic; rB = pearsonr(yt, pB).statistic; rM = pearsonr(yt, paper).statistic

result = {
 'question': 'does the published blind performance survive replacing the authors '
             'features with a local Boltz re-fold (scaled_v1_recycle3)?',
 'rows': {'cohort_after_coverage_filter': int(len(c)),
          'train': int((c['set'] == 'train').sum()), 'test': int(len(te)),
          'feature_rows_needed': int(len(need)), 'feature_rows_covered': int(len(covered)),
          'feature_rows_dropped': len(dropped),
          'dropped_wt_ids': sorted({wt_id[r] for r in dropped}),
          'test_complexes': int(len(groups))},
 'feature_agreement': {'per_row_cosine_median': float(np.median(cos)),
                       'p05': float(np.percentile(cos, 5)),
                       'p95': float(np.percentile(cos, 95))},
 'authors_features': {'selected': selA, **metrics(yt, pA)},
 'local_features':   {'selected': selB, **metrics(yt, pB)},
 'mmCSM_same_rows':  metrics(yt, paper),
 'paired_pearson_deltas_over_complexes': {
     'local_minus_authors': ci(d_la, float(rB - rA)),
     'authors_minus_mmCSM': ci(d_am, float(rA - rM)),
     'local_minus_mmCSM':   ci(d_lm, float(rB - rM))},
 'elapsed_seconds': time.time() - T0}
(OUT / 'RESULTS.json').write_text(json.dumps(result, indent=2))
pd.DataFrame({'pdb': pdbs, 'mutation': c.iloc[te].mutation.to_numpy(), 'y': yt,
              'authors': pA, 'local': pB, 'mmCSM': paper}).to_csv(OUT / 'predictions.csv', index=False)
print(json.dumps(result, indent=2))
