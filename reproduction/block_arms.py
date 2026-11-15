"""Per-block scores for both feature producers, on identical rows.

If the loss is caused by cancellation in the mut-minus-wt difference, it should
be worst in the block where the difference is the smallest fraction of the
vector - `wz`, the 128-wide pair block that is also the only block that
generalises.
"""
from paths import (EXP, AUTHOR, LOCAL, VARIANT, BASE, PREP, DELTA, SKEMPI,
                   BLIND_COHORT, WIDTHS, D, OFFSETS, N_PREPARED_ROWS,
                   join_map, out_dir, missing)
VAR = VARIANT
import os, json, csv, time
os.environ.update(OPENBLAS_NUM_THREADS='8', OMP_NUM_THREADS='8', MKL_NUM_THREADS='8')
from pathlib import Path

import numpy as np, pandas as pd, joblib
from scipy.stats import spearmanr, pearsonr
from sklearn.model_selection import GroupKFold

W = [147456, 147456, 128]; D = sum(W); NAMES = ['wsi', 'ws', 'wz']
T0 = time.time()
def log(m): print(f"[{time.time()-T0:6.0f}s] {m}", flush=True)

jm = json.load(open(join_map())); MUT, WT = jm['mut'], jm['wt']
prep = list(csv.DictReader(open(PREP)))
c = pd.read_csv(BLIND_COHORT)
def fv(s):
    out = []
    for k, shp in [('wsi', (384, 384)), ('ws', (384, 384)), ('wz', (128,))]:
        v = s[k]
        if hasattr(v, 'detach'): v = v.detach().cpu().float().numpy()
        out.append(np.asarray(v, dtype=np.float32).reshape(-1))
    return np.concatenate(out)
log("loading local pkls")
lws, _ = joblib.load(LOCAL / f'{VAR}_wt.pkl'); lw = {s['subpath']: fv(s) for s in lws}; del lws
lms, _ = joblib.load(LOCAL / f'{VAR}_mut.pkl'); lmi = {s['subpath']: i for i, s in enumerate(lms)}

need = np.unique(c.feature_row.to_numpy())
cov = [int(r) for r in need if prep[r]['mut_id'] in MUT and MUT[prep[r]['mut_id']] in lmi
       and prep[r]['wt_id'] in WT and WT[prep[r]['wt_id']] in lw]
c = c[c.feature_row.isin(set(cov))].reset_index(drop=True)
rows = np.unique(c.feature_row.to_numpy()); compact = {int(r): i for i, r in enumerate(rows)}
fi = np.array([compact[int(r)] for r in c.feature_row]); n = len(rows)
log(f"{len(c)} cohort rows, {n} feature rows")

A = np.empty((n, D), 'float32'); src = np.memmap(DELTA, mode='r', dtype='float32', shape=(5729, D))
for i, r in enumerate(rows): A[i] = src[int(r)]
del src
B = np.empty((n, D), 'float32')
for i, r in enumerate(rows):
    r = int(r); B[i] = fv(lms[lmi[MUT[prep[r]['mut_id']]]]) - lw[WT[prep[r]['wt_id']]]
del lms, lmi, lw
log("both matrices built")

tr = np.flatnonzero(c['set'] == 'train'); te = np.flatnonzero(c['set'] == 'test')
y = c.y.to_numpy(); sg = c.sign.to_numpy()
ALPHAS = np.logspace(-2, 7, 37)
OFF = [0, W[0], W[0]+W[1]]

def gram(x, a, b, bi):
    off, width = OFF[bi], W[bi]
    k = np.zeros((len(a), len(a))); v = np.zeros((len(b), len(a)))
    for l in range(off, off+width, 2048):
        r = min(l+2048, off+width)
        z = np.asarray(x[fi[a], l:r], 'float64') * sg[a, None]
        mu = z.mean(0); sd = np.maximum(z.std(0), 1e-8); z = (z-mu)/sd
        t = (np.asarray(x[fi[b], l:r], 'float64') * sg[b, None] - mu)/sd
        k += z@z.T/width; v += t@z.T/width
    return k, v

def arm(x, blocks):
    ip = np.full((len(ALPHAS), len(c)), np.nan)
    for a_, b_ in GroupKFold(3).split(tr, groups=c.iloc[tr].pdb):
        a = tr[a_]; b = tr[b_]; b = b[sg[b] == 1]
        k = np.zeros((len(a), len(a))); v = np.zeros((len(b), len(a)))
        for bi in blocks:
            kk, vv = gram(x, a, b, bi); k += kk; v += vv
        ev, u = np.linalg.eigh(k); ev = np.maximum(ev, 0); ym = y[a].mean(); rhs = u.T@(y[a]-ym)
        for ai, al in enumerate(ALPHAS): ip[ai, b] = v@(u@(rhs/(ev+al)))+ym
    fwd = tr[sg[tr] == 1]
    sc = [spearmanr(y[fwd], p[fwd]).statistic for p in ip]
    ai = int(np.argmax(sc)); al = float(ALPHAS[ai])
    k = np.zeros((len(tr), len(tr))); v = np.zeros((len(te), len(tr)))
    for bi in blocks:
        kk, vv = gram(x, tr, te, bi); k += kk; v += vv
    ev, u = np.linalg.eigh(k); ev = np.maximum(ev, 0); ym = y[tr].mean()
    p = v@(u@((u.T@(y[tr]-ym))/(ev+al)))+ym
    return p, al, ai in (0, len(ALPHAS)-1)

ARMS = {'all three': [0, 1, 2], 'wz only (128)': [2], 'wsi only': [0], 'ws only': [1], 'wsi+ws': [0, 1]}
yt = y[te]; paper = c.iloc[te].paper.to_numpy()
pdbs = c.iloc[te].pdb.to_numpy(); groups = [np.flatnonzero(pdbs == g) for g in pd.unique(pdbs)]
rng = np.random.default_rng(20261115)
boot = [np.concatenate([groups[j] for j in rng.integers(len(groups), size=len(groups))]) for _ in range(3000)]
boot = [ix for ix in boot if len(np.unique(yt[ix])) > 2]

out = {}
for name, blocks in ARMS.items():
    pA, alA, edgeA = arm(A, blocks); pB, alB, edgeB = arm(B, blocks)
    d = [pearsonr(yt[ix], pB[ix]).statistic - pearsonr(yt[ix], pA[ix]).statistic for ix in boot]
    rA, rB = pearsonr(yt, pA).statistic, pearsonr(yt, pB).statistic
    out[name] = {'authors': {'pearson': float(rA), 'spearman': float(spearmanr(yt, pA).statistic), 'alpha': alA, 'grid_edge': edgeA},
                 'local':   {'pearson': float(rB), 'spearman': float(spearmanr(yt, pB).statistic), 'alpha': alB, 'grid_edge': edgeB},
                 'local_minus_authors_pearson': [float(rB-rA), float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))]}
    log(f"{name:16s} authors r={rA:.4f}  local r={rB:.4f}  delta {rB-rA:+.4f} "
        f"[{np.percentile(d,2.5):+.4f},{np.percentile(d,97.5):+.4f}]")
out['mmCSM_same_rows'] = {'pearson': float(pearsonr(yt, paper).statistic), 'spearman': float(spearmanr(yt, paper).statistic)}
out['rows'] = {'train': int(len(tr)), 'test': int(len(te)), 'test_complexes': int(len(groups))}
json.dump(out, open(out_dir()/'BLOCK_ARMS.json', 'w'), indent=1)
