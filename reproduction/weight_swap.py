"""Is the loss in the features, or in the block weight the inner CV picked?

Score both producers across the whole (block weight, alpha) surface, test rows,
selecting nothing. If the two surfaces are similar and only the selected point
differs, the features are not the problem.
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

W = [147456, 147456, 128]; D = sum(W); OFF = [0, W[0], W[0]+W[1]]
jm = json.load(open(join_map())); MUT, WT = jm['mut'], jm['wt']
prep = list(csv.DictReader(open(PREP)))
c = pd.read_csv(BLIND_COHORT)
def fv(s):
    out = []
    for k in ('wsi', 'ws', 'wz'):
        v = s[k]
        if hasattr(v, 'detach'): v = v.detach().cpu().float().numpy()
        out.append(np.asarray(v, dtype=np.float32).reshape(-1))
    return np.concatenate(out)
lws, _ = joblib.load(LOCAL/f'{VAR}_wt.pkl'); lw = {s['subpath']: fv(s) for s in lws}; del lws
lms, _ = joblib.load(LOCAL/f'{VAR}_mut.pkl'); lmi = {s['subpath']: i for i, s in enumerate(lms)}
need = np.unique(c.feature_row.to_numpy())
cov = [int(r) for r in need if prep[r]['mut_id'] in MUT and MUT[prep[r]['mut_id']] in lmi
       and prep[r]['wt_id'] in WT and WT[prep[r]['wt_id']] in lw]
c = c[c.feature_row.isin(set(cov))].reset_index(drop=True)
rows = np.unique(c.feature_row.to_numpy()); cm = {int(r): i for i, r in enumerate(rows)}
fi = np.array([cm[int(r)] for r in c.feature_row]); n = len(rows)
A = np.empty((n, D), 'float32'); src = np.memmap(DELTA, mode='r', dtype='float32', shape=(5729, D))
for i, r in enumerate(rows): A[i] = src[int(r)]
del src
B = np.empty((n, D), 'float32')
for i, r in enumerate(rows):
    r = int(r); B[i] = fv(lms[lmi[MUT[prep[r]['mut_id']]]]) - lw[WT[prep[r]['wt_id']]]
del lms, lmi, lw
tr = np.flatnonzero(c['set'] == 'train'); te = np.flatnonzero(c['set'] == 'test')
y = c.y.to_numpy(); sg = c.sign.to_numpy(); yt = y[te]
def gram(x, a, b, bi):
    off, width = OFF[bi], W[bi]
    k = np.zeros((len(a), len(a))); v = np.zeros((len(b), len(a)))
    for l in range(off, off+width, 2048):
        r = min(l+2048, off+width)
        z = np.asarray(x[fi[a], l:r], 'float64')*sg[a, None]
        mu = z.mean(0); sd = np.maximum(z.std(0), 1e-8); z = (z-mu)/sd
        t = (np.asarray(x[fi[b], l:r], 'float64')*sg[b, None]-mu)/sd
        k += z@z.T/width; v += t@z.T/width
    return k, v
WEIGHTS = [0, .125, .25, .5, 1, 2, 4, 8]; ALPHAS = np.logspace(-3, 3, 13)
surf = {}
for label, x in (('authors', A), ('local', B)):
    g = [gram(x, tr, te, bi) for bi in range(3)]
    S = np.zeros((len(WEIGHTS), len(ALPHAS)))
    for wi, w in enumerate(WEIGHTS):
        k = w*(g[0][0]+g[1][0])+g[2][0]; v = w*(g[0][1]+g[1][1])+g[2][1]
        ev, u = np.linalg.eigh(k); ev = np.maximum(ev, 0); ym = y[tr].mean(); rhs = u.T@(y[tr]-ym)
        for ai, al in enumerate(ALPHAS):
            S[wi, ai] = pearsonr(yt, v@(u@(rhs/(ev+al)))+ym).statistic
    surf[label] = S
    print(f"\n{label}: test Pearson over the (block weight, alpha) surface")
    print("   w \\ a " + "".join(f"{a:7.2g}" for a in ALPHAS))
    for wi, w in enumerate(WEIGHTS):
        print(f"  {w:6g} " + "".join(f"{S[wi,ai]:7.3f}" for ai in range(len(ALPHAS))))
    i, j = np.unravel_index(np.argmax(S), S.shape)
    print(f"  best {S[i,j]:.4f} at w={WEIGHTS[i]} alpha={ALPHAS[j]:.4g}")
sel = {'authors': (2, 10.0), 'local': (0.5, 10.0)}
print("\ncross-applying the weight each arm's inner CV selected:")
for label in ('authors', 'local'):
    for who, (w, al) in sel.items():
        wi = WEIGHTS.index(w); ai = int(np.argmin(abs(ALPHAS-al)))
        print(f"  {label:8s} features at {who}'s selection (w={w}, alpha={al:g}): r={surf[label][wi,ai]:.4f}")
json.dump({'weights': WEIGHTS, 'alphas': ALPHAS.tolist(),
           'authors': surf['authors'].tolist(), 'local': surf['local'].tolist()},
          open(out_dir()/'WEIGHT_SURFACE.json', 'w'), indent=1)
