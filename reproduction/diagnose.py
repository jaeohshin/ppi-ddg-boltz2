"""Where does the authors' / local disagreement live?

Step 0 validates the join itself against the labels stored inside both pkls -
if that fails, nothing else here means anything.
Then the disagreement is split three ways: the wild-type vector, the mutant
vector, and the mut-minus-wt difference the ridge actually reads; and per block.
"""
from paths import (EXP, AUTHOR, LOCAL, VARIANT, BASE, PREP, DELTA, SKEMPI,
                   BLIND_COHORT, WIDTHS, D, OFFSETS, N_PREPARED_ROWS,
                   join_map, out_dir, missing)
VAR = VARIANT
import os, sys, json, csv, time
os.environ.update(OPENBLAS_NUM_THREADS='8', OMP_NUM_THREADS='8', MKL_NUM_THREADS='8')
from pathlib import Path

import numpy as np, pandas as pd, joblib

AUTH = AUTHOR
W = [147456, 147456, 128]; D = sum(W); NAMES = ['wsi', 'ws', 'wz']
SL = [slice(0, W[0]), slice(W[0], W[0]+W[1]), slice(W[0]+W[1], D)]
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
lws, lwy = joblib.load(LOCAL / f'{VAR}_wt.pkl')
lw = {s['subpath']: fv(s) for s in lws}; lwl = {s['subpath']: float(y) for s, y in zip(lws, lwy)}
del lws, lwy
lms, lmy = joblib.load(LOCAL / f'{VAR}_mut.pkl')
lmi = {s['subpath']: i for i, s in enumerate(lms)}
lml = {s['subpath']: float(y) for s, y in zip(lms, lmy)}

rows = np.unique(c.feature_row.to_numpy())
rows = np.array([r for r in rows if prep[r]['mut_id'] in MUT and prep[r]['wt_id'] in WT])
log(f"{len(rows)} joined feature rows under test")

# ---- step 0: does the join agree with the labels stored in both pkls? --------
pa = np.array([float(prep[r]['mut_logkd']) for r in rows])
pl = np.array([lml[MUT[prep[r]['mut_id']]] for r in rows])
wa = np.array([float(prep[r]['wt_logkd']) for r in rows])
wl = np.array([lwl[WT[prep[r]['wt_id']]] for r in rows])
log(f"mutant label: max |authors - local| = {np.max(abs(pa-pl)):.3e}; "
    f"exact within 1e-6: {int((abs(pa-pl)<1e-6).sum())}/{len(rows)}")
log(f"wild-type label: max |authors - local| = {np.max(abs(wa-wl)):.3e}; "
    f"exact within 1e-6: {int((abs(wa-wl)<1e-6).sum())}/{len(rows)}")
bad = np.flatnonzero(abs(pa-pl) > 1e-6)
if len(bad): log(f"  label mismatches on {len(bad)} rows, e.g. {[prep[rows[i]]['mut_id'] for i in bad[:3]]}")

# ---- step 1: where is the disagreement - wt, mut, or the difference? ---------
src = np.memmap(DELTA, mode='r', dtype='float32', shape=(5729, D))
log("loading the authors' wt pkl for the raw-vector comparison")
aws, _ = joblib.load(AUTH / 'skempi_wt.pkl')
aw = {s['subpath']: fv(s) for s in aws}; del aws

def cos(a, b): return float(np.dot(a, b) / (np.linalg.norm(a)*np.linalg.norm(b) + 1e-12))
stats = {k: {n: [] for n in NAMES + ['all']} for k in ['wt', 'mut', 'delta']}
norms = {k: {n: [[], []] for n in NAMES} for k in ['wt', 'delta']}
sub = rows if len(rows) <= 400 else rows[np.linspace(0, len(rows)-1, 400).astype(int)]
for r in sub:
    r = int(r)
    aw_v = aw[prep[r]['wt_id']]; lw_v = lw[WT[prep[r]['wt_id']]]
    ad_v = np.asarray(src[r], dtype=np.float32)
    lm_v = fv(lms[lmi[MUT[prep[r]['mut_id']]]])
    ld_v = lm_v - lw_v
    am_v = ad_v + aw_v
    for key, (a, b) in [('wt', (aw_v, lw_v)), ('mut', (am_v, lm_v)), ('delta', (ad_v, ld_v))]:
        for n, s in zip(NAMES, SL): stats[key][n].append(cos(a[s], b[s]))
        stats[key]['all'].append(cos(a, b))
    for n, s in zip(NAMES, SL):
        norms['wt'][n][0].append(float(np.linalg.norm(aw_v[s]))); norms['wt'][n][1].append(float(np.linalg.norm(lw_v[s])))
        norms['delta'][n][0].append(float(np.linalg.norm(ad_v[s]))); norms['delta'][n][1].append(float(np.linalg.norm(ld_v[s])))

print()
log(f"median cosine, authors vs local, over {len(sub)} rows")
print(f"  {'':8s} {'wsi':>8s} {'ws':>8s} {'wz':>8s} {'all':>8s}")
for key in ['wt', 'mut', 'delta']:
    print(f"  {key:8s}" + "".join(f" {np.median(stats[key][n]):8.4f}" for n in NAMES + ['all']))
print()
log("median L2 norm (authors / local)")
for key in ['wt', 'delta']:
    print(f"  {key:8s}" + "".join(
        f" {n}: {np.median(norms[key][n][0]):8.3f}/{np.median(norms[key][n][1]):8.3f}" for n in NAMES))
json.dump({'cosine_median': {k: {n: float(np.median(v)) for n, v in d.items()} for k, d in stats.items()},
           'cosine_p05': {k: {n: float(np.percentile(v, 5)) for n, v in d.items()} for k, d in stats.items()},
           'norm_median': {k: {n: [float(np.median(v[0])), float(np.median(v[1]))] for n, v in d.items()}
                           for k, d in norms.items()},
           'label_check': {'mut_max_abs_diff': float(np.max(abs(pa-pl))),
                           'wt_max_abs_diff': float(np.max(abs(wa-wl))),
                           'mut_exact': int((abs(pa-pl)<1e-6).sum()), 'n': len(rows)}},
          open(out_dir() / 'DIAGNOSE.json', 'w'), indent=1)
