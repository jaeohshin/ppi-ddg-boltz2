#!/usr/bin/env python3
"""How many of the 128 pair-block dimensions carry the signal, and are they distinct?

    python pair_block_dimensions.py                       # ./pair_block.npz
    python pair_block_dimensions.py --cache path/to/x.npz

Ranks dimensions by their own within-complex PCC, recomputed inside every
training fold, refits the top k, and bootstraps k=4 and k=1 against all 128.
In-fold ranking is the point: ranking on the full cohort inflates k=1 by 0.016
through selection over 128 candidates.  Also prints the correlation structure of
the strongest dimensions, which is where the participation ratio of 1.9 comes
from.  See docs/PAIR_BLOCK.md.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from scipy import stats
from prefold_probe.evaluate import per_structure_pcc, paired_bootstrap
from prefold_probe.reader import BlockRidge, _block_gram
from prefold_probe.splits import load_split
_ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
_ap.add_argument('--cache', type=Path, default=Path('pair_block.npz'),
                 help='the 128-wide pair block written by pair_block.py --cache')
_ap.add_argument('--out', type=Path, default=Path('nested128.json'))
_args = _ap.parse_args()
if not _args.cache.exists():
    raise SystemExit(
        f'{_args.cache} not found.\n'
        f'  Build it first:  python pair_block.py --pickles <dir> --cache {_args.cache}')
d=np.load(_args.cache,allow_pickle=True)
wz,y_all,pdbs_all=d['wz'],d['y'],d['pdbs']
ALPHAS=tuple(np.logspace(-2,5,29)); KS=(1,2,4,8,16,32,64,128)
out={}
for split in ('hom99','holdoutprot'):
    f=load_split(split,pdbs_all); keep=f>=0
    X=wz[keep]; y=y_all[keep]; pdbs=pdbs_all[keep]; folds=f[keep]
    oof={k:{a:np.zeros_like(y) for a in ALPHAS} for k in KS}
    picks={k:[] for k in KS}
    for fo in sorted(set(folds)):
        te=folds==fo; tr=~te
        # rank on TRAINING rows only: within-complex PCC over training complexes
        score=np.array([abs(per_structure_pcc(X[tr][:,j],y[tr],pdbs[tr])[0]) for j in range(X.shape[1])])
        order=np.argsort(-np.nan_to_num(score))
        for k in KS:
            sub=order[:k]; picks[k].append(sorted(int(j) for j in sub))
            K=_block_gram(X[:,sub],tr)
            fit=BlockRidge(alphas=ALPHAS).fit(K[np.ix_(tr,tr)],y[tr])
            for a,p in fit.predict_grid(K[np.ix_(te,tr)]).items(): oof[k][a][te]=p
    print(f"\n=== {split} ===  ranking recomputed per fold")
    best={}
    for k in KS:
        b=max(((per_structure_pcc(oof[k][a],y,pdbs)[0]),a) for a in ALPHAS)
        best[k]=b
        stab = len(set(map(tuple,picks[k])))
        print("  k=%3d  pcc %.4f (alpha %-8.4g)  distinct top-k sets across the 3 folds: %d"%(k,b[0],b[1],stab))
    # is k=4 distinguishable from k=128?
    a4=best[4][1]; a128=best[128][1]
    bt=paired_bootstrap(oof[4][a4],oof[128][a128],y,pdbs,n_boot=10000,seed=20261115)
    print("  k=4 minus k=128: %+.4f [%+.4f, %+.4f] over %d complexes"%(bt['delta'],bt['lo'],bt['hi'],bt['n_groups']))
    a1=best[1][1]
    bt1=paired_bootstrap(oof[1][a1],oof[128][a128],y,pdbs,n_boot=10000,seed=20261115)
    print("  k=1 minus k=128: %+.4f [%+.4f, %+.4f]"%(bt1['delta'],bt1['lo'],bt1['hi']))
    out[split]={'best':{str(k):[float(v[0]),float(v[1])] for k,v in best.items()},
                'top4_per_fold':picks[4],'k4_vs_k128':bt,'k1_vs_k128':bt1}
# which dimensions, and are they duplicates of each other?
f=load_split('holdoutprot',pdbs_all); keep=f>=0
X=wz[keep]
# Derived, not hard-coded: naming a fixed index would contradict the finding
# that the identity of the strongest dimension is not stable across folds.
rank=sorted(range(X.shape[1]),key=lambda j:-abs(per_structure_pcc(X[:,j],y_all[keep],pdbs_all[keep])[0]))
top=rank[:8]
Z=(X[:,top]-X[:,top].mean(0))/np.maximum(X[:,top].std(0),1e-8)
C=np.corrcoef(Z.T)
print("\n=== correlation among the 8 strongest dimensions ===")
print("      "+" ".join("%6d"%t for t in top))
for i,t in enumerate(top): print("%5d "%t+" ".join("%+6.2f"%C[i,j] for j in range(len(top))))
sv=np.linalg.svd(Z,compute_uv=False)
print("participation ratio of just these 8: %.2f"%(((sv**2).sum()**2)/(sv**4).sum()))
json.dump(out,open(_args.out,'w'),indent=1,default=float)
print(f'wrote {_args.out}')
