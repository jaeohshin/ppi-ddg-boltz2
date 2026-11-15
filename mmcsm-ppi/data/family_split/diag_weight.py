"""DIAGNOSTIC ONLY.  Rule A chose w=0.25 for all 34 families -- the bottom of the
grid.  This recomputes the leave-one-family-out predictions over an EXTENDED
weight range, including 0 (the pair block alone), to see whether 0.25 is a real
optimum or a boundary artefact.  It selects nothing; it scores the held-out
families across a grid, exactly as penalty_sensitivity.json does for the blind
set, and is reported as a diagnostic."""
import os
os.environ.update(OPENBLAS_NUM_THREADS='8',OMP_NUM_THREADS='8',MKL_NUM_THREADS='8')
from pathlib import Path
import numpy as np, pandas as pd, json, time, collections
from scipy.stats import spearmanr, pearsonr
O=Path(__file__).resolve().parent; E=O.parent; t0=time.time()
WIDTHS=[147456,147456,128]; P=sum(WIDTHS)
WEIGHTS=[0.0,0.03125,0.0625,0.125,0.25,0.5,1,2,4]; ALPHAS=np.logspace(-3,3,13)
FAM=json.load(open('/data/work/prefold-dg-linear-probe/mmcsm-ppi/data/families_0.99.json'))['family_of']
c=pd.read_csv(E/'mmcsm_ridge_blind_20261115/cohort.csv')
tr=np.flatnonzero(c['set']=='train'); te=np.flatnonzero(c['set']=='test')
y=c.y.to_numpy(); sg=c.sign.to_numpy(); fi=c.feature_row.to_numpy()
pdbl=np.array([p.lower() for p in c.pdb.to_numpy()]); fam=np.array([FAM.get(p,p) for p in pdbl])
x=np.memmap(E/'boltz_clean_baseline_20261101/ridge_full_v1/delta_features.f32',mode='r',dtype='float32',shape=(5729,P))
rows=np.concatenate([fi[tr],fi[te]]); F=np.empty((len(rows),P),dtype='float32')
for i,f in enumerate(rows): F[i]=x[f]
del x
F*=np.concatenate([sg[tr],sg[te]]).astype('float32')[:,None]
TR=np.arange(len(tr)); TE=np.arange(len(tr),len(rows)); ytr=y[tr]; yte=y[te]
ftr=fam[tr]; fte=fam[te]; paper=c.iloc[te].paper.to_numpy()
def grams(a,b):
 ks=[];vs=[];off=0
 for w in WIDTHS:
  k=np.zeros((len(a),len(a))); v=np.zeros((len(b),len(a)))
  for l in range(off,off+w,2048):
   r=min(l+2048,off+w)
   z=np.asarray(F[a,l:r],dtype='float64'); mu=z.mean(0); sd=np.maximum(z.std(0),1e-8); z=(z-mu)/sd
   k+=z@z.T/w; v+=((np.asarray(F[b,l:r],dtype='float64')-mu)/sd)@z.T/w
  ks.append(k); vs.append(v); off+=w
 return ks,vs
fams=sorted(set(fte)); pred=np.full((len(WEIGHTS),len(ALPHAS),len(te)),np.nan)
for n,f in enumerate(fams):
 h=np.flatnonzero(fte==f); pool=TR[ftr[TR]!=f]
 ks,vs=grams(pool,TE[h])
 for wi,w in enumerate(WEIGHTS):
  k=w*(ks[0]+ks[1])+ks[2]; v=w*(vs[0]+vs[1])+vs[2]
  ev,u=np.linalg.eigh(k); ev=np.maximum(ev,0); ym=ytr[pool].mean(); rhs=u.T@(ytr[pool]-ym)
  for ai,a in enumerate(ALPHAS): pred[wi,ai,h]=v@(u@(rhs/(ev+a)))+ym
 print(f"[{time.time()-t0:5.0f}s] {n+1}/{len(fams)} {f}",flush=True)
assert np.isfinite(pred).all()
R=np.array([[pearsonr(yte,pred[i,j]).statistic for j in range(len(ALPHAS))] for i in range(len(WEIGHTS))])
S=np.array([[spearmanr(yte,pred[i,j]).statistic for j in range(len(ALPHAS))] for i in range(len(WEIGHTS))])
np.savez(O/'lofo_weight_surface.npz',pearson=R,spearman=S,weights=WEIGHTS,alphas=ALPHAS)
print("\nfamily-disjoint PEARSON (rows=weight, cols=alpha)")
print("          "+" ".join("%7.3g"%a for a in ALPHAS))
for i,w in enumerate(WEIGHTS): print("w=%-7g "%w+" ".join("%7.4f"%v for v in R[i]))
bi=np.unravel_index(np.nanargmax(R),R.shape); bj=np.unravel_index(np.nanargmax(S),S.shape)
print("\nbest pearson  w=%g alpha=%.4g -> %.4f (mmCSM 0.7019)"%(WEIGHTS[bi[0]],ALPHAS[bi[1]],R[bi]))
print("best spearman w=%g alpha=%.4g -> %.4f (mmCSM 0.6503)"%(WEIGHTS[bj[0]],ALPHAS[bj[1]],S[bj]))
print("rule A chose  w=0.25 alpha=3.162 -> pearson %.4f spearman %.4f"%(R[WEIGHTS.index(0.25),8],S[WEIGHTS.index(0.25),8]))
print("pair block alone w=0 best pearson %.4f at alpha=%.4g"%(R[0].max(),ALPHAS[R[0].argmax()]))
