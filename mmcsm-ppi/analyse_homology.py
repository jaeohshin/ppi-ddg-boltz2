"""Cluster the 147 mmCSM-PPI complexes into protein families, then score both
methods on the blind complexes whose family the training set never contained."""
import csv, glob, os, sys, time, collections, json
import numpy as np
from scipy import stats

USAGE = """usage: analyse_homology.py <release dir> <fasta dir> <identity> [predictions.csv] [out dir]

Clusters the complexes by chain identity and grades the blind set by what the
training set contained.  The clustering is the slow part (all-pairs global
alignment), so it is cached as <out dir>/families_<identity>.json and reused."""
if len(sys.argv) < 4: sys.exit(USAGE)
HERE = os.path.dirname(os.path.abspath(__file__))
D, FA, TH = sys.argv[1], sys.argv[2], float(sys.argv[3])
PRED = sys.argv[4] if len(sys.argv) > 4 else f"{HERE}/data/predictions_boltz_block_ridge_cohort1237.csv"
OUT = sys.argv[5] if len(sys.argv) > 5 else f"{HERE}/data"
os.makedirs(OUT, exist_ok=True)
CACHE = f"{OUT}/families_{TH}.json"
T0 = time.time()
def log(m): print(f"[{time.time()-T0:5.0f}s] {m}", flush=True)
def pid(s): return s.split("-")[0].strip().lower()

def fasta(p):
    out, cur = [], []
    for ln in open(p):
        if ln.startswith(">"):
            if cur: out.append("".join(cur)); cur = []
        else: cur.append(ln.strip())
    if cur: out.append("".join(cur))
    return [s for s in out if len(s) >= 20]

def cluster():
    from Bio import Align
    from Bio.Align import substitution_matrices
    chains = {os.path.basename(p)[:-6]: fasta(p) for p in glob.glob(f"{FA}/*.fasta")}
    uniq = sorted({c for v in chains.values() for c in v})
    idx = {c: i for i, c in enumerate(uniq)}
    par = list(range(len(uniq)))
    def find(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    log(f"{len(chains)} entries, {len(uniq)} distinct chain sequences")

    al = Align.PairwiseAligner(mode="global")
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score, al.extend_gap_score = -11, -1
    n = 0
    for i in range(len(uniq)):
        for j in range(i + 1, len(uniq)):
            if find(i) == find(j): continue
            a, b = uniq[i], uniq[j]
            if not 0.5 <= len(a) / len(b) <= 2.0: continue
            n += 1
            aln = al.align(a, b)[0]
            m = sum(1 for x, y in zip(aln[0], aln[1]) if x == y and x != "-")
            if m / min(len(a), len(b)) >= TH:
                par[find(i)] = find(j)
        if i % 50 == 0: log(f"  {i}/{len(uniq)} ({n} alignments)")
    log(f"{n} alignments at identity >= {TH}")
    return chains, uniq, idx, find

if os.path.exists(CACHE):
    fam = json.load(open(CACHE))["family_of"]
    log(f"{len(set(fam.values()))} families among {len(fam)} entries (cached: {CACHE})")
else:
    chains, uniq, idx, find = cluster()
    ep = {p: p for p in chains}
    def root(x):
        while ep[x] != x: ep[x] = ep[ep[x]]; x = ep[x]
        return x
    bc = collections.defaultdict(list)
    for p, cs in chains.items():
        for c in cs: bc[find(idx[c])].append(p)
    for mem in bc.values():
        for m in mem[1:]:
            a, b = root(mem[0]), root(m)
            if a != b: ep[a] = b
    fam = {p: root(p) for p in chains}
    log(f"{len(set(fam.values()))} families among {len(fam)} entries")

def rows(p): return [r for r in list(csv.reader(open(p)))[1:] if r]
train_pdb = {pid(r[0]) for r in rows(f"{D}/train_set.raw")}
train_fam = {fam[p] for p in train_pdb if p in fam}
blind = [dict(pdb=pid(r[0]), muts=r[1], y=float(r[2]), mm=float(r[5])) for r in rows(f"{D}/blindtest.raw")]
mine = {}
for r in csv.DictReader(open(PRED)):
    mine[(pid(r["pdb"]), tuple(sorted(m.strip() for m in r["mutations"].split(","))))] = float(r["prediction"])
for b in blind:
    b["ridge"] = mine.get((b["pdb"], tuple(sorted(m.strip() for m in b["muts"].split(",")))), np.nan)

def sc(rs, label):
    y = np.array([r["y"] for r in rs]); mm = np.array([r["mm"] for r in rs]); rg = np.array([r["ridge"] for r in rs])
    g = np.array([r["pdb"] for r in rs]); k = np.isfinite(rg)
    out = [f"  {label:34s} n={len(rs):3d}/{k.sum():3d} cplx={len(set(g)):2d}"]
    out.append(f"mmCSM r={stats.pearsonr(y,mm)[0]:+.3f} rho={stats.spearmanr(y,mm)[0]:+.3f} RMSE={np.sqrt(((y-mm)**2).mean()):.2f}")
    if k.sum() >= 5:
        out.append(f"| ridge r={stats.pearsonr(y[k],rg[k])[0]:+.3f} rho={stats.spearmanr(y[k],rg[k])[0]:+.3f} RMSE={np.sqrt(((y[k]-rg[k])**2).mean()):.2f}")
        rng = np.random.default_rng(0); gs = np.unique(g[k]); ds = []
        for _ in range(4000):
            m = np.concatenate([np.where(g[k] == q)[0] for q in rng.choice(gs, len(gs))])
            if len(np.unique(y[k][m])) > 2:
                ds.append(stats.pearsonr(y[k][m], rg[k][m])[0] - stats.pearsonr(y[k][m], mm[k][m])[0])
            if len(ds) > 1500: break
        if ds:
            lo, hi = np.percentile(ds, [2.5, 97.5])
            d = stats.pearsonr(y[k], rg[k])[0] - stats.pearsonr(y[k], mm[k])[0]
            out.append(f"| delta={d:+.3f} [{lo:+.3f},{hi:+.3f}]")
    print("  ".join(out))

seen_pdb = [b for b in blind if b["pdb"] in train_pdb]
new_pdb = [b for b in blind if b["pdb"] not in train_pdb]
new_fam = [b for b in blind if b["pdb"] in fam and fam[b["pdb"]] not in train_fam]
print(f"\nBLIND TEST graded by what training contained (identity >= {TH})")
sc(blind, "all blind")
sc(seen_pdb, "same PDB entry as training")
sc(new_pdb, "new PDB entry")
sc(new_fam, "NEW FAMILY")
print(f"\n  blind complexes {len({b['pdb'] for b in blind})}; new PDB entry "
      f"{len({b['pdb'] for b in new_pdb})}; new family {len({b['pdb'] for b in new_fam})}")
print(f"  training: {len(train_pdb)} entries -> {len(train_fam)} families")
json.dump({"threshold": TH, "family_of": fam,
           "train_families": sorted(train_fam),
           "blind_new_family_pdbs": sorted({b['pdb'] for b in new_fam})},
          open(CACHE, "w"), indent=1)
