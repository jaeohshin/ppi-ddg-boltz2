"""Exact join between the authors' feature records and the local Boltz re-fold.

The two producers key their records differently in two ways:
  * numbering  - the authors' prepared cohort uses SKEMPI's `Mutation(s)_cleaned`
                 positions, the local re-fold uses `Mutation(s)_PDB` positions;
  * token order- the authors sort the mutation list, the local re-fold keeps
                 SKEMPI's original order.
SKEMPI lists both numbering columns in the same token order for every row, so
zipping them gives a per-(complex, chain) position dictionary.  Here the local
keys are normalised into the authors' convention (translate each position, then
sort), which makes the join a dictionary lookup.

Eight antibody complexes also differ only in how the chain groups are spelled
(1DAN_HL_TU vs 1DAN_HL_UT), handled by chain_variants().
"""
from paths import (EXP, AUTHOR, LOCAL, VARIANT, BASE, PREP, DELTA, SKEMPI,
                   BLIND_COHORT, WIDTHS, D, OFFSETS, N_PREPARED_ROWS,
                   join_map, out_dir, missing)
VAR = VARIANT
import csv, json, re, collections, joblib

TOK = re.compile(r'^([ACDEFGHIKLMNPQRSTVWY])([A-Za-z0-9])(-?\d+[A-Za-z]?)([ACDEFGHIKLMNPQRSTVWY])$')

to_cleaned = {}   # (pdb, chain, pdb_pos) -> cleaned_pos
conflict = 0
for r in csv.DictReader(open(SKEMPI), delimiter=';'):
    a, b = r['Mutation(s)_PDB'].split(','), r['Mutation(s)_cleaned'].split(',')
    if len(a) != len(b): continue
    for x, y in zip(a, b):
        mx, my = TOK.match(x), TOK.match(y)
        if not mx or not my: continue
        if (mx.group(1), mx.group(2), mx.group(4)) != (my.group(1), my.group(2), my.group(4)):
            continue
        k = (r['#Pdb'], mx.group(2), mx.group(3))
        if k in to_cleaned and to_cleaned[k] != my.group(3): conflict += 1
        to_cleaned[k] = my.group(3)
print(f"position dictionary: {len(to_cleaned)} entries, {conflict} conflicts")

def chain_variants(prefix):
    p = prefix.split('_'); out = [prefix]
    if len(p) >= 5:
        for i in (3, 4):
            q = list(p); q[i] = q[i][::-1]; out.append('_'.join(q))
        q = list(p); q[3], q[4] = q[4], q[3]; out.append('_'.join(q))
        q = list(p); q[3], q[4] = q[4][::-1], q[3][::-1]; out.append('_'.join(q))
    return out

def normalise(subpath):
    """Local key -> the authors' convention, or None if a position is unknown."""
    prefix, _, tail = subpath.rpartition('_')
    pdb = prefix[len('boltz_results_'):]
    for cand_pdb, cand_prefix in [(pdb, prefix)] + [(p[len('boltz_results_'):], p)
                                                    for p in chain_variants(prefix)[1:]]:
        toks = []
        for t in tail.split(','):
            m = TOK.match(t)
            if not m: return None
            p = to_cleaned.get((cand_pdb, m.group(2), m.group(3)))
            if p is None: toks = None; break
            toks.append(m.group(1) + m.group(2) + p + m.group(4))
        if toks: return f"{cand_prefix}_{','.join(sorted(toks))}", cand_prefix
    return None

for kind in ('mut', 'wt'):
    recs, _ = joblib.load(LOCAL / f'{VAR}_{kind}.pkl')
    local = [s['subpath'] for s in recs]; del recs
    if kind == 'wt':
        index = {}
        for s in local:
            for v in chain_variants(s): index[v] = s
    else:
        index = {}; bad = 0
        for s in local:
            n = normalise(s)
            if n is None: bad += 1; continue
            key, pre = n
            index[key] = s
            # also register the chain-group spellings the authors might use
            for v in chain_variants(pre):
                index[key.replace(pre, v, 1)] = s
        print(f"local mut: {len(local)} records, {bad} unnormalisable, {len(index)} keys")
    globals()[kind + '_index'] = index

mut_map, wt_map = {}, {}
unresolved = collections.Counter()
for r in csv.DictReader(open(PREP)):
    if r['mut_id'] in mut_index: mut_map[r['mut_id']] = mut_index[r['mut_id']]
    else: unresolved['mut'] += 1
    if r['wt_id'] not in wt_map:
        if r['wt_id'] in wt_index: wt_map[r['wt_id']] = wt_index[r['wt_id']]
        else: unresolved['wt'] += 1
print(f"joined: {len(mut_map)} mutant records, {len(wt_map)} wild-type records")
print(f"unresolved: {dict(unresolved)}")
json.dump({'mut': mut_map, 'wt': wt_map}, open(join_map(), 'w'))
