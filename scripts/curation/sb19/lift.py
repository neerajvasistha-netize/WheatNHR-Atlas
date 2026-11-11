"""Anchor-based liftover reproducing the atlas pipeline for new records.
Anchors: (1) atlas records already lifted by the pipeline (v1.0 position -> each assembly);
(2) for RefSeq v2.1 additionally the gene anchors of sr_v21 (v2.1 position) joined to Ensembl
RefSeq v1.0 gene coordinates. Interpolates between the nearest anchors on each side."""
import json, gzip, bisect, csv
from collections import defaultdict
import os
W = os.path.dirname(os.path.abspath(__file__)) + '/'
import re as _re
_h = open(W + '../../../web/index.html').read()   # the built page carries every record's pipeline liftover
seed = json.loads(_re.search(r'<script type="application/json" id="seed">(.*?)</script>', _h, _re.S).group(1).replace('<\\/', '</'))
num = lambda v: None if v in (None, '') else float(v)
ASMS = ['RefSeq v2.1', 'Alchemy', 'ArinaLrFor', 'Jagger', 'Julius', 'Kariega', 'Lancer', 'Landmark', 'Mace',
        'SY Mattis', 'Norin 61', 'Paragon', 'Renan', 'Stanley']
anc = {a: defaultdict(list) for a in ASMS}
for r in seed['reports']:
    if str(r.get('ref')) not in ('RefSeq v1.0', 'RefSeq v1.x', 'RefSeq v1.1'): continue
    for a in ASMS:
        t = (r.get('asm') or {}).get(a)
        if not t or t.get('flag') not in ('ok',) or t.get('chr') != r['chr']: continue
        for k in ('start', 'end', 'peak'):
            if num(r.get(k)) is not None and t.get(k) is not None:
                anc[a][r['chr']].append((num(r[k]), float(t[k])))
v10 = {}
for l in gzip.open(W + 'genes_v10.tsv.gz', 'rt'):
    g, c, s, e = l.split('\t'); v10[g] = (c, (int(s) + int(e)) / 2e6)
S = json.load(gzip.open('' + W + '../../../data/sr_v21.json.gz'))
_rec = {c: sorted(v) for c, v in anc['RefSeq v2.1'].items()}
def _interp(c, x):
    L = _rec.get(c) or []; xs = [p[0] for p in L]; i = bisect.bisect_left(xs, x)
    if i == 0 or i == len(L): return None
    lo, hi = L[i - 1], L[i]
    if hi[0] - lo[0] <= 0 or hi[1] < lo[1]: return None
    return lo[1] + (hi[1] - lo[1]) * (x - lo[0]) / (hi[0] - lo[0])
GENE_KEPT = GENE_DROPPED = 0
for g in S['genes']:
    x = v10.get(g['gene_v1.1'])
    if x and x[0] == g['chr']:
        y = _interp(g['chr'], x[1])
        if y is not None and abs(y - g['pos_Mb_v2.1']) <= 1.5:
            anc['RefSeq v2.1'][g['chr']].append((x[1], g['pos_Mb_v2.1'])); GENE_KEPT += 1
        else: GENE_DROPPED += 1
OUTLIERS = []
for a in ASMS:
    for c in anc[a]:
        d = {}
        for x, y in anc[a][c]: d.setdefault(round(x, 4), y)
        L = sorted(d.items()); keep = []
        for i, (x, y) in enumerate(L):   # drop anchors off the local collinear trend
            nb = [q[1] - q[0] for q in L[max(0, i - 6):i] + L[i + 1:i + 7]]
            if nb:
                nb.sort(); med = nb[len(nb) // 2]
                if abs((y - x) - med) > 5: OUTLIERS.append((a, c, x, y)); continue
            keep.append((x, y))
        anc[a][c] = keep

def lift_point(a, c, x, maxgap=8.0):
    L = anc[a].get(c) or []
    if not L: return None, 'no_anchor'
    xs = [p[0] for p in L]
    i = bisect.bisect_left(xs, x)
    lo = L[i - 1] if i > 0 else None
    hi = L[i] if i < len(L) else None
    if lo and hi and hi[0] - lo[0] > 0:
        gap = min(x - lo[0], hi[0] - x)
        slope = (hi[1] - lo[1]) / (hi[0] - lo[0])
        y = lo[1] + slope * (x - lo[0])
        if slope >= 0 and gap > maxgap:   # sparse anchors, but consistent on both sides
            return round(y, 3), 'rechecked'
        if slope < 0:   # local inversion between assemblies: take the nearest anchor's offset
            near = lo if x - lo[0] <= hi[0] - x else hi
            return round(near[1] + (x - near[0]), 3), 'reoriented_checked'
        return round(y, 3), 'ok' if gap <= maxgap else 'wide_gap'
    near = lo or hi
    if abs(near[0] - x) > maxgap: return round(near[1] + (x - near[0]), 3), 'wide_gap'
    return round(near[1] + (x - near[0]), 3), 'ok'

def lift_record(r):
    """r has chr,start,end,peak,ref (v1.0 or v2.x). Returns asm dict like the pipeline."""
    c = r['chr']; s, e, p = num(r.get('start')), num(r.get('end')), num(r.get('peak'))
    if s is None and p is None: return {'_note': 'no physical position in source'}
    asm = {}
    v2 = str(r.get('ref', '')).startswith('RefSeq v2')
    if v2:
        # back-convert to v1.0 using the inverse of the v2.1 anchors
        inv = sorted((y, x) for x, y in anc['RefSeq v2.1'][c])
        def back(y):
            ys = [q[0] for q in inv]; i = bisect.bisect_left(ys, y)
            lo = inv[i - 1] if i > 0 else inv[0]; hi = inv[i] if i < len(inv) else inv[-1]
            if hi[0] == lo[0]: return round(lo[1] + (y - lo[0]), 3)
            return round(lo[1] + (hi[1] - lo[1]) * (y - lo[0]) / (hi[0] - lo[0]), 3)
        o = {'chr': c, 'flag': 'as_reported' if r['ref'] == 'RefSeq v2.1' else 'v2.0_as_reported'}
        for k, v in (('start', s), ('end', e), ('peak', p)):
            if v is not None: o[k] = v
        asm['RefSeq v2.1'] = o
        s1 = back(s) if s is not None else None; e1 = back(e) if e is not None else None; p1 = back(p) if p is not None else None
        o1 = {'chr': c, 'flag': 'converted from ' + r['ref']}
        for k, v in (('start', s1), ('end', e1), ('peak', p1)):
            if v is not None: o1[k] = v
        asm['RefSeq v1.0'] = o1
        s, e, p = s1, e1, p1
    for a in ASMS:
        if a == 'RefSeq v2.1' and v2: continue
        o = {'chr': c}; flags = []
        for k, v in (('start', s), ('end', e), ('peak', p)):
            if v is None: continue
            y, f = lift_point(a, c, v); flags.append(f)
            if y is not None: o[k] = y
        if 'start' in o and 'end' in o and o['start'] > o['end']:
            o['start'], o['end'] = o['end'], o['start']; flags.append('inverted')
        if flags and all(f in ('ok', 'rechecked') for f in flags) and 'rechecked' in flags:
            flags = ['rechecked']; o['qc'] = 'interpolated between the nearest consistent anchors on both sides (anchor gap above 8 Mb)'
        o['flag'] = 'rechecked' if flags == ['rechecked'] else 'reoriented_checked' if 'reoriented_checked' in flags and 'wide_gap' not in flags and 'inverted' not in flags else 'inverted' if 'inverted' in flags else ('wide_gap' if 'wide_gap' in flags else ('ok' if flags and all(f == 'ok' for f in flags) else (flags[0] if flags else 'no_anchor')))
        asm[a] = o
    return asm

if __name__ == '__main__':
    # self-test: re-lift existing pipeline records and compare
    import random
    random.seed(1); errs = []
    recs = [r for r in seed['reports'] if r.get('ref') == 'RefSeq v1.0' and (r.get('asm') or {}).get('RefSeq v2.1', {}).get('flag') == 'ok' and r.get('peak') is not None]
    for r in random.sample(recs, 300):
        y, f = lift_point('RefSeq v2.1', r['chr'], float(r['peak']))
        errs.append(abs(y - r['asm']['RefSeq v2.1']['peak']))
    errs.sort(); print('v2.1 abs err median %.3f  p95 %.3f  max %.3f Mb' % (errs[150], errs[285], errs[-1]))
    for a in ASMS: print(a, sum(len(v) for v in anc[a].values()))
