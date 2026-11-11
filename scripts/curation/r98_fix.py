"""Atlas-wide re-check of loci with GWM markers using the primer sequences of Roder et al. 1998
(Genetics 149:2007-2023, Appendix; cross-checked with Roder et al. 1998 Genome 41:278-283, Table 2).
Every marker in a record's marker field is placed on RefSeq v1.0 with the same resolver as sb26_fix.py
(gene > exact-match in-silico PCR on the stated chromosome > marker dictionary). A record is corrected when
its current position lies more than TOL Mb outside the span of its own placed markers.
Run from the repository root: python scripts/curation/r98_fix.py  (writes work/r98_changes.json)"""
import json, re, sys, io, contextlib
sys.path.insert(0, 'scripts/curation'); sys.path.insert(0, 'scripts/curation/sb19')
with contextlib.redirect_stdout(io.StringIO()):
    import sb26_fix as F
import lift
TOL, FAR = 5.0, 20.0
RODER = F.RODER
MQTL = {'Liu2020MQTL_TS', 'Sharma2024FPS', 'Z_1228', 'Pundir2023BMC'}   # meta-QTL projections: not re-placed from single markers
base = lambda m: re.sub(r'[a-z]$', '', F.nk(m))
num = lambda x: float(x) if x not in (None, '') else None

def markers(txt):
    txt = re.sub(r'\((?:peak|also|closest|near)[^)]*\)', ' ', str(txt or ''), flags=re.I)
    out = []
    for t in re.split(r'[–—/;,|]|\s-\s|\s+and\s+|\s+', txt):
        t = t.strip(' .()[]')
        if re.match(r'^X?(gwm|wms|wmc|barc|cfd|cfa|gdm|psp|wPt|tPt|IWB|IWA|BS0|Ex_|Excalibur|Kukri|RAC875|BobWhite|Tdurum|wsnp)', t, re.I):
            out.append(re.sub(r'\.\d+$', '', t))
    return out

def v10(r):
    if r.get('ref') == 'RefSeq v1.0': return r
    return (r.get('asm') or {}).get('RefSeq v1.0') or lift.lift_record(r).get('RefSeq v1.0') or {}

def run(reports):
    ch = []
    for r in reports:
        if r['id'].startswith(('S19-', 'S26-')) or r.get('deleted') or r['source'] in MQTL: continue
        ms = markers(r.get('marker'))
        if not any(re.match(r'^X?(gwm|wms)', m, re.I) for m in ms): continue
        c = r.get('chr'); v = v10(r)
        cur = [x for x in (num(v.get('start')), num(v.get('end')), num(v.get('peak'))) if x is not None]
        got = []
        for m in ms:
            p, how = F.pos(m, c)
            if p is None: continue
            if re.match(r'^X?(gwm|wms)', m, re.I) and how == 'in-silico PCR' and c not in RODER.get(base(m), {}).get('chrs', []):
                continue   # product on a chromosome Roder et al. 1998 did not map this primer set to
            got.append((m, round(p, 4), how))
        g = [x for x in got if re.match(r'^X?(gwm|wms)', x[0], re.I) and x[2] == 'in-silico PCR']
        if not g or not cur: continue
        if len({base(x[0]) for x in got}) == 1 and len(got) > 1: got = got[:1]   # gwm639a/gwm639c: one primer set
        lo, hi = min(x[1] for x in got), max(x[1] for x in got)
        if len(got) >= 2 and hi - lo > 150: continue          # flanks too far apart to judge
        clo, chi = min(cur), max(cur)
        if len(got) >= 2:
            if clo >= lo - TOL and chi <= hi + TOL: continue
            new = dict(start=lo, end=hi, grade='B')
        else:
            if clo - TOL <= lo <= chi + TOL or abs((clo + chi) / 2 - lo) <= FAR: continue
            new = dict(peak=lo, grade='C')
        if r.get('name') == 'Yr46' and c == '4D':   # cloned (Lr67/Yr46 = TraesCS4D02G364400, Moore et al. 2015)
            new = dict(start=509.2823, end=509.2843, grade='B')
            got = got + [('Lr67/Yr46 gene TraesCS4D02G364400', 509.2823, 'Ensembl RefSeq v1.0 gene model; position used')]
        ch.append(dict(id=r['id'], source=r['source'], name=r.get('name'), chr=c, ref=r.get('ref'), marker=r.get('marker'),
                       old=[clo, chi], new=new, placed=got,
                       note='Re-checked with Roder et al. 1998 GWM primers (exact-match in-silico PCR on RefSeq v1.0), Nov 2026: '
                            + '; '.join('%s %.3f Mb (%s)' % x for x in got)
                            + '. Previous position %s Mb (v1.0) lay outside the span of these markers.' % (('%.3f' % clo) if clo == chi else '%.3f-%.3f' % (clo, chi))))
    return ch

if __name__ == '__main__':
    M = json.load(open('data/merged_v13f.json'))
    C = run(M['reports'])
    json.dump(C, open('/home/claude/work/r98_changes.json', 'w'), indent=0)
    for x in C: print(x['id'], x['source'], x['name'], x['chr'], x['old'], '->', x['new'], [(a, b) for a, b, h in x['placed']])
    print(len(C))
