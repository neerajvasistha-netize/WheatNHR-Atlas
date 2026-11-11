"""Batch 5 (Nov 2026): tan spot papers from Archive_1.zip (plus Liu et al. 2006 SNB). Records are defined in ts52_recs.py and placed
with the sb26 resolver through the stb44 marker wrapper. Writes data/merged_v13f.json, web/index.html (seed) and /home/claude/work/ts52_db.json.
Run after stb44.py."""
import json, re, sys, io, contextlib
sys.path.insert(0, 'scripts/curation'); sys.path.insert(0, 'scripts/curation/sb19')
with contextlib.redirect_stdout(io.StringIO()):
    import sb26_fix as F
import lift, ts52_recs as RR

RM_SOURCES = ('A_Tan_spot_024 A_Tan_spot_025 Z_1097 A_Tan_spot_022 Z_1110 V_Kalia2018 Z_1198 V_Kariyawasam2018 V_Li2010b V_Liu2017 '
              'V_Shankar2017 V_Zou2017').split()
KEYMAP = {}
def S(cit, url, short): return dict(citation=cit, url=url, short=short)
NEW_SOURCES = {
 'Abeysekara2010PG': S('Abeysekara NS, Friesen TL, Liu ZH, McClean PE, Faris JD (2010) Marker development and saturation mapping of the tan spot Ptr ToxB sensitivity locus Tsc2 in hexaploid wheat. The Plant Genome 3:179-189', 'https://doi.org/10.3835/plantgenome2010.07.0017', 'Abeysekara et al. (2010)'),
 'Chu2010MBts': S('Chu CG, Chao S, Friesen TL, Faris JD, Zhong S, Xu SS (2010) Identification of novel tan spot resistance QTLs using an SSR-based linkage map of tetraploid wheat. Molecular Breeding 25:327-338', 'https://doi.org/10.1007/s11032-009-9335-2', 'Chu et al. (2010)'),
 'Cheong2004AJAR': S('Cheong J, Wallwork H, Williams KJ (2004) Identification of a major QTL for yellow leaf spot resistance in the wheat varieties Brookton and Cranbrook. Australian Journal of Agricultural Research 55:315-319', 'https://doi.org/10.1071/AR03140', 'Cheong et al. (2004)'),
 'Corsi2020TAG': S('Corsi B, Percival-Alwyn L, Downie RC, Venturini L, Iagallo EM, Mantello CC, McCormick-Barnes C, See PT, Oliver RP, Moffat CS, Cockram J (2020) Genetic analysis of wheat sensitivity to the ToxB fungal effector from Pyrenophora tritici-repentis, the causal agent of tan spot. Theoretical and Applied Genetics 133:935-950', 'https://doi.org/10.1007/s00122-019-03517-8', 'Corsi et al. (2020)'),
 'Effertz2002Phyto': S('Effertz RJ, Meinhardt SW, Anderson JA, Jordahl JG, Francl LJ (2002) Identification of a chlorosis-inducing toxin from Pyrenophora tritici-repentis and the chromosomal location of an insensitivity locus in wheat. Phytopathology 92:527-533', 'https://doi.org/10.1094/PHYTO.2002.92.5.527', 'Effertz et al. (2002)'),
 'Faris2020TAG': S('Faris JD, Overlander ME, Kariyawasam GK, Carter A, Xu SS, Liu Z (2020) Identification of a major dominant gene for race-nonspecific tan spot resistance in wild emmer wheat. Theoretical and Applied Genetics 133:829-841', 'https://doi.org/10.1007/s00122-019-03509-8', 'Faris et al. (2020)'),
 'Li2011TAG': S('Li HB, Yan W, Liu GR, Wen SM, Liu CJ (2011) Identification and validation of quantitative trait loci conferring tan spot resistance in the bread wheat variety Ernie. Theoretical and Applied Genetics 122:395-403', 'https://doi.org/10.1007/s00122-010-1455-x', 'Li et al. (2011)'),
 'Liu2006Genome': S('Liu Z, Friesen TL, Ling H, Meinhardt SW, Oliver RP, Rasmussen JB, Faris JD (2006) The Tsn1-ToxA interaction in the wheat-Stagonospora nodorum pathosystem parallels that of the wheat-tan spot system. Genome 49:1265-1273', 'https://doi.org/10.1139/g06-088', 'Liu et al. (2006)'),
 'Liu2020TAGts': S('Liu Y, Zhang Q, Salsman E, Fiedler JD, Hegstad JB, Liu Z, Faris JD, Xu SS, Li X (2020) QTL mapping of resistance to tan spot induced by race 2 of Pyrenophora tritici-repentis in tetraploid wheat. Theoretical and Applied Genetics 133:433-442', 'https://doi.org/10.1007/s00122-019-03474-2', 'Liu et al. (2020)'),
 'Singh2008MB': S('Singh PK, Mergoum M, Gonzalez-Hernandez JL, Ali S, Adhikari TB, Kianian SF, Elias EM, Hughes GR (2008) Genetics and molecular mapping of resistance to necrosis inducing race 5 of Pyrenophora tritici-repentis in tetraploid wheat. Molecular Breeding 21:293-304', 'https://doi.org/10.1007/s11032-007-9129-3', 'Singh et al. (2008)'),
 'Sun2010Phyto': S('Sun XC, Bockus W, Bai GH (2010) Quantitative trait loci for resistance to Pyrenophora tritici-repentis race 1 in a Chinese wheat. Phytopathology 100:468-473', 'https://doi.org/10.1094/PHYTO-100-5-0468', 'Sun et al. (2010)'),
}

POS0 = F.pos
def P(m, c):
    if m.lower().lstrip('x') in ('barc008', 'barc8') and c == '1B':   # cosegregates with Xwmc230 on the Dreisigacker et al. 2015 map
        p, h = POS0('Xwmc230', '1B')
        return p, 'position of Xwmc230, with which Xbarc008 cosegregates (Dreisigacker et al. 2015 map)'
    for v in (m, ('wsnp_' + m) if m.startswith(('Ku_', 'Ex_', 'Ra_')) else None, m.replace('Tag_', ''), m.replace('wPt_', 'wPt-')):
        if not v: continue
        p, h = POS0(v, c)
        if p is not None: return p, h
    return None, ''

def make(x):
    rec = dict(name=x['name'], chr=x['chr'], source=KEYMAP.get(x['source'], x['source']), pve=x['pve'], stat=x['stat'], note=x['note'],
               evidence=x['evidence'], ref='RefSeq v1.0', curation='hand', cls=x['extra'].get('cls', 'QTL'),
               disease=x['extra']['disease'], marker=', '.join(dict.fromkeys(x['markers'] + ([x['peak_marker']] if x['peak_marker'] else []))),
               start=None, end=None, peak=None, grade=None)
    fx = x['fixed']
    if fx:
        for k in ('start', 'end', 'peak', 'grade'):
            if fx.get(k) is not None: rec[k] = fx[k]
        if x['markers']:
            got = [(m,) + P(m, x['chr']) for m in x['markers']]
            got = [g for g in got if g[1] is not None]
            if got: rec['note'] += ' Marker positions (RefSeq v1.0): ' + '; '.join('%s %.3f Mb (%s)' % g for g in got) + '.'
    else:
        F.pos = P
        try: F.place(rec, x['markers'], x['chr'], peak_marker=x['peak_marker'])
        finally: F.pos = POS0
    for k in ('start', 'end', 'peak'):
        if rec[k] is not None: rec[k] = round(float(rec[k]), 4)
    if rec['peak'] is not None and rec['start'] is not None:
        rec['start'], rec['end'] = min(rec['start'], rec['peak']), max(rec['end'], rec['peak'])
    rec['grade'] = rec['grade'] or 'D'
    return rec

R = [make(x) for x in RR.R]
SID = lambda i: 'S52-%03d' % i

PAIRS = {}
def wide_fix(r):
    """Fallback past the last liftover anchor: 90K markers with both a RefSeq v1.0 (Sun 2020) and a v2.1 (Hu 2022) position
    (as sb26_apply.pair_lift). Interval: v2.1 extent of markers inside it; point: offset from the nearest marker."""
    import gzip
    c = r['chr']
    if c not in PAIRS:
        M = json.load(gzip.open('data/markerpos_all.json.gz'))
        L = []
        for val in M.values():
            a = [x[1] for x in val if x[0] == c and x[2] == 'RefSeq v1.0']
            b = [x[1] for x in val if x[0] == c and x[2] == 'RefSeq v2.1']
            if len(a) == 1 and len(b) == 1: L.append((a[0], b[0]))
        PAIRS[c] = L
    L = PAIRS[c]
    if r.get('start') is not None:
        ys = [b for a, b in L if r['start'] - 0.05 <= a <= r['end'] + 0.05]
        assert len(ys) >= 5, (r['name'], len(ys))
        return {'start': round(min(ys), 3), 'end': round(max(ys), 3), 'flag': 'rechecked', 'method': '90K marker pairs (v1.0 Sun 2020, v2.1 Hu 2022), n=%d' % len(ys)}
    a, b = min(L, key=lambda t: abs(t[0] - r['peak']))
    return {'peak': round(b + r['peak'] - a, 3), 'flag': 'rechecked', 'method': 'offset from the nearest 90K marker pair (v1.0 %.3f, v2.1 %.3f)' % (a, b)}

TRIM = {'S52-095': (' Flanking markers lie 291 Mb apart; interval is wide.',),   # already said in the record note
        'S52-115': (' Only Xgwm443 could be placed on 5A; position is that marker, not a reported peak.',)}

def build(reports, sources):
    keep = [r for r in reports if r['source'] not in RM_SOURCES and not r['id'].startswith('S52-')]
    removed = [r['id'] for r in reports if r['source'] in RM_SOURCES]
    new = []
    for i, x in enumerate(R, 1):
        r = dict(x); r['id'] = SID(i)
        r['note'] = (r['note'] + ' Read from the paper, Nov 2026 (ts52.py).').strip()
        for a in TRIM.get(r['id'], ()): r['note'] = r['note'].replace(a, '')
        r['asm'] = lift.lift_record(r)
        v = r['asm'].get('RefSeq v2.1')
        if v and v.get('flag') == 'wide_gap':
            v.update(wide_fix(r))
        new.append(r)
    src = dict(sources); src.update(NEW_SOURCES)
    return keep + new, src, new, removed

if __name__ == '__main__':
    M = json.load(open('data/merged_v13f.json'))
    M['reports'], M['sources'], new, rem = build(M['reports'], M['sources'])
    json.dump(M, open('data/merged_v13f.json', 'w'))
    H = open('web/index.html').read()
    m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', H, re.S)
    S = json.loads(m.group(2).replace('<\\/', '</'))
    S['reports'], S['sources'], new2, rem2 = build(S['reports'], S['sources'])
    ser = lambda o: json.dumps(o, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/')
    open('web/index.html', 'w').write(H[:m.start(2)] + ser(S) + H[m.end(2):])
    docs = []
    for r in new2:
        d = {k: v for k, v in r.items() if v not in (None, '') and k != 'asm'}
        d['asm'] = {k: r['asm'][k] for k in ('RefSeq v2.1', 'RefSeq v1.0', '_note') if k in r['asm']}
        docs.append(d)
    json.dump({'docs': docs, 'removed': rem2, 'sources': NEW_SOURCES}, open('/home/claude/work/ts52_db.json', 'w'), ensure_ascii=False)
    import collections
    print('new', len(new), 'removed', len(rem), len(rem2), collections.Counter(r['grade'] for r in new),
          collections.Counter((r['asm'].get('RefSeq v2.1') or {}).get('flag', 'none') for r in new))
    for r in new: print(r['id'], r['grade'], r['chr'], r.get('start'), r.get('end'), r.get('peak'), r['name'][:40])
