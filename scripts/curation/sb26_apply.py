"""Apply sb26_fix.py: replace the records of the 26 re-read SNB / tan spot papers (batch 2, 11 Nov 2026).
Writes data/merged_v13f.json, web/index.html (seed), and work/sb26_db.json (documents for the live database)."""
import json, re, sys
sys.path.insert(0, 'scripts/curation/sb19'); sys.path.insert(0, 'scripts/curation')
import lift
exec(open('scripts/curation/sb26_fix.py').read().split("if __name__")[0])

RM_SOURCES = set('V_RiveraBurgos2022 Z_1239 A_SNB_019 Z_1440 A_SNB_007 A_SNB_010 Z_1021 Z_1366 Arseniuk2004 '
                 'V_Abeysekara2012 V_Hu2019 V_Francki2018 Z_1039 A_SNB_022 Y_166 Z_1090 Z_1360'.split())
NEW_SOURCES = {
    'Liu2004Phyto': {'citation': 'Liu ZH, Friesen TL, Rasmussen JB, Ali S, Meinhardt SW, Faris JD (2004) Quantitative trait loci analysis and mapping of seedling resistance to Stagonospora nodorum leaf blotch in wheat. Phytopathology 94:1061-1067', 'url': 'https://doi.org/10.1094/PHYTO.2004.94.10.1061', 'short': 'Liu et al. (2004)'},
    'GonzalezHernandez2009': {'citation': 'Gonzalez-Hernandez JL, Singh PK, Mergoum M, Adhikari TB, Kianian SF, Simsek S, Elias EM (2009) A quantitative trait locus on chromosome 5B controls resistance of Triticum turgidum (L.) var. diccocoides to Stagonospora nodorum blotch. Euphytica 166:199-206', 'url': 'https://doi.org/10.1007/s10681-008-9825-z', 'short': 'Gonzalez-Hernandez et al. (2009)'},
    'Reszka2007PPB': {'citation': 'Reszka E, Song Q, Arseniuk E, Cregan PB, Ueng PP (2007) The QTL controlling partial resistance to Stagonospora nodorum blotch disease in winter triticale Bogo. Plant Pathology Bulletin 16:161-167', 'url': '', 'short': 'Reszka et al. (2007)'},
    'Zhang2009MB': {'citation': 'Zhang Z, Friesen TL, Simons KJ, Xu SS, Faris JD (2009) Development, identification, and validation of markers for marker-assisted selection against the Stagonospora nodorum toxin sensitivity genes Tsn1 and Snn2 in wheat. Molecular Breeding 23:35-49', 'url': 'https://doi.org/10.1007/s11032-008-9211-5', 'short': 'Zhang et al. (2009)'},
    'Lin2021TAG': {'citation': 'Lin M, Stadlmeier M, Mohler V, Tan KC, Ficke A, Cockram J, Lillemo M (2021) Identification and cross-validation of genetic loci conferring resistance to Septoria nodorum blotch using a German multi-founder winter wheat population. Theoretical and Applied Genetics 134:125-142', 'url': 'https://doi.org/10.1007/s00122-020-03686-x', 'short': 'Lin et al. (2021)'},
    'Shi2016MGG': {'citation': 'Shi G, Zhang Z, Friesen TL, Bansal U, Cloutier S, Wicker T, Rasmussen JB, Faris JD (2016) Marker development, saturation mapping, and high-resolution mapping of the Septoria nodorum blotch susceptibility gene Snn3-B1 in wheat. Molecular Genetics and Genomics 291:107-119', 'url': 'https://doi.org/10.1007/s00438-015-1091-x', 'short': 'Shi et al. (2016)'},
}
SID = lambda i: 'S26-%03d' % i

def pair_lift(c, s, e):
    """Fallback for intervals beyond the last record anchor: use 90K markers that carry both a
    RefSeq v1.0 (Sun 2020) and a RefSeq v2.1 (Hu 2022) position, and take the v2.1 extent of every
    marker inside [s, e]. Handles terminal inversions such as distal 2AL."""
    import gzip
    M = json.load(gzip.open('data/markerpos_all.json.gz'))
    ys = []
    for v in M.values():
        a = [x[1] for x in v if x[0] == c and x[2] == 'RefSeq v1.0']
        b = [x[1] for x in v if x[0] == c and x[2] == 'RefSeq v2.1']
        if len(a) == 1 and len(b) == 1 and s - 0.05 <= a[0] <= e + 0.05: ys.append(b[0])
    assert len(ys) >= 5, (c, s, e, len(ys))
    return {'start': round(min(ys), 3), 'end': round(max(ys), 3), 'flag': 'rechecked',
            'method': '90K marker pairs (v1.0 Sun 2020, v2.1 Hu 2022), n=%d' % len(ys)}

def build(reports, sources):
    keep = [r for r in reports if r['source'] not in RM_SOURCES]
    removed = [r for r in reports if r['source'] in RM_SOURCES]
    new = []
    for i, x in enumerate(R, 1):
        r = dict(x); r['id'] = SID(i)
        r['note'] = (r.get('note', '') + ' Re-read from the paper, Nov 2026 (sb26_fix.py).').strip()
        r['asm'] = lift.lift_record(r)
        v = r['asm'].get('RefSeq v2.1')
        if v and v.get('flag') == 'wide_gap':
            v.update(pair_lift(r['chr'], float(r['start']), float(r['end'])))
        new.append(r)
    src = dict(sources); src.update(NEW_SOURCES)
    return keep + new, src, removed, new

if __name__ == '__main__':
    M = json.load(open('data/merged_v13f.json'))
    M['reports'], M['sources'], rem, new = build(M['reports'], M['sources'])
    json.dump(M, open('data/merged_v13f.json', 'w'))
    H = open('web/index.html').read()
    m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', H, re.S)
    seed = json.loads(m.group(2).replace('<\\/', '</'))
    seed['reports'], seed['sources'], rem2, new2 = build(seed['reports'], seed['sources'])
    ser = lambda o: json.dumps(o, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/')
    open('web/index.html', 'w').write(H[:m.start(2)] + ser(seed) + H[m.end(2):])
    docs = []
    for r in new2:
        d = {k: v for k, v in r.items() if v not in (None, '') and k != 'asm'}
        d['asm'] = {k: r['asm'][k] for k in ('RefSeq v2.1', 'RefSeq v1.0', '_note') if k in r['asm']}
        docs.append(d)
    json.dump({'docs': docs, 'removed': [r['id'] for r in rem2], 'sources': NEW_SOURCES}, open('/home/claude/work/sb26_db.json', 'w'), ensure_ascii=False)
    import collections
    print('master removed', len(rem), 'added', len(new), '| seed removed', len(rem2), 'added', len(new2))
    print('v2.1 flags', collections.Counter((r['asm'].get('RefSeq v2.1') or {}).get('flag', 'none') for r in new2))
