"""Apply sb19_fix.py: replace the records of the 19 re-read papers in the master set and the page seed.
Writes data/merged_v13f.json, web/index.html (seed replaced in place), data/markerpos_all.json.gz,
and work/sb19_edits.json (old/new text pairs for patching the published page)."""
import json, re, gzip, sys
sys.path.insert(0, 'scripts/curation/sb19'); sys.path.insert(0, 'scripts/curation')
import lift
exec(open('scripts/curation/sb19_fix.py').read().split("if __name__")[0])

RM_SOURCES = set('Gahtyari2021 Plants2021SB A_spot_blotch_015 He2020Phyto Kaur2021FPS A_spot_blotch_016 Kumar2009TAG Kumar2010MB '
                 'Kumar2015MB Lillemo2013 A_spot_blotch_003 Lu2016 Z_1329 Z_1330 Roy2021FPS Singh2016MBR V_Singh2018 V_Sharma2022 '
                 'Zhang2020TAG Zhang2020Sb4 Zhu2014MB'.split())
RM_IDS = {'R0960', 'R0961', 'R0962', 'R0963', 'R0964'}   # Bemister 2019 loci taken second-hand from Liu et al. 2020, as tan spot
DROP_SOURCES = ['Plants2021SB', 'A_spot_blotch_015', 'A_spot_blotch_016', 'A_spot_blotch_003', 'Z_1330', 'V_Sharma2022', 'Zhang2020Sb4']
NEW_SOURCES = {
    'Kumar2016MB': {'citation': 'Kumar S, Röder MS, Singh RP, Kumar S, Chand R, Joshi AK, Kumar U (2016) Mapping of spot blotch disease resistance using NDVI as a substitute to visual observation in wheat (Triticum aestivum L.). Molecular Breeding 36:95', 'url': 'https://doi.org/10.1007/s11032-016-0515-6', 'short': 'Kumar et al. (2016)'},
    'Vasistha2024SR': {'citation': 'Vasistha NK, Sharma V, Singh S, Kaur R, Kumar A, Ravat VK, Kumar R, Gupta PK (2024) Meta-QTL analysis and identification of candidate genes for multiple-traits associated with spot blotch resistance in bread wheat. Scientific Reports 14:13083', 'url': 'https://doi.org/10.1038/s41598-024-63924-w', 'short': 'Vasistha et al. (2024)'},
    'Bemister2019CS': {'citation': 'Bemister DH, Semagn K, Iqbal M, Randhawa H, Strelkov SE, Spaner DM (2019) Mapping QTL associated with stripe rust, leaf rust, and leaf spotting in a Canadian spring wheat population. Crop Science 59:650-658', 'url': 'https://doi.org/10.2135/cropsci2018.05.0348', 'short': 'Bemister et al. (2019)'},
    'Iqbal2023CS': {'citation': "Iqbal M, Semagn K, Randhawa H, Aboukhaddour R, Ciechanowska I, Strenzke K, N'Diaye A, Pozniak C, Spaner D (2023) Identification and characterization of stripe rust, leaf rust, leaf spot, and common bunt resistance in spring wheat. Crop Science 63:2310-2328", 'url': 'https://doi.org/10.1002/csc2.20953', 'short': 'Iqbal et al. (2023)'},
}
FIX_SOURCES = {
    'Lillemo2013': {'citation': 'Lillemo M, Joshi AK, Prasad R, Chand R, Singh RP (2013) QTL for spot blotch resistance in bread wheat line Saar co-locate to the biotrophic disease resistance loci Lr34 and Lr46. Theoretical and Applied Genetics 126:711-719', 'url': 'https://doi.org/10.1007/s00122-012-2012-6', 'short': 'Lillemo et al. (2013)'},
    'Roy2021FPS': {'citation_fix': ('Frontiers in Plant Science 12:64', 'Frontiers in Plant Science 12:641324')},
    'He2020Phyto': {'citation_fix': (None, None)},
}
# keys whose citation was wrong get new keys (the live page lets the page seed win over database edits of an existing key)
RENAME = {'Lillemo2013': 'Lillemo2013TAG'}
SID = lambda i: 'S19-%03d' % i

def build(reports, sources, with_asm):
    keep = [r for r in reports if r['source'] not in RM_SOURCES and r['id'] not in RM_IDS]
    removed = [r for r in reports if r not in keep]
    new = []
    for i, x in enumerate(R, 1):
        r = {k: v for k, v in x.items() if k != 'conflict'}
        r['id'] = SID(i); r['curation'] = 'hand'
        r['source'] = RENAME.get(r['source'], r['source'])
        r['note'] = (r.get('note', '') + ' Re-read from the paper, Nov 2026 (sb19_fix.py).').strip()
        if with_asm:
            r['asm'] = lift.lift_record(r)
            if 'RefSeq v2.1' in r['asm'] and not str(r['ref']).startswith('RefSeq v2'):
                r['asm']['_method'] = 'anchor interpolation on atlas liftover anchors (sb19)'
        new.append(r)
    for r in keep:   # Gao et al. 2023: checked against Table 2, values unchanged
        if r['source'] == 'Gao2023MB':
            r['note'] = re.sub(r'\s*\|\s*LOD/PVE read from paper table \(auto\); verify\.?', '', r.get('note') or '')
            r['note'] += ' Checked against Gao et al. (2023) Table 2, Nov 2026.'; r['curation'] = 'hand'
    src = dict(sources)
    for k in DROP_SOURCES: src.pop(k, None)
    src.update(NEW_SOURCES)
    src['Lillemo2013TAG'] = FIX_SOURCES['Lillemo2013']
    src.pop('Lillemo2013', None)
    return keep + new, src, removed, new

if __name__ == '__main__':
    # master file
    M = json.load(open('data/merged_v13f.json'))
    M['reports'], M['sources'], rem, new = build(M['reports'], M['sources'], True)
    json.dump(M, open('data/merged_v13f.json', 'w'))
    # page seed (repo copy)
    seed = json.load(open(W + 'seed.json'))
    old_reports, old_sources = seed['reports'], dict(seed['sources'])
    seed['reports'], seed['sources'], rem2, new2 = build(seed['reports'], seed['sources'], True)
    H = open('web/index.html').read()
    m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', H, re.S)
    ser = lambda o: json.dumps(o, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/')
    H2 = H[:m.start(2)] + ser(seed) + H[m.end(2):]
    open('web/index.html', 'w').write(H2)
    # edit list for the published page: delete removed records and changed sources, insert new ones
    seedtxt = m.group(2)
    edits, missing = [], []
    for r in rem2:
        t = ser(r) + ','
        if seedtxt.count(t) == 1: edits.append({'old_text': t, 'new_text': ''})
        else: missing.append(r['id'])
    for k in DROP_SOURCES + ['Lillemo2013']:
        if k in old_sources:
            t = '"%s":%s,' % (k, ser(old_sources[k]))
            nt = ''
            if seedtxt.count(t) == 1: edits.append({'old_text': t, 'new_text': nt})
            else: missing.append('source ' + k)
    json.dump({'edits': edits, 'missing': missing, 'new_reports': [ser(r) for r in new2],
               'new_sources': {k: ser(v) for k, v in NEW_SOURCES.items()},
               'removed': [(r['id'], r['source'], r['name']) for r in rem2]}, open(W + 'sb19_edits.json', 'w'), ensure_ascii=False)
    # marker dictionary: wrong 2B entry for BS00036421
    MK['bs00036421'] = [['4D', 32.3474, 'RefSeq v1.0', 'Sun et al. 2020 PBJ (90K), as BS00036421_51; corrected Nov 2026 (was 2B 1.0666)']]
    json.dump(MK, gzip.open('data/markerpos_all.json.gz', 'wt'))
    print('master: removed', len(rem), 'added', len(new), 'total', len(M['reports']))
    print('seed: removed', len(rem2), 'added', len(new2), 'total', len(seed['reports']))
    print('edits', len(edits), 'missing', missing, 'old_text KiB', round(sum(len(e['old_text'].encode()) for e in edits) / 1024, 1))
