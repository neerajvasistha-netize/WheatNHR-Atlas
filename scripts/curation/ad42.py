"""ad42 (Nov 2026): Herter et al. 2018 TAG 131:1263 FHB QTL placed by decoding the marker codes in ESM2 (Supplementary Table 2);
Capo x Arina 1B and 2A QTL refined with the in-silico positions of the SSRs on the first map of Online Resource 3 (Xgwm11, Xgwm558, Xgwm312, Xgwm294).
Dinglasan et al. 2019 Online Resources 2-3 checked: they give DArTseq IDs with consensus cM only (no sequences), adding no placeable anchors."""
import json, re, sys
exec(open('scripts/curation/lib/insilico_res.py').read())
sys.path.insert(0, 'scripts/curation/sb19'); import lift
KEY = json.load(open('docs/supplements/herter2018_marker_key.json'))
def fix(rs):
    for r in rs:
        if r['id'] in ('S133-001', 'S133-002', 'S133-003'):
            a, b = re.findall(r'M\d+', r['marker'])
            pl = [(c, KEY[c]) + tuple(resolve(KEY[c], r['chr'])) for c in (a, b)]
            r['marker'] = f'{KEY[a]} ({a}), {KEY[b]} ({b})'
            ok = [x for x in pl if x[2]]
            base = r['note'].split(' Flanking markers are given only as codes')[0]
            rest = r['note'].split('cannot be placed.')[1].split(' Read from the paper')[0]
            if len(ok) == 2:
                r['start'], r['end'] = min(x[2] for x in ok), max(x[2] for x in ok); r['peak'] = None; r['grade'] = 'B'
                t = ' Placed in silico: ' + ', '.join(f'{n} {p:.3f} Mb ({h})' for c, n, p, h in ok) + '.'
            else:
                c, n, p, h = ok[0]; r['peak'], r['grade'] = p, 'C'
                t = f' Only {n} could be placed ({p:.3f} Mb, {h}); position is that marker.'
            r['note'] = base + ' Flanking marker codes decoded from Supplementary Table 2 (ESM2).' + rest + t + ' Read from the paper and ESM2, Nov 2026 (ad42.py).'
            r['evidence'] += '; ESM2 marker key'
        if r['id'] == 'S195-001':
            r['marker'] = 'wPt-7721, wPt-8490, wPt-0568, wPt-8068 (QTL); wPt-1657, Xgwm558, Xgwm312, Xgwm294 (map anchors)'
            r['start'], r['end'], r['peak'], r['grade'] = 0.0, 47.1528, None, 'C'
            r['note'] = r['note'].split(' Only wPt-7721 could be placed')[0] + (' The SSRs on the same map orient it: Xgwm558 (about 33 cM) places at 361.2 Mb and Xgwm312/Xgwm294 (85 cM) at 709.0/716.5 Mb (in-silico PCR), so 0 cM is the distal end of 2AS. '
                'The QTL is therefore placed from the 2AS end to wPt-1657 (about 20 cM, 47.153 Mb, DArT tag BLAST); wPt-7721 (0 cM) places at 29.941 Mb. A 2012-only LOD peak near 40-45 cM is not entered. Read from Online Resource 3 (Fig. S1), Nov 2026 (ad40.py, ad42.py).')
        if r['id'] == 'V1433':
            r['marker'] = 'wPt-3103, wPt-6117, wPt-666564; Xgwm11 (map anchor)'
            r['note'] = r['note'].replace(' Read from Online Resource 3 (Fig. S1), Nov 2026 (ad40.py).', ' Xgwm11, just proximal of the LOD support interval (about 29 cM), places at 216.036 Mb (in-silico PCR), consistent with a 1BS QTL distal of it. Read from Online Resource 3 (Fig. S1), Nov 2026 (ad40.py, ad42.py).')
        if r['id'] in ('S133-001', 'S133-002', 'S133-003', 'S195-001', 'V1433'): r['asm'] = lift.lift_record(r)
    return rs
IDS = ['S133-001', 'S133-002', 'S133-003', 'S195-001', 'V1433']
if __name__ == '__main__':
    M = json.load(open('data/merged_v13f.json')); fix(M['reports']); json.dump(M, open('data/merged_v13f.json', 'w'))
    H = open('web/index.html').read(); m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', H, re.S)
    Sd = json.loads(m.group(2).replace('<\\/', '</')); fix(Sd['reports'])
    open('web/index.html', 'w').write(H[:m.start(2)] + json.dumps(Sd, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/') + H[m.end(2):])
    ops = []
    for r in M['reports']:
        if r['id'] in IDS:
            d = {k: v for k, v in r.items() if v not in (None, '') and k != 'asm'}
            if r['id'] == 'S195-001': d['start'] = 0.0
            d['asm'] = {k: r['asm'][k] for k in ('RefSeq v2.1', '_note') if k in r['asm']}
            ops.append({'op': 'set', 'path': 'reports/' + r['id'], 'ifVersion': 1, 'data': d})
            print(r['id'], r['grade'], r.get('start'), r.get('end'), r.get('peak'), d['asm'])
    json.dump(ops, open('/tmp/ops42.json', 'w'), ensure_ascii=False)
